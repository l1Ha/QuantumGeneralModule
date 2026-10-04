"""
Sum-of-products (SOP) Hamiltonians and the POTFIT decomposition.

Python reference implementation of book section 18.1 and the proposed
`mod_sop_hamiltonian` (book section 18.8):

    H = sum_alpha c_alpha prod_kappa h_alpha^(kappa)

The matrix-free action costs O(M D n^(D+1)) for dense local factors and never
builds the full N_tot = prod n_kappa matrix. Potential energy surfaces enter as
rank-1 (vector) factors, which are applied elementwise.

Fortran counterpart (roadmap): `mod_sop_hamiltonian` with
`sop_hamiltonian_t`, `sop_apply`, `sop_from_potfit`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np

__all__ = ["SOPTerm", "SOPHamiltonian", "sop_apply", "sop_to_dense",
           "sop_from_potfit"]


@dataclass
class SOPTerm:
    """One product term: c * h^(1) x h^(2) x ... x h^(D).

    Each factor is either
      * an (n_k, n_k) matrix acting on mode k, or
      * an (n_k,) vector representing a diagonal (potential) factor.
    """
    coeff: complex
    factors: List[np.ndarray]


@dataclass
class SOPHamiltonian:
    """H = sum_alpha c_alpha prod_kappa h_alpha^(kappa)."""
    dims: Sequence[int]
    terms: List[SOPTerm] = field(default_factory=list)

    def __post_init__(self):
        self.dims = tuple(int(d) for d in self.dims)
        for term in self.terms:
            if len(term.factors) != len(self.dims):
                raise ValueError("each term needs one factor per mode")

    def apply(self, v: np.ndarray) -> np.ndarray:
        """Matrix-free action v -> H v (book section 18.1).

        Cost: O(M D n^(D+1)) for dense factors; diagonal (vector) factors are
        applied elementwise in O(M n^D).
        """
        v = np.asarray(v, dtype=complex)
        if v.shape != self.dims:
            raise ValueError(f"v has shape {v.shape}, expected {self.dims}")
        out = np.zeros_like(v)
        for term in self.terms:
            tmp = v
            # contract from the last axis backward so the leading axes stay contiguous
            for k in range(len(self.dims) - 1, -1, -1):
                fac = term.factors[k]
                if fac.ndim == 1:
                    shape = [1] * len(self.dims)
                    shape[k] = -1
                    tmp = tmp * fac.reshape(shape)
                else:
                    tmp = np.moveaxis(np.tensordot(fac, tmp, axes=([1], [k])), 0, k)
            out += term.coeff * tmp
        return out

    def to_dense(self) -> np.ndarray:
        """Dense (N, N) matrix; only for testing on small systems."""
        n_tot = int(np.prod(self.dims))
        eye = np.eye(n_tot, dtype=complex).reshape(self.dims + (n_tot,))
        cols = [self.apply(eye[..., i]) for i in range(n_tot)]
        return np.stack(cols, axis=-1).reshape(n_tot, n_tot)

    def to_tt_operator(self):
        """Represent the SOP as a TT-operator (see tensor_train.sop_to_tt_operator)."""
        from .tensor_train import sop_to_tt_operator
        return sop_to_tt_operator(self)


def sop_apply(sop: SOPHamiltonian, v: np.ndarray) -> np.ndarray:
    """Functional form of the matrix-free SOP action."""
    return sop.apply(v)


def sop_to_dense(sop: SOPHamiltonian) -> np.ndarray:
    """Functional form of the dense reconstruction (testing only)."""
    return sop.to_dense()


def sop_from_potfit(v_grid: np.ndarray, eps: float = 1.0e-8,
                    max_terms: Optional[int] = None) -> SOPHamiltonian:
    """POTFIT-type decomposition of a full-grid potential into an SOP.

    The tensor is compressed by successive SVDs (the HOSVD/TT-SVD peeling of
    book section 18.3) and the resulting TT cores are expanded back into
    explicit product terms, so the result is a genuine SOP of 1D vector
    factors with absolute error ~ eps * ||V||_F.

    Parameters
    ----------
    v_grid : (n_1, ..., n_D) potential values on the product grid.
    eps    : target relative error of the decomposition.
    max_terms : optional cap on the number of product terms.
    """
    from .tensor_train import tt_from_dense, tt_to_sop
    v_grid = np.asarray(v_grid, dtype=float)
    tt = tt_from_dense(v_grid, eps=eps, max_rank=max_terms)
    return tt_to_sop(tt)
