"""
Tensor train (TT) / matrix product state (MPS) representation.

Python reference implementation of book section 18.3 and the proposed
`mod_tensor_train` (book section 18.8): TT-SVD, rounding, inner product, and
TT-operator application. Storage is O(sum n_k r_{k-1} r_k) instead of O(n^D).

TT-Cross, DMRG-type eigensolvers, TDVP and TEBD (book section 18.3, items
3-6) remain on the roadmap; the operations here cover the representation
layer of the verification hierarchy in book section 18.9.

Fortran counterpart (roadmap): `mod_tensor_train` with `tt_tensor_t`,
`tt_from_dense`, `tt_round`, `tt_apply`, `tt_dot`.

References:
    I. V. Oseledets, "Tensor-train decomposition", SIAM J. Sci. Comput. 33,
    2295 (2011). DOI: 10.1137/090752286.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

__all__ = ["TTTensor", "tt_from_dense", "tt_to_dense", "tt_round", "tt_dot",
           "tt_norm", "tt_zeros", "tt_to_sop", "sop_to_tt_operator",
           "tt_operator_apply", "tt_operator_to_dense"]


@dataclass
class TTTensor:
    """Tensor train with cores G^(k) of shape (r_{k-1}, n_k, r_k)."""
    cores: List[np.ndarray]

    def __post_init__(self):
        if len(self.cores) == 0:
            raise ValueError("TTTensor needs at least one core")
        if self.cores[0].shape[0] != 1 or self.cores[-1].shape[2] != 1:
            raise ValueError("boundary ranks must be 1")

    @property
    def dims(self) -> tuple:
        return tuple(c.shape[1] for c in self.cores)

    @property
    def ranks(self) -> tuple:
        return tuple(c.shape[2] for c in self.cores[:-1])


def _truncation_rank(s: np.ndarray, delta: float) -> int:
    """Smallest rank with sum s_i^2 >= (1 - delta^2) sum s_i^2 (at least 1)."""
    total = float(np.sum(s ** 2))
    if total <= 0.0:
        return 1
    kept = int(np.searchsorted(np.cumsum(s ** 2), (1.0 - delta ** 2) * total)) + 1
    return max(1, min(kept, s.size))


def tt_from_dense(a: np.ndarray, eps: float = 1.0e-10,
                  max_rank: Optional[int] = None) -> TTTensor:
    """TT-SVD: compress a dense tensor to relative error ~ eps (book 18.3).

    The per-step threshold eps/sqrt(D-1) guarantees the total error bound
    ||A - TT(A)||_F <= eps ||A||_F (Oseledets 2011, Theorem 2.4).
    """
    a = np.asarray(a, dtype=complex)
    if a.ndim == 0:
        return TTTensor([a.reshape(1, 1, 1)])
    if a.ndim == 1:
        return TTTensor([a.reshape(1, -1, 1)])
    d = a.ndim
    norm_a = float(np.linalg.norm(a))
    if norm_a == 0.0:
        return tt_zeros(a.shape)
    delta = eps / max(np.sqrt(d - 1), 1.0)
    c = a
    r_prev = 1
    cores: List[np.ndarray] = []
    for k in range(d - 1):
        m = c.reshape(r_prev * a.shape[k], -1)
        u, s, vh = np.linalg.svd(m, full_matrices=False)
        rank = _truncation_rank(s, delta)
        if max_rank is not None:
            rank = min(rank, max_rank)
        cores.append(u[:, :rank].reshape(r_prev, a.shape[k], rank))
        c = s[:rank, None] * vh[:rank]
        r_prev = rank
    cores.append(c.reshape(r_prev, a.shape[-1], 1))
    return TTTensor(cores)


def tt_to_dense(tt: TTTensor) -> np.ndarray:
    """Reconstruct the dense tensor from TT cores."""
    res = tt.cores[0][0]
    for core in tt.cores[1:]:
        res = np.tensordot(res, core, axes=([-1], [0]))
    return res[..., 0]


def tt_round(tt: TTTensor, eps: float = 1.0e-10,
             max_rank: Optional[int] = None) -> TTTensor:
    """TT rounding: re-orthogonalize, then truncate to relative error ~ eps.

    Left-to-right QR orthogonalizes the cores; right-to-left SVDs truncate
    the ranks while minimizing the Frobenius error (book 18.3, TT rounding).
    """
    cores = [c.copy() for c in tt.cores]
    d = len(cores)
    for k in range(d - 1):
        r_prev, n, r = cores[k].shape
        q, rr = np.linalg.qr(cores[k].reshape(r_prev * n, r))
        cores[k] = q.reshape(r_prev, n, q.shape[1])
        cores[k + 1] = np.tensordot(rr, cores[k + 1], axes=([1], [0]))
    delta = eps / max(np.sqrt(d - 1), 1.0)
    for k in range(d - 1, 0, -1):
        r_prev, n, r = cores[k].shape
        u, s, vh = np.linalg.svd(cores[k].reshape(r_prev, n * r), full_matrices=False)
        rank = _truncation_rank(s, delta)
        if max_rank is not None:
            rank = min(rank, max_rank)
        # Vh rows (length n*r) become the right-orthogonal core; U S stays on
        # the bond and is folded into the previous core
        cores[k] = vh[:rank].reshape(rank, n, r)
        us = u[:, :rank] * s[:rank]                     # (r_prev, rank)
        cores[k - 1] = np.tensordot(cores[k - 1], us, axes=([2], [0]))
    return TTTensor(cores)


def tt_dot(a: TTTensor, b: TTTensor) -> complex:
    """Inner product <a, b> with conjugation on b (book 18.3)."""
    if a.dims != b.dims:
        raise ValueError("dimension mismatch")
    res = np.ones((1, 1), dtype=complex)
    for ca, cb in zip(a.cores, b.cores):
        res = np.einsum("ab,aic,bid->cd", res, ca, cb.conj())
    return complex(res[0, 0])


def tt_norm(a: TTTensor) -> float:
    """Frobenius norm of a TT tensor."""
    return float(np.sqrt(abs(tt_dot(a, a)).real))


def tt_zeros(shape: Sequence[int]) -> TTTensor:
    """Zero tensor in TT format with all ranks 1."""
    return TTTensor([np.zeros((1, int(n), 1), dtype=complex) for n in shape])


def tt_to_sop(tt: TTTensor):
    """Expand a TT tensor into explicit SOP product terms (POTFIT back-end).

    C_i = sum_alpha prod_k G^(k)[alpha_{k-1}, :, alpha_k], so the expansion
    has prod(r_k) terms, each a rank-1 product of 1D vectors.
    """
    from .sop_hamiltonian import SOPHamiltonian, SOPTerm
    dims = tt.dims
    terms: List = []
    if len(dims) == 1:
        terms.append(SOPTerm(coeff=1.0 + 0.0j, factors=[tt.cores[0][0, :, 0]]))
        return SOPHamiltonian(dims=dims, terms=terms)
    import itertools
    for combo in itertools.product(*(range(c.shape[2]) for c in tt.cores[:-1])):
        factors = []
        idx_prev = 0
        for k, core in enumerate(tt.cores):
            idx_next = combo[k] if k < len(combo) else 0
            factors.append(np.asarray(core[idx_prev, :, idx_next]))
            idx_prev = idx_next
        terms.append(SOPTerm(coeff=1.0 + 0.0j, factors=factors))
    return SOPHamiltonian(dims=dims, terms=terms)


def sop_to_tt_operator(sop):
    """Build a TT-operator (4-leg cores) from a SOP Hamiltonian.

    Each SOP term is a rank-1 operator TT; their sum has TT ranks equal to
    the number of terms. Operator cores have shape (r_{k-1}, n_out, n_in, r_k).
    """
    d = len(sop.dims)
    m = len(sop.terms)
    if m == 0:
        return [np.zeros((1, n, n, 1), dtype=complex) for n in sop.dims]

    def factor_matrix(term, k):
        fac = term.factors[k]
        h = np.diag(np.asarray(fac, dtype=complex)) if fac.ndim == 1 else np.asarray(fac, dtype=complex)
        if h.shape != (sop.dims[k], sop.dims[k]):
            raise ValueError("SOP factor dimension mismatch")
        # the term coefficient enters exactly once (in the first core)
        return term.coeff * h if k == 0 else h

    if d == 1:
        return [(sum(factor_matrix(t, 0) for t in sop.terms)).reshape(1, sop.dims[0], sop.dims[0], 1)]
    cores: List[np.ndarray] = []
    # first core: (1, n, n, M) — the branch index selects the term
    first = np.stack([factor_matrix(t, 0) for t in sop.terms], axis=2)[None]
    cores.append(first)
    # interior cores: diagonal blocks h^(k)_alpha
    for k in range(1, d - 1):
        core = np.zeros((m, sop.dims[k], sop.dims[k], m), dtype=complex)
        for a, t in enumerate(sop.terms):
            core[a, :, :, a] = factor_matrix(t, k)
        cores.append(core)
    # last core: (M, n, n, 1)
    last = np.stack([factor_matrix(t, d - 1) for t in sop.terms], axis=0)[..., None]
    cores.append(last)
    return cores


def tt_operator_apply(op_cores: List[np.ndarray], v: TTTensor,
                      eps: float = 1.0e-10, max_rank: Optional[int] = None) -> TTTensor:
    """Apply a TT-operator to a TT vector: w = H v with rounding (book 18.3).

    The exact product has ranks r^H_k * r^v_k; the final TT-SVD + rounding
    controls the rank growth. (On-the-fly rounding between contractions, as
    used in production TT solvers, is on the roadmap; this implementation
    targets the verification hierarchy of book section 18.9.)
    """
    res = np.ones((1, 1), dtype=complex)
    for ocore, vcore in zip(op_cores, v.cores):
        # res (..., r_op, r_vec); ocore (r_op, n_out, n_in, r_op'); vcore (r_vec, n_in, r_vec')
        res = np.einsum("...ab,axyc,byd->...xcd", res, ocore, vcore)
    dense = res[..., 0, 0]
    return tt_round(tt_from_dense(dense, eps=eps, max_rank=max_rank),
                    eps=eps, max_rank=max_rank)


def tt_operator_to_dense(op_cores: List[np.ndarray]) -> np.ndarray:
    """Dense matrix of a TT-operator (testing only; small dimensions)."""
    import itertools
    d = len(op_cores)
    dims_o = [c.shape[1] for c in op_cores]
    dims_i = [c.shape[2] for c in op_cores]
    n_out = int(np.prod(dims_o))
    n_in = int(np.prod(dims_i))
    mat = np.zeros((n_out, n_in), dtype=complex)
    for i_out in itertools.product(*(range(n) for n in dims_o)):
        for i_in in itertools.product(*(range(n) for n in dims_i)):
            vec = op_cores[0][0, i_out[0], i_in[0], :]
            for k in range(1, d):
                vec = np.tensordot(op_cores[k][:, i_out[k], i_in[k], :],
                                   vec, axes=([0], [0]))
            mat[np.ravel_multi_index(i_out, dims_o),
                np.ravel_multi_index(i_in, dims_i)] = vec[0]
    return mat
