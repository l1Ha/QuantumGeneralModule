"""
Minimal MCTDH core: coefficient (A) and single-particle-function (SPF)
propagation for SOP Hamiltonians.

Python reference implementation of book section 18.2 and the proposed
`mod_mctdh_core` (book section 18.8):

    Psi = sum_{j_1..j_D} A_{j_1..j_D}(t) prod_kappa phi_{j_kappa}^(kappa)(t)

Equations of motion (Dirac-Frenkel variational principle in the natural
gauge <phi_i^(k)|phi-dot_j^(k)> = 0, cf. Lubich, "From Quantum to Classical
Molecular Dynamics"; book section 18.2 notes the gauge freedom):

    i hbar A-dot = H_A A,
    i hbar phi-dot^(k) = (1 - P^(k)) [rho^(k)]^+ (H^(k) phi^(k)),

with H_A acting through the SPF matrices eta^(k)_alpha = phi h_alpha phi^dagger,
the mean-field operator H^(k) = sum_alpha G^(k)_alpha h_alpha^(k), and the
mode-k density matrix rho^(k). rho^(k) can be singular; we use a
pseudo-inverse with an explicit rcond (book section 18.2 regularization note).
Integration is classical RK4; SPFs are re-orthonormalized after each step
(with the QR factors absorbed into A) to control round-off drift.

Fortran counterpart (roadmap): `mod_mctdh_core` with `mctdh_config_t`,
`mctdh_state_t`, `mctdh_propagate_step`, `mctdh_mean_field`.

References:
    U. Manthe, H.-D. Meyer, L. S. Cederbaum, J. Chem. Phys. 97, 3199 (1992).
    M. H. Beck, A. Jackle, G. A. Worth, H.-D. Meyer, Phys. Rep. 324, 1 (2000).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np

__all__ = ["MCTDHConfig", "MCTDHState", "mctdh_mean_field",
           "mctdh_propagate_step", "mctdh_propagate"]


@dataclass
class MCTDHConfig:
    """Mirrors `mctdh_config_t` of book section 18.8.

    rho_floor: relative occupation floor for the regularized inverse of
    rho^(k) (book section 18.2 regularization note). Eigenvalues of rho below
    rho_floor * lambda_max are zeroed in the inverse, which suppresses the
    gauge singularity of the (A, SPF) parametrization for weakly occupied
    SPFs; step sizes must respect dt * ||H|| << 1 in any case.
    """
    dt: float = 1.0              # RK4 step (a.u.)
    n_steps: int = 100
    rho_floor: float = 1.0e-8    # relative occupation floor for rho^(k)^+
    reorthonormalize: bool = True


@dataclass
class MCTDHState:
    """Mirrors `mctdh_state_t` of book section 18.8.

    A    : coefficient tensor (n_spf_1, ..., n_spf_D), complex
    spf  : per mode, array (n_spf_k, n_prim_k); rows are the SPF vectors in the
           primitive (DVR/FGH) basis and are kept orthonormal.
    """
    A: np.ndarray
    spf: List[np.ndarray]
    time: float = 0.0

    def norm(self) -> float:
        return float(np.linalg.norm(self.A))


def _eta_matrix(spf: np.ndarray, h: np.ndarray) -> np.ndarray:
    """eta[i, j] = <phi_i| h |phi_j> for one mode (SPF vectors are rows)."""
    return spf.conj() @ h @ spf.T


def _apply_along_axis(mat: np.ndarray, tensor: np.ndarray, axis: int) -> np.ndarray:
    """Contract mat (n_out, n_in) with tensor along `axis` (output keeps place)."""
    return np.moveaxis(np.tensordot(mat, tensor, axes=([1], [axis])), 0, axis)


def _spf_apply(spf: np.ndarray, fac: np.ndarray) -> np.ndarray:
    """SPF-space matrix of one SOP factor (matrix or diagonal vector)."""
    if fac.ndim == 1:
        return (spf.conj() * fac[None, :]) @ spf.T
    return _eta_matrix(spf, fac)


def _regularized_inverse(rho: np.ndarray, floor: float) -> np.ndarray:
    """rho^+ with a relative occupation floor (book section 18.2).

    rho = U diag(w) U^dagger is inverted as U diag(w_inv) U^dagger with
    w_inv = 1/w for w > floor * w_max and 0 otherwise. Directions with
    occupations below the floor carry no dynamical information (the A tensor
    does not weight them), so their mean-field contribution is dropped
    instead of being amplified by 1/w.
    """
    w, u = np.linalg.eigh(0.5 * (rho + rho.conj().T))
    w_max = float(w.max()) if w.size else 0.0
    w_inv = np.zeros_like(w)
    keep = w > floor * max(w_max, 1.0e-300)
    w_inv[keep] = 1.0 / w[keep]
    return (u * w_inv) @ u.conj().T


def _matricize(a: np.ndarray, k: int) -> np.ndarray:
    """Move axis k to the end and flatten the rest: (R, n_k)."""
    return np.moveaxis(np.asarray(a), k, -1).reshape(-1, np.asarray(a).shape[k])


def _mean_field_operators(state: MCTDHState, sop, etas=None,
                          rho_floor: float = 1.0e-8):
    """Mean-field data for the SPF equations (book section 18.2).

    For mode k and SOP term alpha the mean-field coefficient matrix is the
    single-hole expectation

        G^(k,alpha)[m, l] = sum_r conj(A)_{r m} B^(alpha)_{r l},
        B^(alpha) = A with eta^(lambda)_alpha applied on every mode lambda != k,

    and the mean-field operator acts as
    (H^(k) phi)_m = sum_alpha c_alpha sum_l G^(k,alpha)[m, l] h^(k)_alpha phi_l.

    Returns (G, W) with G[k][alpha] the (n_spf, n_spf) matrices and W[k] the
    (n_spf, n_prim) projected, rho^{+}-weighted derivative columns
    (i hbar phi-dot = (1 - P) rho^+ H^(k) phi).

    ``etas`` may carry the SPF-space matrices eta^(k)_alpha shared with the
    A-equation so they are built only once per time derivative.
    """
    d = len(sop.dims)
    if etas is None:
        etas = [[_spf_apply(state.spf[k], term.factors[k]) for k in range(d)]
                for term in sop.terms]
    rhos = _density_matrices(state)
    G, W = [], []
    for k in range(d):
        a_mat = _matricize(state.A, k)                    # (R, n_k)
        g_list = []
        w = np.zeros((sop.dims[k], state.spf[k].shape[0]), dtype=complex)
        for alpha, term in enumerate(sop.terms):
            c = state.A
            for m in range(d):
                if m != k:
                    c = _apply_along_axis(etas[alpha][m], c, m)
            g_mat = a_mat.conj().T @ _matricize(c, k)     # (n_k, n_k)
            g_list.append(g_mat)
            fac = term.factors[k]
            if fac.ndim == 1:
                hphi = (state.spf[k] * np.asarray(fac)[None, :]).T
            else:
                hphi = np.asarray(fac) @ state.spf[k].T
            w += term.coeff * (hphi @ g_mat.T)
        # (1 - P) projection, then rho^(k)^{+} on the SPF index
        w = w - state.spf[k].T @ (state.spf[k].conj() @ w)
        rho_p = _regularized_inverse(rhos[k], rho_floor)
        w = rho_p @ w.T                                   # (n_spf, n_prim)
        G.append(g_list)
        W.append(w)
    return G, W


def mctdh_mean_field(state: MCTDHState, sop, rho_floor: float = 1.0e-8):
    """Public mean-field hook (book 18.8 `mctdh_mean_field`).

    Returns the per-mode/per-term single-hole matrices G^(k)_alpha and the
    projected derivative columns W^(k) (i hbar phi-dot = (1-P) rho^+ H phi).
    """
    return _mean_field_operators(state, sop, rho_floor)


def _density_matrices(state: MCTDHState) -> List[np.ndarray]:
    """rho^(k)[j, l] = sum over all other A indices of conj(A) A."""
    d = state.A.ndim
    letters = "abcdefghijklmnopqrstuv"
    rhos = []
    for k in range(d):
        idx_a = [letters[m] for m in range(d)]
        idx_b = [letters[m] for m in range(d)]
        idx_a[k] = "z"
        idx_b[k] = "y"
        subscripts = "".join(idx_a) + "," + "".join(idx_b) + "->zy"
        rhos.append(np.einsum(subscripts, state.A.conj(), state.A))
    return rhos


def _derivative(state: MCTDHState, sop, rho_floor: float = 1.0e-8):
    """Time derivatives (A-dot, phi-dot list) in the natural gauge.

    i hbar A-dot = (eta-form contraction of H with A),
    i hbar phi-dot^(k) = (1 - P^(k)) [rho^(k)]^+ H^(k) phi^(k).
    """
    d = len(sop.dims)
    etas = []
    for term in sop.terms:
        etas.append([_spf_apply(state.spf[k], term.factors[k]) for k in range(d)])
    # A-dot = -i H_A A  (H_A acts via the eta matrices)
    a_dot = np.zeros_like(state.A)
    for term, eta in zip(sop.terms, etas):
        c = state.A
        for k in range(d):
            c = _apply_along_axis(eta[k], c, k)
        a_dot += term.coeff * c
    a_dot = -1j * a_dot

    G, W = _mean_field_operators(state, sop, etas=etas, rho_floor=rho_floor)
    spf_dots = [-1j * W[k] for k in range(d)]      # W rows = SPF index
    return a_dot, spf_dots


def _reorthonormalize(state: MCTDHState) -> None:
    """QR the SPFs and absorb the R factors into A (gauge-conserving)."""
    d = state.A.ndim
    for k in range(d):
        q, rr = np.linalg.qr(state.spf[k].T, mode="reduced")   # q: (n_prim, n_spf)
        state.spf[k] = q.T
        # old spf rows = R^T maps new rows: phi_old_j = sum_i R[j, i] phi_new_i
        state.A = _apply_along_axis(rr, state.A, k)


def mctdh_propagate_step(state: MCTDHState, sop, cfg: MCTDHConfig) -> MCTDHState:
    """One classical RK4 step of the coupled A/SPF equations (book 18.8)."""
    h = cfg.dt

    def add(s, da, dph, scale):
        return MCTDHState(A=s.A + scale * da,
                          spf=[s.spf[k] + scale * dph[k] for k in range(len(s.spf))],
                          time=s.time)

    fl = cfg.rho_floor
    k1_a, k1_p = _derivative(state, sop, fl)
    k2_a, k2_p = _derivative(add(state, k1_a, k1_p, 0.5 * h), sop, fl)
    k3_a, k3_p = _derivative(add(state, k2_a, k2_p, 0.5 * h), sop, fl)
    k4_a, k4_p = _derivative(add(state, k3_a, k3_p, h), sop, fl)
    state.A = state.A + (h / 6.0) * (k1_a + 2.0 * k2_a + 2.0 * k3_a + k4_a)
    for k in range(len(state.spf)):
        state.spf[k] = state.spf[k] + (h / 6.0) * (
            k1_p[k] + 2.0 * k2_p[k] + 2.0 * k3_p[k] + k4_p[k])
    state.time += h
    if cfg.reorthonormalize:
        _reorthonormalize(state)
    return state


def mctdh_propagate(state: MCTDHState, sop, cfg: MCTDHConfig) -> MCTDHState:
    """Propagate for cfg.n_steps RK4 steps (book 18.8)."""
    for _ in range(cfg.n_steps):
        mctdh_propagate_step(state, sop, cfg)
    return state
