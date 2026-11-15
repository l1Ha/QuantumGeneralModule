"""
Krylov (Lanczos) Short-Iterative Propagator — Python reference implementation.

Mirrors `mod_krylov_propagator` (src/mod_krylov_propagator.f90):
Lanczos tridiagonalization with two-pass full reorthogonalization, exact
matrix exponential of the tridiagonal T within the Krylov subspace, and the
a-posteriori error estimate beta_{m+1} |[exp(-i T dt)]_{m,1}|.

Requires a real-symmetric Hamiltonian (no magnetic field, no complex
absorbing potential); for complex potentials use the Chebyshev /
split-operator paths.
"""
import numpy as np


def lanczos_tridiag(r_mult, v_start, m_max, tol_break=1.0e-14):
    """m-step Lanczos tridiagonalization of a real-symmetric operator.

    Parameters
    ----------
    r_mult : callable
        Real matvec: ``r_mult(v) -> H @ v`` (real symmetric H).
    v_start : ndarray
        Real starting vector (nonzero).
    m_max : int
        Maximum Krylov dimension.
    tol_break : float
        Breakdown tolerance (relative to ||v_start||); beta below this
        signals an invariant subspace and terminates early.

    Returns
    -------
    alpha : ndarray (m_eff,)
    beta : ndarray (m_eff + 1,)
        beta[0] = 0; beta[m_eff] is the final residual norm (0 on exact
        invariant-subspace termination).
    v_basis : ndarray (n, m_eff)
        Orthonormal Krylov basis.
    m_eff : int
    """
    v_start = np.asarray(v_start, dtype=float)
    n = v_start.size
    nrm = np.linalg.norm(v_start)
    if nrm <= np.finfo(float).tiny:
        raise ValueError("zero starting vector")

    alpha = np.zeros(m_max)
    beta = np.zeros(m_max + 1)
    v_basis = np.zeros((n, m_max))
    v_basis[:, 0] = v_start / nrm
    v_prev = np.zeros(n)
    tol = tol_break * nrm

    m_eff = m_max
    for j in range(m_max):
        w = np.asarray(r_mult(v_basis[:, j]), dtype=float)
        alpha[j] = float(v_basis[:, j] @ w)
        w = w - alpha[j] * v_basis[:, j]
        if j > 0:
            w = w - beta[j] * v_prev
        # two-pass classical Gram-Schmidt (V^T V = I to machine precision)
        coeff = v_basis[:, :j + 1].T @ w
        w = w - v_basis[:, :j + 1] @ coeff
        coeff = v_basis[:, :j + 1].T @ w
        w = w - v_basis[:, :j + 1] @ coeff
        nrm = float(np.linalg.norm(w))
        beta[j + 1] = nrm
        if nrm <= tol:
            m_eff = j + 1
            break
        if j < m_max - 1:
            v_basis[:, j + 1] = w / nrm
            v_prev = v_basis[:, j]

    return alpha[:m_eff], beta[:m_eff + 1], v_basis[:, :m_eff].copy(), m_eff


def _lanczos_apply_expm(r_mult, v, nrm, dt, m_max):
    """out = nrm * V_m exp(-i T_m dt) e_1 for a real vector v (complex out)."""
    alpha, beta, v_basis, m_eff = lanczos_tridiag(r_mult, v, m_max)
    tri = np.diag(alpha)
    if m_eff > 1:
        off = beta[1:m_eff]
        tri += np.diag(off, 1) + np.diag(off, -1)
    lam, umat = np.linalg.eigh(tri)
    coeff = np.exp(-1j * lam * dt) * umat[0, :]
    # a-posteriori error estimate beta_{m+1} |[exp(-i T dt)]_{m,1}|
    err_est = float(beta[m_eff] * abs(np.dot(umat[m_eff - 1, :], coeff)))
    if beta[m_eff] <= 0.0:
        err_est = 0.0
    c = umat @ coeff
    out = nrm * (v_basis @ c)
    return out, err_est, m_eff


def krylov_propagate_step(psi, dt, m_krylov, r_mult):
    """One Lanczos-Krylov step: psi <- exp(-i H dt) psi.

    Complex psi is split into real and imaginary parts (H real symmetric).

    Returns
    -------
    err_est : float
        A-posteriori error estimate.
    m_used : int
        Krylov dimension actually used.
    """
    psi = np.asarray(psi, dtype=complex)
    n = psi.size
    mm = max(2, int(m_krylov))
    v_re = psi.real.copy()
    v_im = psi.imag.copy()
    nrm_re = np.linalg.norm(v_re)
    nrm_im = np.linalg.norm(v_im)
    err_re = err_im = 0.0
    m_re = m_im = 0

    out = np.zeros(n, dtype=complex)
    if nrm_re > np.finfo(float).tiny:
        part, err_re, m_re = _lanczos_apply_expm(r_mult, v_re, nrm_re, dt, mm)
        out += part
    if nrm_im > np.finfo(float).tiny:
        part, err_im, m_im = _lanczos_apply_expm(r_mult, v_im, nrm_im, dt, mm)
        out += 1j * part

    psi[...] = out
    return max(err_re, err_im), max(m_re, m_im)
