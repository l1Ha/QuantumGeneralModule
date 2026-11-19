"""
Non-symmetric Arnoldi — Python mirror of `mod_arnoldi`
(src/mod_arnoldi.f90).

Restarted Arnoldi for k target eigenpairs of a general complex matrix
(ARPACK naming LM/SM/LR/SR/LI/SI). Design mirrors the Fortran module:

* m-step Arnoldi factorization (modified Gram-Schmidt with a second
  reorthogonalization pass);
* the small Hessenberg eigenproblem is solved by ``numpy.linalg.eig`` —
  the independent reference for the Fortran shifted-Hessenberg-QR engine
  (Wilkinson shift + exceptional shift with Schur-vector accumulation);
* residual a-posteriori estimate ``r_i = |h_{m+1,m}| |y_i(m)|`` where
  ``y_i`` is the eigenvector of the Hessenberg matrix in the ARNOLDI
  basis (the Schur back-substitution vector has ``y(m) = 0`` for i < m —
  using it directly is a classic mistake; the Fortran side transforms
  through the Schur vectors Z first);
* locking + restart orthogonalization: converged pairs enter a locked
  table (cosine-deduplicated at |<x_i, x_kept>| > 0.999 — scale-free
  detection of the asymptotically duplicated Ritz pairs a restart Krylov
  space produces for the same eigenvalue); the restart vector is the
  unconverged best Ritz vector plus geometrically decaying noise,
  orthogonalized against the locked directions (implicit deflation).
"""
import numpy as np

__all__ = ["arnoldi_factorize", "arnoldi_general_eigenpairs"]


def arnoldi_factorize(c_mult, n, m_max, v_start):
    """m-step Arnoldi factorization A V = V H + h_sub v_{m+1} e_m^T.

    Returns (v_basis (n, m_eff), h_hess (m_eff, m_eff), h_sub, n_matvec).
    """
    v_start = np.asarray(v_start, dtype=complex)
    nrm0 = np.linalg.norm(v_start)
    if nrm0 <= 1e-300:
        raise ValueError("zero starting vector")
    tol_break = 1.0e-13 * np.sqrt(max(n, 1))

    v_basis = np.zeros((n, m_max), dtype=complex)
    h_hess = np.zeros((m_max, m_max), dtype=complex)
    v_basis[:, 0] = v_start / nrm0
    n_mv = 0
    h_sub = 0.0
    for j in range(m_max):
        w = np.asarray(c_mult(v_basis[:, j]), dtype=complex)
        n_mv += 1
        # MGS + second pass when the first-pass projection is large
        for pass2 in range(2):
            reorth_max = 0.0
            for i in range(j + 1):
                proj = np.vdot(v_basis[:, i], w)
                reorth_max = max(reorth_max, abs(proj))
                h_hess[i, j] += proj
                w = w - proj * v_basis[:, i]
            if pass2 == 0 and reorth_max > 1.0e-9:
                continue
            break
        h_sub = np.linalg.norm(w)
        if h_sub <= tol_break:
            # invariant subspace: legal early exit
            return v_basis[:, :j + 1], h_hess[:j + 1, :j + 1], 0.0, n_mv
        if j == m_max - 1:
            return v_basis, h_hess, h_sub, n_mv
        v_basis[:, j + 1] = w / h_sub
        h_hess[j + 1, j] = h_sub
    raise RuntimeError("unreachable")


def _target_key(lam, which):
    if which == 'LM':
        return -abs(lam)
    if which == 'SM':
        return abs(lam)
    if which == 'SR':
        return lam.real
    if which == 'LI':
        return -lam.imag
    if which == 'SI':
        return lam.imag
    return -lam.real            # 'LR'


def arnoldi_general_eigenpairs(c_mult, n, k, m_max, n_restart, which='LR',
                               tol=1.0e-10, seed=42):
    """Restarted Arnoldi for k target eigenpairs (ARPACK naming).

    Returns (evals (k,), rvecs (n, k), res_max, n_matvec, stat) with
    stat 0 = all k pairs converged to residual < tol.
    """
    rng = np.random.default_rng(seed)
    grid = np.arange(1, n + 1)
    v_start = rng.random(n) - 0.5 + 1j * (rng.random(n) - 0.5) \
        + 0.13 * np.sin(grid) + 0.11j * np.cos(grid)
    v_start /= np.linalg.norm(v_start)

    lock_lam, lock_vec, lock_res = [], [], []
    open_lam, open_vec, open_res = [], [], []
    n_matvec = 0
    stat = 1

    for pass_i in range(max(n_restart, 1)):
        v_basis, h_hess, h_sub, n_mv = arnoldi_factorize(c_mult, n, m_max,
                                                         v_start)
        n_matvec += n_mv
        m_eff = h_hess.shape[0]
        lam_h, y_h = np.linalg.eig(h_hess)
        resid = abs(h_sub) * np.abs(y_h[m_eff - 1, :]) if h_sub > 0 \
            else np.zeros(m_eff)
        ritz = v_basis @ y_h
        order = np.argsort([_target_key(l, which) for l in lam_h])

        for q in range(min(k, m_eff)):
            i = order[q]
            cv = ritz[:, i] / np.linalg.norm(ritz[:, i])
            # parallel to a locked pair? upgrade it if better
            dup = False
            for p in range(len(lock_lam)):
                if abs(np.vdot(lock_vec[p], cv)) > 0.999:
                    if resid[i] < lock_res[p]:
                        lock_lam[p] = lam_h[i]
                        lock_vec[p] = cv
                        lock_res[p] = resid[i]
                    dup = True
                    break
            if dup:
                continue
            if resid[i] < tol and len(lock_lam) < k:
                lock_lam.append(lam_h[i])
                lock_vec.append(cv)
                lock_res.append(resid[i])
                continue
            # unconverged: merge into the open set (cosine-upgrade, key order)
            merged = False
            for p in range(len(open_lam)):
                if abs(np.vdot(open_vec[p], cv)) > 0.999:
                    if resid[i] < open_res[p]:
                        open_lam[p] = lam_h[i]
                        open_vec[p] = cv
                        open_res[p] = resid[i]
                    merged = True
                    break
            if not merged:
                if len(open_lam) < k:
                    open_lam.append(lam_h[i])
                    open_vec.append(cv)
                    open_res.append(resid[i])
                else:
                    worst = max(range(len(open_lam)),
                                key=lambda t: _target_key(open_lam[t], which))
                    if _target_key(lam_h[i], which) < \
                            _target_key(open_lam[worst], which):
                        open_lam[worst] = lam_h[i]
                        open_vec[worst] = cv
                        open_res[worst] = resid[i]
            idx = sorted(range(len(open_lam)),
                         key=lambda t: _target_key(open_lam[t], which))
            open_lam = [open_lam[t] for t in idx]
            open_vec = [open_vec[t] for t in idx]
            open_res = [open_res[t] for t in idx]

        if len(lock_lam) >= k:
            stat = 0
            break

        # purge open entries parallel to locked ones (stale residuals)
        keep = [t for t in range(len(open_lam))
                if all(abs(np.vdot(lock_vec[p], open_vec[t])) <= 0.999
                       for p in range(len(lock_lam)))]
        open_lam = [open_lam[t] for t in keep]
        open_vec = [open_vec[t] for t in keep]
        open_res = [open_res[t] for t in keep]

        # restart vector: unconverged best + decaying noise, orthogonalized
        # against locked directions (implicit deflation)
        v_start = np.zeros(n, dtype=complex)
        for p in range(len(open_lam)):
            v_start += open_vec[p]
        noise_amp = 1.0e-2 * (0.1 ** pass_i)
        v_start = v_start + noise_amp * (rng.random(n) - 0.5
                                         + 1j * (rng.random(n) - 0.5))
        for p in range(len(lock_lam)):
            v_start = v_start - np.vdot(lock_vec[p], v_start) * lock_vec[p]
        nrm = np.linalg.norm(v_start)
        if nrm <= 1e-300:
            v_start = rng.random(n) - 0.5 + 1j * (rng.random(n) - 0.5)
            for p in range(len(lock_lam)):
                v_start = v_start - np.vdot(lock_vec[p], v_start) * lock_vec[p]
            nrm = np.linalg.norm(v_start)
        v_start = v_start / nrm

    n_out = min(k, len(lock_lam) + len(open_lam))
    evals = np.array(lock_lam + open_lam[:max(k - len(lock_lam), 0)][:n_out],
                     dtype=complex)[:n_out]
    rvecs = np.array(lock_vec + open_vec[:max(k - len(lock_lam), 0)][:n_out],
                     dtype=complex).T if n_out else \
        np.zeros((n, 0), dtype=complex)
    # res_max 只统计实际返回的对 (锁定 + 填充用的 open 前段)
    n_fill = max(k - len(lock_lam), 0)
    returned_res = lock_res + open_res[:n_fill]
    res_max = max(returned_res) if returned_res else 1e30
    if res_max <= tol and len(lock_lam) >= k:
        stat = 0
    return evals, rvecs, res_max, n_matvec, stat
