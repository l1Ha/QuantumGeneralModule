"""
Dense linear algebra reference implementations.

`diag_hermitian_matrix` mirrors `mod_linear_algebra::diag_hermitian_matrix`
(src/mod_linear_algebra.f90): the complex Hermitian eigenproblem is reduced to
a REAL symmetric one of twice the size via  A = B + iC  (B symmetric,
C antisymmetric)  ->  M = [[B, -C], [C, B]],  whose spectrum contains every
eigenvalue of A exactly twice; eigenvectors are reconstructed as w = u + iv
from the real eigenvector [u; v] of M. For the Fortran side this reuses the
verified Householder-QL path with zero new numerical machinery; here we call
``numpy.linalg.eigh`` on the doubled real matrix so the reduction itself is
what gets validated against the direct complex eigensolver in the tests.
"""
import numpy as np


def diag_hermitian_matrix(a):
    """Eigen-decompose a complex Hermitian matrix by doubling reduction.

    Parameters
    ----------
    a : (n, n) ndarray (complex)
        Hermitian matrix (only the Hermitian part is used; the real part is
        symmetrized and the imaginary part antisymmetrized).

    Returns
    -------
    d : (n,) ndarray, real, ascending eigenvalues
    z : (n, n) ndarray, complex eigenvectors as columns, unit norm,
        phase fixed so that the largest-magnitude component is real positive.
    """
    a = np.asarray(a, dtype=complex)
    n = a.shape[0]
    if a.shape != (n, n):
        raise ValueError("square matrix required")

    b = 0.5 * (a.real + a.real.T)          # symmetric real part
    c = 0.5 * (a.imag - a.imag.T)          # antisymmetric imaginary part

    m = np.zeros((2 * n, 2 * n))
    m[:n, :n] = b
    m[n:, n:] = b
    m[:n, n:] = -c
    m[n:, :n] = c                          # M = [[B, -C], [C, B]] is symmetric
    # (the (2,1) block is +C because (-C)^T = +C for antisymmetric C)

    d2, z2 = np.linalg.eigh(m)

    # Robust selection of n independent eigenvectors: ANY eigenvector [u; v]
    # of the doubled matrix M maps to a valid eigenvector w = u + i v of A
    # (the eigen-equations hold for every vector of the eigenspace), so:
    #   1) build candidates w_j from ALL 2n columns of z2;
    #   2) group by degenerate eigenvalue (group size g = 2m, with m the
    #      multiplicity in A);
    #   3) modified Gram-Schmidt within each group in order, discarding
    #      collapsed candidates and keeping the first m independent ones.
    # This does not rely on LAPACK's ordering of degenerate eigenvectors
    # (a naive "every other" selection failed for 4-fold degeneracies where
    # the basis ordering yields parallel vectors).
    cand = z2[:n, :] + 1j * z2[n:, :]      # (n, 2n)
    d = np.empty(n)
    z = np.empty((n, n), dtype=complex)
    n_keep = 0
    j = 0
    while j < 2 * n:
        tol_deg = 1.0e-10 * max(1.0, abs(d2[j]))
        kend = j
        while kend + 1 < 2 * n and abs(d2[kend + 1] - d2[j]) <= tol_deg:
            kend += 1
        m_keep = (kend - j + 1) // 2
        n_surv = 0
        for jj in range(j, kend + 1):
            if n_surv >= m_keep:
                break
            v = cand[:, jj].copy()
            for kk in range(n_keep):
                v = v - np.vdot(z[:, kk], v) * z[:, kk]
            nrm = np.linalg.norm(v)
            if nrm > 1.0e-8:
                z[:, n_keep] = v / nrm
                d[n_keep] = d2[j]
                n_keep += 1
                n_surv += 1
        j = kend + 1

    # unit norm + phase convention (largest-magnitude component real positive)
    for k in range(n):
        nrm = np.linalg.norm(z[:, k])
        if nrm > 0.0:
            z[:, k] /= nrm
        imax = int(np.argmax(np.abs(z[:, k])))
        amax = abs(z[imax, k])
        if amax > 0.0:
            ph = z[imax, k] / amax
            z[:, k] *= np.conj(ph)
    return d, z


def peierls_ring_spectrum(n_sites, flux, t_hop=1.0):
    """Nearest-neighbour ring with uniform Peierls phase (closed form helper).

    Returns the exact eigenvalues ``-2 t cos(2 pi (k + flux) / N)`` (ascending)
    of the N-site ring Hamiltonian  H_{i,i+1} = -t exp(i phi)  with
    ``phi = 2 pi flux / N``  (flux in units of the flux quantum).
    """
    lam = -2.0 * t_hop * np.cos(2.0 * np.pi * (np.arange(n_sites) + flux) / n_sites)
    return np.sort(lam)


def peierls_ring_hamiltonian(n_sites, flux, t_hop=1.0):
    """Dense N-site Peierls ring Hamiltonian (Hermitian, complex)."""
    phi = 2.0 * np.pi * flux / n_sites
    h = np.zeros((n_sites, n_sites), dtype=complex)
    hop = -t_hop * np.exp(1j * phi)
    for i in range(n_sites):
        j = (i + 1) % n_sites
        h[i, j] = hop
        h[j, i] = np.conj(hop)
    return h


# ---------------------------------------------------------------------------
# Cholesky / QR / SVD reference wrappers.
#
# Unlike the from-scratch Fortran implementations
# (src/mod_linear_algebra.f90: cholesky_decompose, qr_decompose,
# least_squares_qr, svd_decompose, svd_pseudoinverse), the Python side calls
# numpy/LAPACK: it serves as the INDEPENDENT third-party reference for
# cross-validation (the Fortran tests assert closed forms and the classic
# Hilbert conditioning story; the Python tests assert agreement with LAPACK).
# ---------------------------------------------------------------------------


def cholesky_decompose(a):
    """Lower Cholesky factor L with A = L L^T (numpy/LAPACK reference)."""
    a = np.asarray(a, dtype=float)
    return np.linalg.cholesky(a)


def cholesky_solve(a, b):
    """Solve A x = b for SPD A via the Cholesky factor."""
    l = np.linalg.cholesky(np.asarray(a, dtype=float))
    y = np.linalg.solve(l, np.asarray(b, dtype=float))
    return np.linalg.solve(l.T, y)


def qr_decompose(a):
    """Thin QR (reduced mode): returns Q (m, n), R (n, n) with A = Q R."""
    a = np.asarray(a, dtype=float)
    q, r = np.linalg.qr(a, mode="reduced")
    # sign convention: make the diagonal of R non-negative (unique QR)
    d = np.sign(np.diag(r))
    d[d == 0] = 1.0
    q = q * d
    r = r * d[:, None]
    return q, r


def least_squares_qr(a, b):
    """Least-squares solution via QR (returns x and the residual norm)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    q, r = qr_decompose(a)
    x = np.linalg.solve(r, q.T @ b)
    return x, float(np.linalg.norm(a @ x - b))


def svd_decompose(a):
    """SVD A = U diag(s) V^T with descending non-negative singular values."""
    a = np.asarray(a, dtype=float)
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    return u, s, vt


def svd_pseudoinverse(a, b, rtol=1e-12):
    """Minimum-norm least-squares solution via the SVD pseudoinverse."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    u, s, vt = svd_decompose(a)
    keep = s > rtol * max(1.0, float(s[0]) if s.size else 1.0)
    x = vt[keep].T @ ((u[:, keep].T @ b) / s[keep]) if np.any(keep) else np.zeros(a.shape[1])
    return x, float(np.linalg.norm(a @ x - b))


def hilbert_matrix(m, n):
    """Hilbert matrix H_ij = 1/(i+j-1) (1-based indexing convention)."""
    i, j = np.meshgrid(np.arange(1, m + 1), np.arange(1, n + 1), indexing="ij")
    return 1.0 / (i + j - 1.0)
