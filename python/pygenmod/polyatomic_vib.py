"""Polyatomic N-mode vibrations: product-DVR variational solve and coupled
normal-mode analytics (multi-atomic expansion step 3, book Sec 3.4).

Hamiltonian in mass-weighted normal coordinates (unit mass in the DVR):

    H = sum_k [ T_k (+) I ] + V(Q),
    V(Q) = sum_k 1/2 w_k^2 Q_k^2 + sum_{k<l} lam_kl Q_k Q_l,

with Colbert-Miller Sinc-DVR kinetic matrices (reusing pygenmod.dvr) and a
potential/bilinear-coupling term diagonal on the product grid. The two-mode
bilinear coupling admits an analytic normal-mode solution (eigenvalues of
the dynamical matrix), serving as a machine-precision acceptance reference
for the variational product-DVR solve.
"""
import numpy as np

from pygenmod.dvr import dvr_sinc_init

POLY_MAX_MODES = 3


def coupled2_normal_freqs(w1: float, w2: float, lam: float) -> tuple[float, float]:
    """Analytic normal-mode frequencies of two bilinearly coupled oscillators.

    The dynamical matrix is K = [[w1^2, lam], [lam, w2^2]]; returns
    (Om_lo, Om_hi) with Om^2 = (w1^2 + w2^2)/2 -/+ sqrt((w1^2-w2^2)^2 + 4 lam^2)/2.
    """
    s = 0.5 * (w1 * w1 + w2 * w2)
    disc = 0.5 * np.sqrt(max((w1 * w1 - w2 * w2) ** 2 + 4.0 * lam * lam, 0.0))
    return (float(np.sqrt(max(s - disc, 1e-12))), float(np.sqrt(s + disc)))


def poly_vib_product_levels(omegas, lam, x_half, n_pts: int) -> np.ndarray:
    """Variational levels of an N-mode (N = 1..3) harmonic + bilinear-coupling
    vibrational Hamiltonian on a product Sinc-DVR grid.

    Parameters
    ----------
    omegas : (N,) harmonic frequencies of the uncoupled modes (a.u.)
    lam    : (N, N) bilinear couplings lam_kl (upper triangle used; a.u.)
    x_half : (N,) half box width per mode (a.u.; ~5-8 harmonic lengths)
    n_pts  : sinc-DVR points per mode (n_pts**N <= 3400 for dense diag)

    Returns
    -------
    e_levels : (n_pts**N,) ascending variational levels (a.u.)
    """
    omegas = np.asarray(omegas, dtype=float)
    lam = np.asarray(lam, dtype=float)
    x_half = np.asarray(x_half, dtype=float)
    n_modes = len(omegas)
    if not 1 <= n_modes <= POLY_MAX_MODES:
        raise ValueError("n_modes must be 1..3")
    if n_pts ** n_modes > 3400:
        raise ValueError("product basis too large for dense diagonalization")
    if np.any(omegas <= 0.0) or np.any(x_half <= 0.0):
        raise ValueError("omegas and x_half must be positive")

    dvrs = [dvr_sinc_init(-x_half[k], x_half[k], n_pts, 1.0)
            for k in range(n_modes)]

    # Kinetic energy = Kronecker SUM over modes (not a product!):
    #   H_kin = sum_k kron(I_slow, T_k, I_fast),
    # with mode 0 on the fastest index, matching the potential mapping below
    # (x_vals[k] uses t[k] since np.ndindex iterates the last axis fastest).
    n_tot = n_pts ** n_modes
    h_mat = np.zeros((n_tot, n_tot))
    # mode 0 = fastest (last ndindex axis), mode N-1 = slowest (first kron factor)
    for k in range(n_modes):
        n_slow = n_pts ** (n_modes - 1 - k)
        n_fast = n_pts ** k
        term = dvrs[k].t_mat
        if n_fast > 1:
            term = np.kron(term, np.eye(n_fast))
        if n_slow > 1:
            term = np.kron(np.eye(n_slow), term)
        h_mat += term

    # Diagonal potential + bilinear couplings on the product grid.
    # Axis mapping: kron(T1, kron(T2, ...)) places mode 1 on the FASTEST
    # flat index, while np.ndindex yields the LAST tuple element fastest —
    # so mode k's grid index is t[n_modes - 1 - k].
    for t in np.ndindex(*([n_pts] * n_modes)):
        pot = 0.0
        x_vals = [dvrs[k].x[t[n_modes - 1 - k]] for k in range(n_modes)]
        for k in range(n_modes):
            pot += 0.5 * omegas[k] ** 2 * x_vals[k] ** 2
        for k in range(n_modes):
            for l in range(k + 1, n_modes):
                pot += lam[k, l] * x_vals[k] * x_vals[l]
        flat = np.ravel_multi_index(t, (n_pts,) * n_modes)
        h_mat[flat, flat] += pot

    return np.linalg.eigvalsh(h_mat)


def morse_bound_energy(v: float, d_e: float, a_morse: float) -> float:
    """Exact Morse bound level measured from the well bottom (mass-weighted).

    G(v) = w0 (v + 1/2) - w0^2 (v + 1/2)^2 / (4 D_e),  w0 = a sqrt(2 D_e);
    the quadratic coefficient w0^2/(4 D_e) is the Dunham constant w_e x_e.
    Fortran mirror: `morse_bound_energy` (src/mod_polyatomic_vib.f90).
    """
    w0 = a_morse * np.sqrt(2.0 * d_e)
    vh = v + 0.5
    return w0 * vh - w0 * w0 * vh * vh / (4.0 * d_e)


def morse_bound_count(d_e: float, a_morse: float) -> int:
    """Number of Morse bound levels: floor(sqrt(2 D_e)/a - 1/2) + 1."""
    return int(np.floor(np.sqrt(2.0 * d_e) / a_morse - 0.5)) + 1


def poly_vib_morse_levels(d_e, a_morse, lam, x_half, n_pts: int) -> np.ndarray:
    """Variational levels of an N-mode (N = 1..3) Morse-anharmonic vibrational
    Hamiltonian on a product Sinc-DVR grid (multi-atomic step 4).

    V(Q) = sum_k D_k (1 - exp(-a_k Q_k))^2 + sum_{k<l} lam_kl Q_k Q_l.

    Acceptance: separable case (lam = 0) reproduces sum_k G_k(v_k) exactly
    (morse_bound_energy); a single mode shows the anharmonic red shift
    E_1 - E_0 < w0 and the bound-level count morse_bound_count.

    Warning: a bilinear coupling lam_kl is incompatible with the saturated
    Morse walls (V is unbounded below along Q_k = -Q_l), so lam != 0 is only
    meaningful for localized states near the well bottom on small boxes;
    near-dissociation states require lam = 0 (or a full quartic force field,
    book roadmap 12.2). Morse tails decay at kappa = sqrt(2 (D - E)) which
    saturates near dissociation, so box-truncation errors of the top levels
    converge with box width, not with grid density.
    """
    d_e = np.asarray(d_e, dtype=float)
    a_morse = np.asarray(a_morse, dtype=float)
    lam = np.asarray(lam, dtype=float)
    x_half = np.asarray(x_half, dtype=float)
    n_modes = len(d_e)
    if not 1 <= n_modes <= POLY_MAX_MODES:
        raise ValueError("n_modes must be 1..3")
    if n_pts ** n_modes > 3400:
        raise ValueError("product basis too large for dense diagonalization")
    if np.any(d_e <= 0.0) or np.any(a_morse <= 0.0) or np.any(x_half <= 0.0):
        raise ValueError("d_e, a_morse and x_half must be positive")

    dvrs = [dvr_sinc_init(-x_half[k], x_half[k], n_pts, 1.0)
            for k in range(n_modes)]
    n_tot = n_pts ** n_modes
    h_mat = np.zeros((n_tot, n_tot))
    for k in range(n_modes):
        n_slow = n_pts ** k
        n_fast = n_pts ** (n_modes - 1 - k)
        term = dvrs[k].t_mat
        if n_fast > 1:
            term = np.kron(term, np.eye(n_fast))
        if n_slow > 1:
            term = np.kron(np.eye(n_slow), term)
        h_mat += term

    for t in np.ndindex(*([n_pts] * n_modes)):
        x_vals = [dvrs[k].x[t[k]] for k in range(n_modes)]
        pot = 0.0
        for k in range(n_modes):
            pot += d_e[k] * (1.0 - np.exp(-a_morse[k] * x_vals[k])) ** 2
        for k in range(n_modes):
            for l in range(k + 1, n_modes):
                pot += lam[k, l] * x_vals[k] * x_vals[l]
        flat = np.ravel_multi_index(t, (n_pts,) * n_modes)
        h_mat[flat, flat] += pot

    return np.linalg.eigvalsh(h_mat)
