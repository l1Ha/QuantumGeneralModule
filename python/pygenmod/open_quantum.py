"""Open-quantum-system spectral methods for the thermalized two-level
Lindblad equation (multi-atomic step 9 / open-system spectral methods).

The vectorized Liouvillian supermatrix, the steady state (augmented solve
cross-checked against the closed form), the damping basis (analytic
eigenmodes verified against the supermatrix action), and the spectral
reconstruction rho(t) = rho_ss + sum_k c_k exp(lambda_k t) R_k.

Fortran mirror: src/mod_open_quantum.f90 routines
`lindblad_supermatrix_2level`, `lindblad_steady_state_2level`,
`lindblad_eigenmodes_2level`, `lindblad_spectral_reconstruct`.
Conventions: H = delta |e><e|; collapse rates gamma_down = gamma (n_th + 1),
gamma_up = gamma n_th; dephasing gamma_phi in the sigma_z convention
(coherence decay rate = Gamma/2 + 2 gamma_phi with Gamma = gamma_down +
gamma_up).
"""

import numpy as np

__all__ = [
    "lindblad_supermatrix_2level",
    "lindblad_steady_state_2level",
    "lindblad_eigenmodes_2level",
    "lindblad_spectral_reconstruct",
]


def lindblad_supermatrix_2level(gamma_down: float, gamma_up: float,
                                gamma_phi: float, delta: float) -> np.ndarray:
    """4x4 complex vectorized Liouvillian in the (ee, eg, ge, gg) ordering."""
    s = 0.5 * (gamma_down + gamma_up) + 2.0 * gamma_phi
    l_super = np.zeros((4, 4), dtype=complex)
    l_super[0, 0] = -gamma_down
    l_super[0, 3] = gamma_up
    l_super[1, 1] = -s - 1j * delta
    l_super[2, 2] = -s + 1j * delta
    l_super[3, 0] = gamma_down
    l_super[3, 3] = -gamma_up
    return l_super


def lindblad_steady_state_2level(gamma_down: float, gamma_up: float) -> np.ndarray:
    """Thermal two-level steady state, closed form rho_ee = gamma_up/Gamma.

    The augmented-solve route (Fortran) reproduces this exactly; here the
    closed form is returned and the eigen-relation L vec(rho_ss) = 0 is
    verified by the caller through `lindblad_eigenmodes_2level`.
    """
    gam_t = gamma_down + gamma_up
    rho_ss = np.zeros((2, 2), dtype=complex)
    rho_ss[0, 0] = gamma_up / gam_t
    rho_ss[1, 1] = gamma_down / gam_t
    return rho_ss


def lindblad_eigenmodes_2level(gamma_down: float, gamma_up: float,
                               gamma_phi: float, delta: float,
                               rho_ss: np.ndarray):
    """Damping basis: analytic eigenpairs verified against the supermatrix.

    Returns (evals, modes, stat) with stat = 0 iff
    || L vec(R_k) - lambda_k vec(R_k) || < 1e-12 for all four modes
    (row-major vectorization (ee, eg, ge, gg), matching the supermatrix).
    """
    l_super = lindblad_supermatrix_2level(gamma_down, gamma_up, gamma_phi, delta)
    s = 0.5 * (gamma_down + gamma_up) + 2.0 * gamma_phi
    evals = np.array([0.0, -(gamma_down + gamma_up),
                      -s - 1j * delta, -s + 1j * delta], dtype=complex)
    modes = np.zeros((2, 2, 4), dtype=complex)
    modes[:, :, 0] = rho_ss
    modes[0, 0, 1] = 1.0
    modes[1, 1, 1] = -1.0
    modes[0, 1, 2] = 1.0          # R_eg = |e><g|
    modes[1, 0, 3] = 1.0          # R_ge = |g><e|
    stat = 0
    for k in range(4):
        vec = np.array([modes[0, 0, k], modes[0, 1, k],
                        modes[1, 0, k], modes[1, 1, k]])
        resid = l_super @ vec - evals[k] * vec
        if np.sum(np.abs(resid)) > 1e-12:
            stat = -(k + 1)
    return evals, modes, stat


def lindblad_spectral_reconstruct(rho0: np.ndarray, rho_ss: np.ndarray,
                                  gamma_down: float, gamma_up: float,
                                  gamma_phi: float, delta: float,
                                  t_now: float) -> np.ndarray:
    """rho(t) = rho_ss + c_pop e^{-Gamma t} R_pop + c_eg e^{lambda_eg t} R_eg
    + c_ge e^{lambda_ge t} R_ge with c_pop = rho_ee(0) - rho_ee^ss,
    c_eg = rho_eg(0), c_ge = rho_ge(0) (R_ss has trace 1, R_pop traceless)."""
    s = 0.5 * (gamma_down + gamma_up) + 2.0 * gamma_phi
    c_pop = float(np.real(rho0[0, 0]) - np.real(rho_ss[0, 0]))
    rho_t = rho_ss.copy()
    decay_p = np.exp(-(gamma_down + gamma_up) * t_now)
    rho_t[0, 0] += c_pop * decay_p
    rho_t[1, 1] -= c_pop * decay_p
    rho_t[0, 1] += rho0[0, 1] * np.exp(-(s + 1j * delta) * t_now)
    rho_t[1, 0] += rho0[1, 0] * np.exp(-(s - 1j * delta) * t_now)
    return rho_t
