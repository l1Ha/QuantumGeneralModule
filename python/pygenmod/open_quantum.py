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




def lindblad_supermatrix_driven(omega: float, delta: float, gamma: float) -> np.ndarray:
    """4x4 complex vectorized Liouvillian of the DRIVEN dissipative two-level
    system (rotating frame H = -delta sigma_z/2 + Omega sigma_x/2, collapse
    Gamma via sigma^-). Row-major vectorization (ee, eg, ge, gg).
    Fortran mirror: `lindblad_supermatrix_driven`.
    """
    l_super = np.zeros((4, 4), dtype=complex)
    h_rot = np.array([[-delta / 2, omega / 2], [omega / 2, delta / 2]],
                     dtype=complex)
    s_m = np.array([[0, 0], [1, 0]], dtype=complex)
    s_dag = s_m.conj().T
    for j in range(2):
        for i in range(2):
            e = np.zeros((2, 2), dtype=complex)
            e[i, j] = 1.0
            lm = -1j * (h_rot @ e - e @ h_rot) + \
                gamma * (s_m @ e @ s_dag - 0.5 * (s_dag @ s_m @ e + e @ s_dag @ s_m))
            for k in range(2):
                for l in range(2):
                    l_super[k * 2 + l, i * 2 + j] = lm[k, l]
    return l_super


def optical_bloch_steady_state_driven(omega: float, delta: float,
                                      gamma: float) -> np.ndarray:
    """Driven-dissipative two-level steady state (rotating frame).

    H = -delta sigma_z/2 + Omega sigma_x/2, collapse Gamma (sigma^-).
    Saturation closed form rho_ee = s/(2(1+s)) with
    s = 2 Omega^2 / (Gamma^2 + 4 delta^2); solved here by the augmented
    linear system (L vec = 0, trace = 1) with the closed-form check.
    Fortran mirror: `optical_bloch_steady_state_driven`.
    """
    l_super = lindblad_supermatrix_driven(omega, delta, gamma)
    aug = np.zeros((5, 5), dtype=complex)
    aug[:4, :4] = l_super
    aug[4, :4] = 1.0
    aug[:4, 4] = 1.0
    x = np.linalg.solve(aug, np.array([0, 0, 0, 0, 1], dtype=complex))
    s_sat = 2.0 * omega**2 / (gamma**2 + 4.0 * delta**2)
    rho_ss = np.array([[x[0], x[1]], [x[2], x[3]]])
    if abs(rho_ss[0, 0].real - s_sat / (2.0 * (1.0 + s_sat))) > 1e-10:
        raise ValueError("driven steady state inconsistent with closed form")
    return rho_ss


def first_order_correlation_2level(omega: float, delta: float, gamma: float,
                                   rho_ss: np.ndarray, tau_max: float,
                                   n_tau: int):
    """Quantum regression theorem: g1(tau) = <sigma+(tau) sigma-(0)>.

    Propagates M(0) = sigma^- rho_ss under the driven Liouvillian by
    vectorized RK4 (no Hermitization/trace normalization — M(tau) is not
    a density matrix). g1(tau) = tr[sigma+ M(tau)] = M[1, 0] (0-based).
    Returns (tau_grid, g1_grid).
    Fortran mirror: `first_order_correlation_2level`.
    """
    l_super = lindblad_supermatrix_driven(omega, delta, gamma)
    tau_grid = tau_max * np.arange(n_tau) / (n_tau - 1.0)
    m0 = np.array([[0, 0], [1, 0]], dtype=complex) @ rho_ss
    m_vec = np.array([m0[0, 0], m0[0, 1], m0[1, 0], m0[1, 1]])
    dt = tau_max / (n_tau - 1.0) / 40.0
    g1 = [m_vec[2]]
    for _ in range(n_tau - 1):
        for _ in range(40):
            k1 = l_super @ m_vec
            k2 = l_super @ (m_vec + 0.5 * dt * k1)
            k3 = l_super @ (m_vec + 0.5 * dt * k2)
            k4 = l_super @ (m_vec + dt * k3)
            m_vec = m_vec + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        g1.append(m_vec[2])
    return tau_grid, np.array(g1)


def mollow_spectrum_grid(omega: float, delta: float, gamma: float,
                         rho_ss: np.ndarray, tau_max: float, n_tau: int,
                         n_om: int):
    """Resonance-fluorescence spectrum via the Wiener-Khinchin Fourier
    transform of g1(tau): S(omega) = Re int_0^inf g1 e^{i omega tau} dtau / pi.

    Weak drive: single Lorentzian at omega = delta (HWHM Gamma/2).
    Strong drive: the Mollow triplet (central + sidebands at
    +/- sqrt(Omega^2 - Gamma^2/4)).
    Returns (om_grid, s_grid).
    Fortran mirror: `mollow_spectrum_grid`.
    """
    tau_grid, g1 = first_order_correlation_2level(omega, delta, gamma,
                                                  rho_ss, tau_max, n_tau)
    d_tau = tau_grid[1] - tau_grid[0]
    om_lo = -2.0 * (omega + 5.0 * gamma)
    om_grid = om_lo + (om_lo * -2.0 + 2.0 * (omega + 5.0 * gamma)) * \
        np.arange(n_om) / (n_om - 1.0)
    phase = np.exp(1j * np.outer(om_grid, tau_grid))
    weights = np.ones(n_tau)
    weights[0] = weights[-1] = 0.5
    s_grid = np.real(phase @ (g1 * weights)) * d_tau / np.pi
    return om_grid, s_grid


def lindblad_supermatrix_general(h_mat, l_ops, gamma_rates):
    """General n-level Lindblad supermatrix (n^2 x n^2 complex).

    Row-major vectorization: vec(rho)[(i-1)*n+j] = rho[i,j].
    Fortran mirror: `lindblad_supermatrix_general`.
    """
    n = h_mat.shape[0]
    n_ops = len(gamma_rates)
    l_super = np.zeros((n * n, n * n), dtype=complex)
    for j in range(n):
        for i in range(n):
            e_mat = np.zeros((n, n), dtype=complex)
            e_mat[i, j] = 1.0
            l_m = -1j * (h_mat @ e_mat - e_mat @ h_mat)
            for ch in range(n_ops):
                if gamma_rates[ch] > 0:
                    l_op = l_ops[:, :, ch]
                    l_m += gamma_rates[ch] * (l_op @ e_mat @ l_op.conj().T
                                              - 0.5 * (l_op.conj().T @ l_op @ e_mat
                                                       + e_mat @ l_op.conj().T @ l_op))
            for k in range(n):
                for l in range(n):
                    l_super[k * n + l, i * n + j] = l_m[k, l]
    return l_super


def lindblad_steady_state_general(h_mat, l_ops, gamma_rates):
    """General n-level Lindblad steady state via the augmented system.

    Solves [[L, c], [r^T, 0]] [x; nu] = [0; 1] where c/r implement the
    Tr(rho) = 1 constraint (only diagonal positions in vec space).
    Returns rho_ss (n x n complex).
    Fortran mirror: `lindblad_steady_state_general`.
    """
    n = h_mat.shape[0]
    l_super = lindblad_supermatrix_general(h_mat, l_ops, gamma_rates)
    n_aug = n * n + 1
    aug = np.zeros((n_aug, n_aug), dtype=complex)
    aug[:n*n, :n*n] = l_super
    for i in range(n):
        aug[i * (n + 1), n_aug - 1] = 1.0
        aug[n_aug - 1, i * (n + 1)] = 1.0
    x = np.linalg.solve(aug, np.concatenate([np.zeros(n*n, dtype=complex), [1.0]]))
    rho_ss = x[:n * n].reshape(n, n)
    return rho_ss


def partial_trace_b(rho_ab: np.ndarray, dim_a: int, dim_b: int) -> np.ndarray:
    """Partial trace over subsystem B: rho_A = Tr_B(rho_AB).

    Row-major convention: rho_AB[(ia-1)*dim_b+ib, (ja-1)*dim_b+jb].
    Fortran mirror: `partial_trace_b`.
    """
    return np.trace(rho_ab.reshape(dim_a, dim_b, dim_a, dim_b),
                    axis1=1, axis2=3)


def entanglement_entropy(psi: np.ndarray, dim_a: int, dim_b: int) -> float:
    """Entanglement entropy of a pure bipartite state: -Tr(rho_A ln rho_A).

    Fortran mirror: `entanglement_entropy`.
    """
    rho_a = partial_trace_b(np.outer(psi, psi.conj()), dim_a, dim_b)
    evals = np.linalg.eigvalsh(rho_a)
    return float(-sum(p * np.log(p) for p in evals if p > 1e-14))
