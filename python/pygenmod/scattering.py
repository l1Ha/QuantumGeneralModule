"""
pygenmod: Scattering and Ultracold Collision Module
Provides non-time-dependent and time-dependent quantum scattering tools:
- Scattering length (Numerov & Johnson Log-Derivative)
- Effective range expansion (ERE)
- Van der Waals mean scattering length & Gribakin-Flambaum formula
- Time-dependent transmission & S-matrix analysis
- Publication-quality visualization of scattering wave functions and asymptotic extrapolation
"""
import numpy as np
import matplotlib.pyplot as plt


def van_der_waals_mean_length(mass: float, c6_au: float) -> float:
    """Compute mean van der Waals scattering length a_bar."""
    return 0.47798881258618 * ((2.0 * mass * c6_au) ** 0.25)


def gribakin_flambaum_length(mass: float, c6_au: float, phase_phi: float) -> float:
    """Compute s-wave scattering length via Gribakin-Flambaum formula: a_s = a_bar * [1 - tan(Phi - pi/8)]."""
    a_bar = van_der_waals_mean_length(mass, c6_au)
    return a_bar * (1.0 - np.tan(phase_phi - np.pi / 8.0))


def square_well_scattering_length_exact(r_well: float, v0_well: float, mass: float) -> float:
    """Exact analytical scattering length for attractive square well: a_s = R * [1 - tan(kappa*R)/(kappa*R)]."""
    kappa = np.sqrt(2.0 * mass * v0_well)
    kr = kappa * r_well
    return r_well * (1.0 - np.tan(kr) / kr)


def calc_scattering_length_numerov(r_grid: np.ndarray, v_pot: np.ndarray, mass: float) -> tuple[float, np.ndarray]:
    """Compute s-wave scattering length a_s and zero-energy wave function u(r) via Numerov algorithm."""
    n_pts = len(r_grid)
    dr = r_grid[1] - r_grid[0]
    dr2_12 = (dr * dr) / 12.0

    u = np.zeros(n_pts)
    u[0] = 0.0
    u[1] = 1e-7

    for i in range(1, n_pts - 1):
        q_prev = -2.0 * mass * v_pot[i - 1]
        q_curr = -2.0 * mass * v_pot[i]
        q_next = -2.0 * mass * v_pot[i + 1]

        c_prev = 1.0 + dr2_12 * q_prev
        c_curr = 2.0 * (1.0 - 5.0 * dr2_12 * q_curr)
        c_next = 1.0 + dr2_12 * q_next

        u[i + 1] = (c_curr * u[i] - c_prev * u[i - 1]) / c_next
        if abs(u[i + 1]) > 1e15:
            u[:i + 2] *= 1e-10

    # Asymptotic derivative: (3 u_N - 4 u_{N-1} + u_{N-2}) / (2 dr)
    du = (3.0 * u[-1] - 4.0 * u[-2] + u[-3]) / (2.0 * dr)
    if abs(du) > 1e-14:
        a_s = r_grid[-1] - u[-1] / du
    else:
        a_s = 1e30
    return a_s, u


def plot_scattering_length_wavefunction(r_grid: np.ndarray, v_pot: np.ndarray,
                                       u_wf: np.ndarray, a_s: float,
                                       title: str = "Zero-Energy Wave Function and Scattering Length",
                                       filename: str = "result_scattering_length.png"):
    """
    Generate publication-quality plot demonstrating the zero-energy radial wave function u(r)
    and its asymptotic linear extrapolation C*(r - a_s) crossing the r-axis at the scattering length a_s.
    """
    plt.rcParams.update({
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'lines.linewidth': 2.0,
        'figure.dpi': 300
    })

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                   gridspec_kw={'height_ratios': [1, 1.6]})

    # Top panel: Interaction Potential V(r)
    ax1.plot(r_grid, v_pot, color='#D9534F', lw=2.2, label=r'Potential $V(r)$')
    ax1.axhline(0.0, color='gray', linestyle='--', alpha=0.6)
    ax1.set_ylabel(r'$V(r)$ (a.u.)')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    ax1.set_title(title, fontweight='bold')

    # Bottom panel: Zero-energy wave function u(r) and tangent/asymptote
    # Normalize u for visualization
    u_norm = u_wf / np.max(np.abs(u_wf))
    slope = (u_norm[-1] - u_norm[-2]) / (r_grid[-1] - r_grid[-2])
    asymptote = slope * (r_grid - a_s)

    ax2.plot(r_grid, u_norm, color='#0275D8', lw=2.5, label=r'Zero-Energy Wave Function $u(r)$')
    ax2.plot(r_grid, asymptote, color='#F0AD4E', linestyle='--', lw=2.0,
             label=rf'Asymptote $\propto (r - a_s)$ ($a_s = {a_s:.3f}$ a.u.)')
    ax2.axvline(a_s, color='#5CB85C', linestyle=':', lw=2.0, label=rf'Scattering Length $a_s = {a_s:.3f}$')
    ax2.axhline(0.0, color='black', linestyle='-', lw=0.8, alpha=0.5)

    ax2.set_xlabel(r'Interatomic Distance $r$ (Bohr / a.u.)')
    ax2.set_ylabel(r'$u(r)$ (arb. units)')
    ax2.set_ylim(-1.5, 1.5)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[Scattering] High-resolution figure saved to: {filename}")


def calc_differential_cross_section(energy: float, mass: float, delta_arr: np.ndarray,
                                    theta_grid: np.ndarray) -> np.ndarray:
    """Compute differential scattering cross section dsigma/dOmega(theta)."""
    k = np.sqrt(2.0 * mass * max(1e-14, energy))
    l_max = len(delta_arr) - 1
    f_theta = np.zeros(len(theta_grid), dtype=complex)
    for l in range(l_max + 1):
        pl = np.polynomial.legendre.Legendre.basis(l)(np.cos(theta_grid))
        f_theta += (2 * l + 1) * np.exp(1j * delta_arr[l]) * np.sin(delta_arr[l]) * pl
    f_theta /= k
    return np.abs(f_theta) ** 2


def calc_differential_cross_section_identical(energy: float, mass: float, delta_arr: np.ndarray,
                                              theta_grid: np.ndarray,
                                              particle_stat: str = 'distinguishable') -> np.ndarray:
    """
    Compute differential cross section considering quantum statistics of identical particles:
    - 'distinguishable': |f(theta)|^2
    - 'boson': |f(theta) + f(pi - theta)|^2
    - 'fermion': |f(theta) - f(pi - theta)|^2
    - 'fermion_unpolarized': 0.25*|f(th)+f(pi-th)|^2 + 0.75*|f(th)-f(pi-th)|^2
    """
    k = np.sqrt(2.0 * mass * max(1e-14, energy))
    l_max = len(delta_arr) - 1
    f_th = np.zeros(len(theta_grid), dtype=complex)
    f_pi = np.zeros(len(theta_grid), dtype=complex)
    for l in range(l_max + 1):
        pl = np.polynomial.legendre.Legendre.basis(l)(np.cos(theta_grid))
        pl_pi = ((-1.0) ** l) * pl
        term = (2 * l + 1) * np.exp(1j * delta_arr[l]) * np.sin(delta_arr[l])
        f_th += term * pl
        f_pi += term * pl_pi
    f_th /= k
    f_pi /= k

    if particle_stat == 'boson':
        return np.abs(f_th + f_pi) ** 2
    elif particle_stat == 'fermion':
        return np.abs(f_th - f_pi) ** 2
    elif particle_stat == 'fermion_unpolarized':
        return 0.25 * np.abs(f_th + f_pi) ** 2 + 0.75 * np.abs(f_th - f_pi) ** 2
    else:
        return np.abs(f_th) ** 2


def plot_differential_cross_sections(theta_grid: np.ndarray, dsigma_dict: dict,
                                     filename: str = "result_differential_cross_section.png"):
    """
    Plot differential cross sections in both polar and Cartesian coordinates.
    dsigma_dict maps label to dsigma/dOmega array.
    """
    fig = plt.figure(figsize=(11, 5))
    ax_cart = fig.add_subplot(1, 2, 1)
    ax_polar = fig.add_subplot(1, 2, 2, projection='polar')

    colors = ['#0275D8', '#D9534F', '#5CB85C', '#F0AD4E']
    theta_deg = np.degrees(theta_grid)

    for (label, ds), col in zip(dsigma_dict.items(), colors):
        ax_cart.plot(theta_deg, ds, lw=2.2, label=label, color=col)
        ax_polar.plot(theta_grid, ds, lw=2.0, label=label, color=col)

    ax_cart.set_xlabel(r'Scattering Angle $\theta$ (deg)')
    ax_cart.set_ylabel(r'$d\sigma/d\Omega$ (a.u.)')
    ax_cart.set_title("Differential Cross Section (Cartesian)", fontweight='bold')
    ax_cart.grid(True, alpha=0.3)
    ax_cart.legend(loc='upper right')

    ax_polar.set_title("Angular Distribution (Polar)", fontweight='bold', va='bottom')
    ax_polar.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[Scattering] Differential cross section plot saved to: {filename}")


def calc_multichannel_close_coupling(r_grid: np.ndarray, v_mat: np.ndarray,
                                     mass: float, total_energy: float,
                                     thresholds: np.ndarray, l_channels: np.ndarray = None) -> dict:
    """
    Multichannel close-coupling solver using Johnson's matrix log-derivative method.
    v_mat has shape (n_chan, n_chan, n_pts).
    Returns a dictionary with:
    - 'n_open': number of open channels
    - 'n_closed': number of closed channels
    - 'open_channels': list of open channel indices
    - 'closed_channels': list of closed channel indices
    - 'k_open': wavenumbers
    - 'k_matrix': K_oo reaction matrix
    - 's_matrix': S_oo unitary scattering matrix
    - 'prob_matrix': P(i->j) = |S_ij|^2 transition probabilities
    - 'eigenphase_sum': sum of eigenphases
    """
    n_chan = len(thresholds)
    n_pts = len(r_grid)
    if l_channels is None:
        l_channels = np.zeros(n_chan, dtype=int)
    dr = r_grid[1] - r_grid[0]
    dr2_12 = (dr * dr) / 12.0

    # Classify open / closed
    e_kin = total_energy - thresholds
    open_idx = np.where(e_kin > 1e-12)[0]
    closed_idx = np.where(e_kin <= 1e-12)[0]
    n_open = len(open_idx)
    n_closed = len(closed_idx)

    if n_open == 0:
        return {'n_open': 0, 'n_closed': n_closed}

    k_open = np.sqrt(2.0 * mass * e_kin[open_idx])
    kappa_closed = np.sqrt(2.0 * mass * np.maximum(0.0, -e_kin[closed_idx]))

    # Propagate ratio matrix R
    def build_w(step):
        r_c = r_grid[step]
        w = 2.0 * mass * v_mat[:, :, step].copy()
        for i in range(n_chan):
            w[i, i] += -2.0 * mass * e_kin[i] + l_channels[i] * (l_channels[i] + 1) / (r_c * r_c)
        return w

    w1 = build_w(0)
    q1 = np.eye(n_chan) - dr2_12 * w1
    q1_inv = np.linalg.inv(q1)
    m1 = 12.0 * q1_inv - 10.0 * np.eye(n_chan)
    r_curr = m1.copy()
    q_prev = q1.copy()

    for step in range(1, n_pts - 1):
        w_step = build_w(step)
        q_step = np.eye(n_chan) - dr2_12 * w_step
        q_inv = np.linalg.inv(q_step)
        m_step = 12.0 * q_inv - 10.0 * np.eye(n_chan)

        r_inv = np.linalg.inv(r_curr)
        r_next = m_step - r_inv
        r_next = 0.5 * (r_next + r_next.T)
        if step == n_pts - 2:
            q_curr = q_step.copy()
        r_curr = r_next

    # Outer boundary log-derivative
    w_end = build_w(n_pts - 1)
    q_end = np.eye(n_chan) - dr2_12 * w_end
    r_inv = np.linalg.inv(r_curr)
    p1 = np.linalg.inv(q_curr) @ r_inv @ q_end
    p2 = p1 @ p1
    y_mat = (3.0 * np.eye(n_chan) - 4.0 * p1 + p2) / (2.0 * dr)
    y_mat = 0.5 * (y_mat + y_mat.T)

    # Open-closed Schur complement
    if n_closed > 0:
        y_oo = y_mat[np.ix_(open_idx, open_idx)]
        y_oc = y_mat[np.ix_(open_idx, closed_idx)]
        y_co = y_mat[np.ix_(closed_idx, open_idx)]
        y_cc = y_mat[np.ix_(closed_idx, closed_idx)]
        a_cc = y_cc + np.diag(kappa_closed)
        a_cc_inv = np.linalg.inv(a_cc)
        y_eff = y_oo - y_oc @ a_cc_inv @ y_co
        y_eff = 0.5 * (y_eff + y_eff.T)
    else:
        y_eff = y_mat.copy()

    # Riccati matching
    r_match = r_grid[-1]
    j_mat = np.zeros((n_open, n_open))
    n_mat = np.zeros((n_open, n_open))
    dj_mat = np.zeros((n_open, n_open))
    dn_mat = np.zeros((n_open, n_open))

    for i in range(n_open):
        k_i = k_open[i]
        l_i = l_channels[open_idx[i]]
        x = k_i * r_match
        if l_i == 0:
            jl = np.sin(x)
            nl = -np.cos(x)
            djl = np.cos(x)
            dnl = np.sin(x)
        elif l_i == 1:
            jl = np.sin(x) / x - np.cos(x)
            nl = -np.cos(x) / x - np.sin(x)
            djl = np.cos(x) / x - np.sin(x) / (x * x) + np.sin(x)
            dnl = np.sin(x) / x + np.cos(x) / (x * x) - np.cos(x)
        else:
            jl = np.sin(x - l_i * np.pi / 2.0)
            nl = -np.cos(x - l_i * np.pi / 2.0)
            djl = np.cos(x - l_i * np.pi / 2.0)
            dnl = np.sin(x - l_i * np.pi / 2.0)
        j_mat[i, i] = jl / np.sqrt(k_i)
        n_mat[i, i] = nl / np.sqrt(k_i)
        dj_mat[i, i] = djl * np.sqrt(k_i)
        dn_mat[i, i] = dnl * np.sqrt(k_i)

    mj = dj_mat - y_eff @ j_mat
    mn = dn_mat - y_eff @ n_mat
    k_mat = np.linalg.inv(mn) @ mj
    k_mat = 0.5 * (k_mat + k_mat.T)

    # S-matrix via Cayley transform
    eye_c = np.eye(n_open, dtype=complex)
    ik = 1j * k_mat
    s_mat = (eye_c + ik) @ np.linalg.inv(eye_c - ik)
    prob_mat = np.abs(s_mat) ** 2
    det_s = np.linalg.det(s_mat)
    eigenphase_sum = 0.5 * np.angle(det_s)

    return {
        'n_open': n_open,
        'n_closed': n_closed,
        'open_channels': open_idx,
        'closed_channels': closed_idx,
        'k_open': k_open,
        'k_matrix': k_mat,
        's_matrix': s_mat,
        'prob_matrix': prob_mat,
        'eigenphase_sum': eigenphase_sum
    }


def plot_multichannel_smatrix(prob_mat: np.ndarray, channel_labels: list = None,
                              filename: str = "result_multichannel_smatrix.png"):
    """
    Plot heatmap of multichannel transition probabilities P(i->j) = |S_ij|^2.
    """
    n = prob_mat.shape[0]
    if channel_labels is None:
        channel_labels = [f"Ch {i+1}" for i in range(n)]

    plt.figure(figsize=(6, 5))
    im = plt.imshow(prob_mat, cmap='viridis', vmin=0.0, vmax=1.0, origin='upper')
    cbar = plt.colorbar(im)
    cbar.set_label(r'Transition Probability $|S_{ij}|^2$')

    plt.xticks(range(n), channel_labels)
    plt.yticks(range(n), channel_labels)
    plt.xlabel("Entrance Channel $i$")
    plt.ylabel("Exit Channel $j$")
    plt.title("Multichannel S-Matrix Transition Probabilities", fontweight='bold')

    for i in range(n):
        for j in range(n):
            color = "white" if prob_mat[i, j] < 0.5 else "black"
            plt.text(j, i, f"{prob_mat[i, j]:.3f}", ha="center", va="center", color=color, fontweight='bold')

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[Scattering] Multichannel S-matrix plot saved to: {filename}")


def plot_feshbach_resonance(energy_grid: np.ndarray, scattering_lengths: np.ndarray,
                            eigenphase_sums: np.ndarray = None,
                            filename: str = "result_feshbach_resonance.png"):
    """
    Plot Feshbach resonance profile showing scattering length pole a_s(E) and eigenphase shift.
    """
    fig, ax1 = plt.subplots(figsize=(7, 4.5))

    color_as = '#0275D8'
    ax1.set_xlabel('Collision Energy $E$ (a.u.)')
    ax1.set_ylabel(r's-Wave Scattering Length $a_s$ (a.u.)', color=color_as)
    ax1.plot(energy_grid, scattering_lengths, lw=2.2, color=color_as, label=r'$a_s(E)$')
    ax1.tick_params(axis='y', labelcolor=color_as)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0, color='gray', ls='--', lw=0.8)

    if eigenphase_sums is not None:
        ax2 = ax1.twinx()
        color_delta = '#D9534F'
        ax2.set_ylabel(r'Eigenphase Sum $\delta_{\rm sum}$ (rad)', color=color_delta)
        ax2.plot(energy_grid, eigenphase_sums, lw=2.0, ls='-.', color=color_delta, label=r'$\delta_{\rm sum}$')
        ax2.tick_params(axis='y', labelcolor=color_delta)

    plt.title("Multichannel Feshbach Resonance Profile", fontweight='bold')
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[Scattering] Feshbach resonance plot saved to: {filename}")


def calc_scattering_wavefunction_ti(r_grid: np.ndarray, v_pot: np.ndarray,
                                    mass: float, energy: float, l: int = 0,
                                    norm_type: str = 'energy') -> tuple:
    """
    Time-independent scattering radial eigenfunction u_{l, E}(r).
    Returns (u_wf, phase_shift).
    norm_type: 'energy' (delta(E-E') normalized, amplitude sqrt(2*mu/(pi*hbar^2*k))),
               'momentum' (delta(k-k') normalized, amplitude sqrt(2/pi)),
               'unit' (asymptotic amplitude 1.0).
    """
    n_pts = len(r_grid)
    dr = r_grid[1] - r_grid[0]
    dr2_12 = (dr * dr) / 12.0
    k = np.sqrt(2.0 * mass * max(1e-14, energy))

    u_wf = np.zeros(n_pts)
    u_wf[0] = 0.0
    u_wf[1] = (dr ** (l + 1)) * 1e-5

    q = 2.0 * mass * (energy - v_pot) - l * (l + 1) / (r_grid ** 2)

    for i in range(1, n_pts - 1):
        c_prev = 1.0 + dr2_12 * q[i - 1]
        c_curr = 2.0 * (1.0 - 5.0 * dr2_12 * q[i])
        c_next = 1.0 + dr2_12 * q[i + 1]
        u_wf[i + 1] = (c_curr * u_wf[i] - c_prev * u_wf[i - 1]) / c_next
        if abs(u_wf[i + 1]) > 1e20:
            u_wf[:i + 2] *= 1e-15

    # Log-derivative and asymptotic matching
    du = (3.0 * u_wf[-1] - 4.0 * u_wf[-2] + u_wf[-3]) / (2.0 * dr)
    y_logder = du / u_wf[-1]

    x = k * r_grid[-1]
    if l == 0:
        jl, nl = np.sin(x), -np.cos(x)
        djl, dnl = np.cos(x), np.sin(x)
    else:
        jl = np.sin(x - l * np.pi / 2.0)
        nl = -np.cos(x - l * np.pi / 2.0)
        djl = np.cos(x - l * np.pi / 2.0)
        dnl = np.sin(x - l * np.pi / 2.0)

    delta = np.arctan2(k * djl - y_logder * jl, k * dnl - y_logder * nl)

    asymp_amp = np.sqrt(u_wf[-1]**2 + (du / k)**2)
    if norm_type == 'energy':
        norm_target = np.sqrt(2.0 * mass / (np.pi * k))
    elif norm_type == 'momentum':
        norm_target = np.sqrt(2.0 / np.pi)
    else:
        norm_target = 1.0

    target_asymp = np.cos(delta) * jl - np.sin(delta) * nl
    scale = norm_target / max(1e-30, asymp_amp)
    if u_wf[-1] * target_asymp < 0:
        scale = -scale
    u_wf *= scale

    return u_wf, delta


def plot_scattering_wavefunction(r_grid: np.ndarray, v_pot: np.ndarray,
                                 u_wf: np.ndarray, energy: float,
                                 phase_shift: float = None,
                                 filename: str = "result_scattering_wavefunction.png"):
    """
    Plot the continuous scattering energy eigenfunction u_E(r) alongside the potential.
    """
    fig, ax1 = plt.subplots(figsize=(8, 5))

    color_wf = '#0275D8'
    color_pot = '#D9534F'

    ax1.set_xlabel('Radial Coordinate $r$ (a.u.)')
    ax1.set_ylabel(r'Scattering Wavefunction $u_E(r)$ (a.u.)', color=color_wf)
    title_str = f"Continuous Scattering Energy Eigenstate ($E = {energy:.4f}$ a.u."
    if phase_shift is not None:
        title_str += f", $\\delta = {phase_shift:.4f}$ rad)"
    else:
        title_str += ")"
    ax1.plot(r_grid, u_wf, lw=2.2, color=color_wf, label=r'$u_E(r)$')
    ax1.tick_params(axis='y', labelcolor=color_wf)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0, color='gray', ls='--', lw=0.8)

    ax2 = ax1.twinx()
    ax2.set_ylabel(r'Potential $V(r)$ (a.u.)', color=color_pot)
    ax2.plot(r_grid, v_pot, lw=1.8, ls='--', color=color_pot, label=r'$V(r)$')
    ax2.tick_params(axis='y', labelcolor=color_pot)

    plt.title(title_str, fontweight='bold')
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[Scattering] Scattering wavefunction plot saved to: {filename}")


def calc_scattering_wavefunction_1d_cartesian(x_grid: np.ndarray, v_pot: np.ndarray,
                                              mass: float, energy: float,
                                              norm_type: str = 'energy') -> tuple:
    """
    1D Cartesian time-independent scattering energy eigenfunction psi_E(x).
    Integrates backward from transmission boundary to determine incident/reflected amplitudes.
    Returns (psi_wf, trans_prob, refl_prob).
    norm_type: 'energy' (delta(E-E') normalized, incident amplitude 1/sqrt(2*pi) * sqrt(m/(hbar^2*k))),
               'momentum' (delta(k-k') normalized, incident amplitude 1/sqrt(2*pi)),
               'unit' (unit incident amplitude 1.0).
    """
    n_pts = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    dx2_12 = (dx * dx) / 12.0
    k = np.sqrt(2.0 * mass * max(1e-14, energy))

    q = 2.0 * mass * (energy - v_pot)
    psi = np.zeros(n_pts, dtype=complex)
    psi[-1] = np.exp(1j * k * x_grid[-1])
    psi[-2] = np.exp(1j * k * x_grid[-2])

    for i in range(n_pts - 2, 0, -1):
        c_curr = 2.0 * (1.0 - 5.0 * dx2_12 * q[i]) * psi[i]
        c_next = (1.0 + dx2_12 * q[i + 1]) * psi[i + 1]
        c_prev = 1.0 + dx2_12 * q[i - 1]
        psi[i - 1] = (c_curr - c_next) / c_prev
        if abs(psi[i - 1]) > 1e20:
            psi[i - 1:] *= 1e-15

    d_psi_left = (-3.0 * psi[0] + 4.0 * psi[1] - psi[2]) / (2.0 * dx)
    a_inc = 0.5 * (psi[0] - 1j * d_psi_left / k) * np.exp(-1j * k * x_grid[0])
    b_ref = 0.5 * (psi[0] + 1j * d_psi_left / k) * np.exp(1j * k * x_grid[0])

    t_prob = 1.0 / (abs(a_inc)**2)
    r_prob = (abs(b_ref)**2) / (abs(a_inc)**2)

    if norm_type == 'energy':
        norm_factor = (1.0 / np.sqrt(2.0 * np.pi)) * np.sqrt(mass / k)
    elif norm_type == 'momentum':
        norm_factor = 1.0 / np.sqrt(2.0 * np.pi)
    else:
        norm_factor = 1.0

    psi_wf = (psi / a_inc) * norm_factor
    return psi_wf, t_prob, r_prob
