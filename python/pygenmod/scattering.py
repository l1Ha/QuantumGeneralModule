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

