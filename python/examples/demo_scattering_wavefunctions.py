#!/usr/bin/env python3
"""
GeneralModule Python Demonstration:
Solving and Comparing Continuous Scattering Energy Eigenfunctions
[Time-Independent (TI) vs Time-Dependent (TD)]
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add pygenmod to python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pygenmod import (
    calc_scattering_wavefunction_1d_cartesian,
    calc_scattering_wavefunction_ti,
    set_publication_style
)


def run_scattering_comparison():
    print("=" * 70)
    print("  GeneralModule Python: Continuous Scattering Wavefunction Demo   ")
    print("=" * 70)

    # --------------------------------------------------------------------------
    # 1. 1D 笛卡尔势垒散射: TI (逆向 Numerov) vs TD (含时波包谱投影)
    # --------------------------------------------------------------------------
    nx = 512
    x = np.linspace(-25.0, 25.0, nx)
    dx = x[1] - x[0]
    mass = 1.0
    hbar = 1.0
    v0 = 0.80
    v_pot = v0 * np.exp(-(x ** 2) / 2.0)
    energy = 0.50
    k0 = np.sqrt(2.0 * mass * energy) / hbar

    print(f">> [1/2] 1D Barrier: V0 = {v0:.2f} a.u., E = {energy:.2f} a.u., k = {k0:.4f} a.u.")

    # (A) 非含时定态求解 (TI)
    psi_ti, trans_ti, refl_ti = calc_scattering_wavefunction_1d_cartesian(
        x, v_pot, mass=mass, energy=energy, norm_type='energy'
    )
    print(f"    TI Method: Transmission T = {trans_ti:.5f}, Reflection R = {refl_ti:.5f}, Sum = {trans_ti + refl_ti:.5f}")

    # (B) 含时波包推进与谱投影求解 (TD)
    x0 = -12.0
    sigma_x = 1.5
    psi_wp = (1.0 / (2.0 * np.pi * sigma_x**2)**0.25) * np.exp(-(x - x0)**2 / (4.0 * sigma_x**2)) * np.exp(1j * k0 * x)

    p_grid = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx) * hbar
    exp_t = np.exp(-1j * (p_grid**2 / (2.0 * mass)) * 0.02)
    exp_v_half = np.exp(-0.5j * v_pot * 0.02)

    dt = 0.02
    nt = 1400
    psi_accum = np.zeros(nx, dtype=complex)
    psi_t = psi_wp.copy()

    for it in range(nt + 1):
        t = it * dt
        psi_accum += psi_t * np.exp(1j * energy * t / hbar) * dt
        # Split-Operator
        psi_t *= exp_v_half
        psi_t = np.fft.ifft(exp_t * np.fft.fft(psi_t))
        psi_t *= exp_v_half

    gk = (2.0 * sigma_x**2 / np.pi)**0.25 * np.exp(-sigma_x**2 * (k0 - k0)**2) * np.exp(-1j * k0 * x0)
    weight_e = np.sqrt(mass / (hbar**2 * k0))
    denom = 2.0 * np.pi * hbar * weight_e * gk
    psi_td = psi_accum / denom

    eval_mask = (x >= -5.0) & (x <= 15.0)
    rel_diff = np.abs(np.abs(psi_td[eval_mask]) - np.abs(psi_ti[eval_mask])) / np.max(np.abs(psi_ti))
    print(f"    TD Method: Max difference in [-5, 15] = {np.max(rel_diff):.4f}, Mean diff = {np.mean(rel_diff)*100:.2f}%")

    # --------------------------------------------------------------------------
    # 2. 3D 径向散射态: u_{l=0, E}(r) 与分波相移 delta_0
    # --------------------------------------------------------------------------
    r = np.linspace(0.01, 20.0, 1000)
    v_radial = -1.2 * np.exp(-(r - 2.5)**2)  # 吸引势阱
    e_rad = 0.35
    u_rad, delta_0 = calc_scattering_wavefunction_ti(
        r, v_radial, mass=mass, energy=e_rad, l=0, norm_type='energy'
    )
    print(f">> [2/2] 3D Radial Potential: E = {e_rad:.2f} a.u. -> s-wave Phase Shift delta_0 = {delta_0:.4f} rad ({np.degrees(delta_0):.2f} deg)")

    # --------------------------------------------------------------------------
    # 3. 绘制出版级高清对比图
    # --------------------------------------------------------------------------
    set_publication_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: 1D Cartesian Barrier (TI vs TD)
    ax1_twin = ax1.twinx()
    ax1_twin.fill_between(x, 0, v_pot, color='gray', alpha=0.18, label="Potential $V(x)$")
    ax1_twin.set_ylabel("Potential $V(x)$ (a.u.)", color='gray')
    ax1_twin.tick_params(axis='y', labelcolor='gray')
    ax1_twin.set_ylim(0, 1.2)

    ax1.plot(x, np.abs(psi_ti)**2, 'b-', lw=2.2, label=r"TI (Numerov): $|\psi_E^{\rm TI}|^2$")
    ax1.plot(x, np.abs(psi_td)**2, 'r--', lw=2.0, label=r"TD (Spectral Proj): $|\psi_E^{\rm TD}|^2$")
    ax1.plot(x, np.real(psi_ti), color='#0275D8', lw=1.0, alpha=0.5, label=r"${\rm Re}[\psi_E^{\rm TI}(x)]$")

    ax1.set_xlim(-15.0, 15.0)
    ax1.set_xlabel("Coordinate $x$ (a.u.)")
    ax1.set_ylabel(r"Probability Density $|\psi_E(x)|^2$ (a.u.)")
    ax1.set_title(f"1D Cartesian Scattering ($E={energy:.2f}$, $T={trans_ti:.3f}$)", fontweight='bold')
    ax1.legend(loc='upper left', frameon=True)
    ax1.grid(True, alpha=0.25)

    # Panel 2: 3D Radial State
    ax2_twin = ax2.twinx()
    ax2_twin.plot(r, v_radial, 'k:', lw=1.5, label="Well $V(r)$")
    ax2_twin.set_ylabel("Potential $V(r)$ (a.u.)", color='k')
    ax2_twin.tick_params(axis='y', labelcolor='k')

    ax2.plot(r, u_rad, color='#D9534F', lw=2.2, label=rf"$u_{{0, E}}(r)$ ($\delta_0 = {delta_0:.3f}$ rad)")
    ax2.axhline(0, color='gray', ls='--', lw=0.8)
    ax2.set_xlim(0, 18.0)
    ax2.set_xlabel("Radial Coordinate $r$ (a.u.)")
    ax2.set_ylabel(r"Radial Scattering Wavefunction $u_E(r)$ (a.u.)")
    ax2.set_title(f"3D Radial Partial Wave ($l=0, E={e_rad:.2f}$ a.u.)", fontweight='bold')
    ax2.legend(loc='lower right', frameon=True)
    ax2.grid(True, alpha=0.25)

    plt.tight_layout()
    out_file = "result_scattering_wavefunctions_comparison.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"\n>> High-resolution comparison plot saved to: {out_file}")
    print("=" * 70)


if __name__ == "__main__":
    run_scattering_comparison()
