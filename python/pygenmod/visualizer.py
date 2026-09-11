"""
Scientific visualization utilities for wavefunctions, pulses, and quantum dynamics.
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional


def set_publication_style():
    """Configure matplotlib for clean publication-quality plots."""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': 150,
        'lines.linewidth': 1.8,
        'grid.alpha': 0.4
    })


def plot_wavefunctions(x: np.ndarray, v_pot: np.ndarray,
                       eig_vals: np.ndarray, wavefunctions: np.ndarray,
                       n_states: int = 5, scale: float = 0.015,
                       title: str = "Bound State Energy Levels & Wavefunctions",
                       filename: Optional[str] = None):
    """
    Plot 1D potential energy curve with overlaid bound state energy levels and wavefunctions.
    """
    set_publication_style()
    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Plot potential curve
    ax.plot(x, v_pot, color='black', linewidth=2.0, label='Potential $V(x)$')

    n_plot = min(n_states, len(eig_vals))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, n_plot))

    for v in range(n_plot):
        e_val = eig_vals[v]
        psi = wavefunctions[:, v]
        # Baseline at energy eigenvalue
        ax.axhline(e_val, color=colors[v], linestyle='--', alpha=0.5)
        # Shifted wavefunction: E_v + scale * psi
        ax.plot(x, e_val + scale * psi, color=colors[v], label=f'$v={v}$ ($E={e_val:.4f}$ a.u.)')
        ax.fill_between(x, e_val, e_val + scale * psi, color=colors[v], alpha=0.15)

    # Automatically set y limits
    e_max = eig_vals[n_plot - 1]
    ax.set_ylim(min(v_pot) - 0.05 * abs(min(v_pot) + 1e-6), e_max * 1.3)
    ax.set_xlabel("Coordinate $R$ / $x$ (a.u.)")
    ax.set_ylabel("Energy (Hartree)")
    ax.set_title(title, fontweight='bold')
    ax.legend(loc='upper right', frameon=True)
    fig.tight_layout()

    if filename:
        fig.savefig(filename, dpi=300)
        print(f"[Visualizer] Saved figure to: {filename}")
    plt.close(fig)


def plot_pulses(t: np.ndarray, fields: List[np.ndarray], labels: List[str],
                title: str = "Laser Pulse Electric Field Syntheses",
                filename: Optional[str] = None):
    """Plot multiple pulse electric fields in time domain."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(9, 4.5))

    for field, label in zip(fields, labels):
        ax.plot(t, field, label=label)

    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Electric Field (a.u.)")
    ax.set_title(title, fontweight='bold')
    ax.legend(loc='upper right', frameon=True)
    fig.tight_layout()

    if filename:
        fig.savefig(filename, dpi=300)
        print(f"[Visualizer] Saved pulse plot to: {filename}")
    plt.close(fig)


def plot_alignment_dynamics(t: np.ndarray, traces: List[np.ndarray], labels: List[str],
                            title: str = "Molecular Alignment Dynamics $\\langle\\cos^2\\theta\\rangle(t)$",
                            filename: Optional[str] = None):
    """Plot molecular orientation or alignment dynamics."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    for trace, label in zip(traces, labels):
        ax.plot(t, trace, label=label)

    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Expectation Value")
    ax.set_title(title, fontweight='bold')
    ax.legend(loc='upper right', frameon=True)
    fig.tight_layout()

    if filename:
        fig.savefig(filename, dpi=300)
        print(f"[Visualizer] Saved alignment plot to: {filename}")
    plt.close(fig)
