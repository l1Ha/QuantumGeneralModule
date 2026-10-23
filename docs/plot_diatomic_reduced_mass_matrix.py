#!/usr/bin/env python3
"""
GeneralModule: Diatomic Reduced Mass Matrix and Global Periodic Trends.

Computes and visualizes the reduced mass matrix:
    mu = (m1 * m2) / (m1 + m2)
across:
1. Representative Diatomic Matrix: 16 key isotopic components from ultra-light (H, D)
   to halogens and heavy alkali metals (Li, Na, K, Rb, Cs, C, N, O, F, Cl, Br, I).
2. Global Periodic Landscape: Elements Z = 1 to 36 (H to Kr) showing the hydride
   saturation boundary (mu <= 1 amu) and the homonuclear ridge (mu = m / 2).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

def main():
    # 1. Representative atomic isotopes with IUPAC/NIST standard atomic masses (amu)
    atoms = [
        (r"$^{1}\mathrm{H}$", 1.0078),
        (r"$^{2}\mathrm{H(D)}$", 2.0141),
        (r"$^{4}\mathrm{He}$", 4.0026),
        (r"$^{6}\mathrm{Li}$", 6.0151),
        (r"$^{7}\mathrm{Li}$", 7.0160),
        (r"$^{12}\mathrm{C}$", 12.0000),
        (r"$^{14}\mathrm{N}$", 14.0031),
        (r"$^{16}\mathrm{O}$", 15.9949),
        (r"$^{19}\mathrm{F}$", 18.9984),
        (r"$^{23}\mathrm{Na}$", 22.9898),
        (r"$^{35}\mathrm{Cl}$", 34.9689),
        (r"$^{40}\mathrm{K}$", 39.9640),
        (r"$^{79}\mathrm{Br}$", 78.9183),
        (r"$^{87}\mathrm{Rb}$", 86.9092),
        (r"$^{127}\mathrm{I}$", 126.9045),
        (r"$^{133}\mathrm{Cs}$", 132.9055),
    ]

    names = [a[0] for a in atoms]
    masses = np.array([a[1] for a in atoms])
    n = len(atoms)

    # 2D grid of reduced mass: mu = m1 * m2 / (m1 + m2)
    M1, M2 = np.meshgrid(masses, masses)
    mu_matrix = (M1 * M2) / (M1 + M2)

    # 2. Periodic table standard atomic weights for Z = 1 to 36
    periodic_masses = [
        1.008, 4.003, 6.94, 9.012, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, # 1-10 (H - Ne)
        22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.95, 39.098, 40.078, # 11-20 (Na - Ca)
        44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, # 21-30 (Sc - Zn)
        69.723, 72.630, 74.922, 78.971, 79.904, 83.798 # 31-36 (Ga - Kr)
    ]
    p_masses = np.array(periodic_masses)
    P1, P2 = np.meshgrid(p_masses, p_masses)
    mu_periodic = (P1 * P2) / (P1 + P2)

    plt.rcParams.update({
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'axes.edgecolor': '#333333',
        'axes.linewidth': 1.0,
    })

    fig = plt.figure(figsize=(19, 8.8), dpi=250)

    # ----------------- Panel (a): Representative Precision Grid -----------------
    ax1 = fig.add_subplot(1, 2, 1)
    norm1 = LogNorm(vmin=0.45, vmax=70.0)
    cmap1 = plt.cm.viridis
    im1 = ax1.imshow(mu_matrix, cmap=cmap1, norm=norm1, origin='lower')

    # Automatic W3C luminance contrast for text labels
    for i in range(n):
        for j in range(n):
            val = mu_matrix[i, j]
            rgba = cmap1(norm1(val))
            luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            text_color = "#111111" if luminance > 0.52 else "#ffffff"
            txt = f"{val:.2f}" if val < 10.0 else f"{val:.1f}"
            ax1.text(j, i, txt, ha="center", va="center", color=text_color,
                     fontsize=7.2, fontweight='bold' if luminance > 0.52 else 'normal')

    ax1.plot([-0.5, n - 0.5], [-0.5, n - 0.5], color='crimson', linestyle='--',
             linewidth=1.5, alpha=0.85, label='Homonuclear Line ($m_1=m_2$)')

    ax1.set_xticks(range(n))
    ax1.set_yticks(range(n))
    ax1.set_xticklabels(names, rotation=45, ha='right', fontsize=9.5)
    ax1.set_yticklabels(names, fontsize=9.5)
    ax1.set_xlabel("Atom 1 Component ($m_1$)", fontsize=11, fontweight='bold', labelpad=6)
    ax1.set_ylabel("Atom 2 Component ($m_2$)", fontsize=11, fontweight='bold', labelpad=6)
    ax1.set_title(r"(a) Representative Diatomic Reduced Mass $\mu = \frac{m_1 m_2}{m_1 + m_2}$ (amu)",
                  fontsize=12, fontweight='bold', pad=12)
    ax1.legend(loc='upper left', frameon=True, facecolor='white', edgecolor='gray', fontsize=8.5)

    cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    cbar1.set_label(r"Reduced Mass $\mu$ (amu, Logarithmic Color Scale)", fontsize=10, fontweight='bold')

    # ----------------- Panel (b): Periodic Table (Z = 1 to 36) -----------------
    ax2 = fig.add_subplot(1, 2, 2)
    cmap2 = plt.cm.plasma
    im2 = ax2.imshow(mu_periodic, cmap=cmap2, origin='lower', extent=[0.5, 36.5, 0.5, 36.5])

    # Contours
    levels = [0.98, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
    cs = ax2.contour(np.arange(1, 37), np.arange(1, 37), mu_periodic, levels=levels,
                     colors='white', linewidths=0.9, alpha=0.85)
    ax2.clabel(cs, inline=True, fontsize=8.5, fmt='%.1f')

    # Homonuclear diagonal ridge
    ax2.plot([1, 36], [1, 36], color='cyan', linestyle='--', linewidth=1.8,
             label=r'Homonuclear Ridge ($\mu = m/2$)')

    # Hydride saturation band
    ax2.axhline(1, color='lime', linestyle=':', linewidth=1.4, alpha=0.9)
    ax2.axvline(1, color='lime', linestyle=':', linewidth=1.4, alpha=0.9)
    ax2.text(18, 2.5, r"Hydride Saturation Band: $\mu \approx m_{\mathrm{H}} \leq 1.0\ \mathrm{amu}$",
             color='yellow', fontsize=9.5, ha='center', fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.25', facecolor='#222222', edgecolor='yellow', alpha=0.85))

    ax2.set_xlim(0.5, 36.5)
    ax2.set_ylim(0.5, 36.5)
    ax2.set_xlabel(r"Atomic Number of Atom 1 ($Z_1$)", fontsize=11, fontweight='bold', labelpad=6)
    ax2.set_ylabel(r"Atomic Number of Atom 2 ($Z_2$)", fontsize=11, fontweight='bold', labelpad=6)
    ax2.set_title(r"(b) Global Periodic Landscape of Reduced Mass ($Z = 1$ to $36$)",
                  fontsize=12, fontweight='bold', pad=12)
    ax2.legend(loc='lower right', frameon=True, facecolor='#222222', edgecolor='white',
               labelcolor='white', fontsize=9)

    cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label(r"Reduced Mass $\mu$ (amu, Linear Color Scale)", fontsize=10, fontweight='bold')

    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, "fig_diatomic_reduced_mass_matrix.png")
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()
    print(f"Generated successfully: {out_path}")

if __name__ == "__main__":
    main()
