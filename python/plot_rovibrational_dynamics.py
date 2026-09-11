#!/usr/bin/env python3
"""
Publication-Quality Plotting Pipeline for Laser Rovibrational Quantum Dynamics
Reads: test_rovibrational_dynamics.dat
Generates: result_rovibrational_dynamics.png
"""
import os
import sys
import numpy as np

def main():
    # 查找数据文件路径
    search_paths = [
        "test_rovibrational_dynamics.dat",
        "../test_rovibrational_dynamics.dat",
        "GeneralModule/test_rovibrational_dynamics.dat"
    ]
    data_file = None
    for p in search_paths:
        if os.path.isfile(p):
            data_file = p
            break

    if not data_file:
        print(f"[-] Data file not found in: {search_paths}")
        print("[*] Please run 'test_laser_rovibrational_control' first to generate data.")
        sys.exit(0)

    print(f"[+] Loading dynamics data from: {data_file}")
    data = np.loadtxt(data_file)
    time_fs = data[:, 0]
    pops = data[:, 1:]  # (nt, 12)

    # 绘制图形 (无显示设备时使用 Agg 后端)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
    })

    fig, axs = plt.subplots(2, 2, figsize=(13, 9), dpi=300)

    # 1. 主转振态布居时间跃迁曲线
    ax1 = axs[0, 0]
    ax1.plot(time_fs, pops[:, 0], label=r"$|v=0, J=0\rangle$ (Initial)", color="#1f77b4", lw=2.2)
    ax1.plot(time_fs, pops[:, 5], label=r"$|v=1, J=1\rangle$ (Target)", color="#d62728", lw=2.2)
    ax1.plot(time_fs, pops[:, 2], label=r"$|v=0, J=2\rangle$", color="#2ca02c", lw=1.6, ls="--")
    ax1.plot(time_fs, pops[:, 10], label=r"$|v=2, J=2\rangle$", color="#9467bd", lw=1.6, ls=":")
    ax1.set_xlabel("Time $t$ (fs)")
    ax1.set_ylabel("Population $P(v, J, t)$")
    ax1.set_title("(a) Coherent Rovibrational Population Transfer", fontweight="bold")
    ax1.set_xlim([time_fs[0], time_fs[-1]])
    ax1.set_ylim([-0.02, 1.05])
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9)

    # 2. 振动态总布居演化 P_v(t)
    ax2 = axs[0, 1]
    pop_v0 = np.sum(pops[:, 0:4], axis=1)
    pop_v1 = np.sum(pops[:, 4:8], axis=1)
    pop_v2 = np.sum(pops[:, 8:12], axis=1)
    tot_norm = pop_v0 + pop_v1 + pop_v2

    ax2.plot(time_fs, pop_v0, label=r"Total $v=0$", color="#1f77b4", lw=2.0)
    ax2.plot(time_fs, pop_v1, label=r"Total $v=1$", color="#d62728", lw=2.0)
    ax2.plot(time_fs, pop_v2, label=r"Total $v=2$", color="#ff7f0e", lw=2.0)
    ax2.plot(time_fs, tot_norm, label=r"Total Norm $\sum P$", color="black", lw=1.4, ls="-.")
    ax2.set_xlabel("Time $t$ (fs)")
    ax2.set_ylabel("Vibrational Manifold Population")
    ax2.set_title("(b) Vibrational Population Manifolds", fontweight="bold")
    ax2.set_xlim([time_fs[0], time_fs[-1]])
    ax2.set_ylim([-0.02, 1.08])
    ax2.grid(True, ls=":", alpha=0.6)
    ax2.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9)

    # 3. 二维 (v, J) 空间末态分布热力图
    ax3 = axs[1, 0]
    final_pop = pops[-1, :].reshape(3, 4)  # (v=0..2, J=0..3)
    im = ax3.imshow(final_pop, cmap="viridis", origin="lower", aspect="auto",
                    extent=[-0.5, 3.5, -0.5, 2.5], vmin=0, vmax=np.max(final_pop))
    cbar = fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
    cbar.set_label("Final Population $P(v, J)$")
    ax3.set_xticks([0, 1, 2, 3])
    ax3.set_yticks([0, 1, 2])
    ax3.set_xticklabels([f"$J={j}$" for j in range(4)])
    ax3.set_yticklabels([f"$v={v}$" for v in range(3)])
    ax3.set_title("(c) 2D Final State Distribution Matrix", fontweight="bold")

    for v in range(3):
        for j in range(4):
            val = final_pop[v, j]
            txt_color = "white" if val < 0.4 else "black"
            ax3.text(j, v, f"{val:.3f}", ha="center", va="center", color=txt_color, fontweight="bold")

    # 4. 各转振态最终布居条形柱状图
    ax4 = axs[1, 1]
    labels = [f"|{v},{j}⟩" for v in range(3) for j in range(4)]
    x_pos = np.arange(12)
    bars = ax4.bar(x_pos, pops[-1, :], color="#2b5c8f", edgecolor="black", alpha=0.85, width=0.65)
    # 高亮目标态
    bars[5].set_color("#d62728")
    bars[5].set_edgecolor("black")
    bars[5].set_alpha(1.0)

    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(labels, rotation=45, ha="right", fontsize=9.5)
    ax4.set_ylabel("Final Population")
    ax4.set_title(r"(d) Final Rovibrational State Yields (Target: $|1,1\rangle$)", fontweight="bold")
    ax4.set_ylim([0, 1.0])
    ax4.grid(True, axis="y", ls=":", alpha=0.6)

    plt.tight_layout()
    out_name = "result_rovibrational_dynamics.png"
    plt.savefig(out_name, dpi=300)
    print(f"[+] Publication-quality plot saved successfully to: {out_name}")

if __name__ == "__main__":
    main()
