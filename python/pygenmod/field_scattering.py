"""
pygenmod: Field-Dressed Ultracold Scattering & Multi-Basis Coupled Channels
Provides tools for:
- Clebsch-Gordan coefficients for half-integers
- Single-atom Breit-Rabi Zeeman-hyperfine structure
- Four representation bases:
  1. Uncoupled Basis: |ms1, mi1, ms2, mi2, L, ML>
  2. f-Coupled Basis: |(s1 i1) f1 mf1, (s2 i2) f2 mf2, L, ML>
  3. Total Spin Coupled Basis: |(s1 s2) S, (i1 i2) I, F MF, L, ML>
  4. Field-Dressed Channel Basis: |\alpha_1(B), \alpha_2(B), L, ML>
- Exact unitary transformations between all bases
- Spin-exchange interaction s1 . s2 and asymptotic field Hamiltonian
- Magnetic Feshbach resonance fitting and visualization
"""

import math
from dataclasses import dataclass
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt

from .constants import AMU2AU

# ------------------------------------------------------------------------------
# 物理常数与基组枚举
# ------------------------------------------------------------------------------
BASIS_UNCOUPLED = 1
BASIS_F_COUPLED = 2
BASIS_TOTAL_SPIN = 3
BASIS_FIELD_DRESSED = 4

GAUSS2AU = 4.254382e-10           # 1 Gauss -> a.u.
AU2GAUSS = 1.0 / GAUSS2AU
MU_B_AU = 0.5                    # 玻尔磁子 (a.u.)
MU_N_AU = 0.5 / 1836.15267343     # 核磁子 (a.u.)
GHZ2AU = 1.51982984600e-7         # 1 GHz -> a.u.
AU2GHZ = 1.0 / GHZ2AU


def clebsch_gordan_half(two_j1: int, two_m1: int,
                        two_j2: int, two_m2: int,
                        two_j: int, two_m: int) -> float:
    """
    计算半整数角动量 Clebsch-Gordan 系数 <j1 m1, j2 m2 | j m>
    所有输入均为 2 倍角动量整数 (如 j=1/2 传入 1, m=-1/2 传入 -1)
    """
    if two_m1 + two_m2 != two_m:
        return 0.0
    if two_j < abs(two_j1 - two_j2) or two_j > two_j1 + two_j2:
        return 0.0
    if abs(two_m1) > two_j1 or abs(two_m2) > two_j2 or abs(two_m) > two_j:
        return 0.0

    # 三角系数 Delta(j1, j2, j) = sqrt((j1+j2-j)! (j1-j2+j)! (-j1+j2+j)! / (j1+j2+j+1)!)
    def fact_half(two_n: int) -> float:
        return float(math.factorial(two_n // 2))

    a = (two_j1 + two_j2 - two_j)
    b = (two_j1 - two_j2 + two_j)
    c = (-two_j1 + two_j2 + two_j)
    d = (two_j1 + two_j2 + two_j + 2)

    delta_sq = (fact_half(a) * fact_half(b) * fact_half(c)) / fact_half(d)
    delta = math.sqrt(delta_sq)

    prefactor = math.sqrt(float(two_j + 1)) * delta * math.sqrt(
        fact_half(two_j1 + two_m1) * fact_half(two_j1 - two_m1) *
        fact_half(two_j2 + two_m2) * fact_half(two_j2 - two_m2) *
        fact_half(two_j + two_m) * fact_half(two_j - two_m)
    )

    k_min = max(0, (two_j2 - two_j - two_m1) // 2, (two_j1 - two_j + two_m2) // 2)
    k_max = min((two_j1 + two_j2 - two_j) // 2, (two_j1 - two_m1) // 2, (two_j2 + two_m2) // 2)

    val_sum = 0.0
    for k in range(k_min, k_max + 1):
        num = (-1.0) ** k
        denom = (
            fact_half(2 * k) *
            fact_half(two_j1 + two_j2 - two_j - 2 * k) *
            fact_half(two_j1 - two_m1 - 2 * k) *
            fact_half(two_j2 + two_m2 - 2 * k) *
            fact_half(two_j - two_j2 + two_m1 + 2 * k) *
            fact_half(two_j - two_j1 - two_m2 + 2 * k)
        )
        val_sum += num / denom

    return prefactor * val_sum


@dataclass
class ColdAtom:
    name: str
    mass_amu: float
    mass_au: float
    two_s: int
    two_i: int
    g_s: float
    g_i: float
    a_hf_ghz: float
    a_hf_au: float
    dipole_debye: float = 0.0
    dipole_au: float = 0.0


def get_cold_atom_preset(name: str) -> ColdAtom:
    """获取常用冷原子/极性分子同位素预设参数."""
    two_s = 1  # 碱金属电子自旋 s=1/2
    g_s = 2.00231930436256
    dipole_debye = 0.0

    if name in ("6Li", "Li6"):
        mass_amu = 6.015122
        two_i = 2  # i = 1
        a_hf_ghz = 0.15213684
        g_i = -0.0004476540
    elif name in ("7Li", "Li7"):
        mass_amu = 7.016004
        two_i = 3  # i = 3/2
        a_hf_ghz = 0.40181204
        g_i = -0.001182213
    elif name in ("23Na", "Na23"):
        mass_amu = 22.989769
        two_i = 3
        a_hf_ghz = 0.88581306
        g_i = -0.000804610
    elif name in ("40K", "K40"):
        mass_amu = 39.963998
        two_i = 8  # i = 4
        a_hf_ghz = -0.2857308
        g_i = 0.00017649
    elif name in ("87Rb", "Rb87"):
        mass_amu = 86.9091805
        two_i = 3
        a_hf_ghz = 3.417341305452
        g_i = -0.0009951414
    elif name in ("133Cs", "Cs133"):
        mass_amu = 132.9054519
        two_i = 7  # i = 7/2
        a_hf_ghz = 2.2981579425
        g_i = -0.0003988539
    elif name in ("KRb", "40K87Rb"):
        mass_amu = 126.873
        two_s = 0
        two_i = 0
        g_s = 0.0
        g_i = 0.0
        a_hf_ghz = 0.0
        dipole_debye = 0.566
    else:
        raise ValueError(f"Unknown atom preset: {name}")

    mass_au = mass_amu * AMU2AU
    a_hf_au = a_hf_ghz * GHZ2AU
    dipole_au = dipole_debye / 2.541746473

    return ColdAtom(
        name=name, mass_amu=mass_amu, mass_au=mass_au,
        two_s=two_s, two_i=two_i, g_s=g_s, g_i=g_i,
        a_hf_ghz=a_hf_ghz, a_hf_au=a_hf_au,
        dipole_debye=dipole_debye, dipole_au=dipole_au
    )


def calc_breit_rabi_energies(atom: ColdAtom, b_gauss: float) -> tuple[np.ndarray, np.ndarray]:
    """
    单原子 Breit-Rabi 塞曼-超精细哈密顿量本征值与本征态求解
    H = a_hf * s . i + (g_s * mu_B * s_z - g_i * mu_N * i_z) * B
    """
    n_s = atom.two_s + 1
    n_i = atom.two_i + 1
    n_states = n_s * n_i

    two_ms = []
    two_mi = []
    for ms in range(atom.two_s, -atom.two_s - 1, -2):
        for mi in range(atom.two_i, -atom.two_i - 1, -2):
            two_ms.append(ms)
            two_mi.append(mi)

    two_ms = np.array(two_ms)
    two_mi = np.array(two_mi)

    h_mat = np.zeros((n_states, n_states), dtype=float)
    b_au = b_gauss * GAUSS2AU

    for idx in range(n_states):
        ms_val = two_ms[idx] / 2.0
        mi_val = two_mi[idx] / 2.0
        z_term = (atom.g_s * MU_B_AU * ms_val - atom.g_i * MU_N_AU * mi_val) * b_au
        h_mat[idx, idx] = z_term + atom.a_hf_au * (ms_val * mi_val)

        for jdx in range(n_states):
            if idx == jdx:
                continue
            # s+ i-
            if two_ms[idx] == two_ms[jdx] + 2 and two_mi[idx] == two_mi[jdx] - 2:
                sp_m = (0.5 * atom.a_hf_au *
                        math.sqrt((atom.two_s * (atom.two_s + 2) - two_ms[jdx] * (two_ms[jdx] + 2)) / 4.0) *
                        math.sqrt((atom.two_i * (atom.two_i + 2) - two_mi[jdx] * (two_mi[jdx] - 2)) / 4.0))
                h_mat[idx, jdx] = sp_m
            # s- i+
            elif two_ms[idx] == two_ms[jdx] - 2 and two_mi[idx] == two_mi[jdx] + 2:
                sm_p = (0.5 * atom.a_hf_au *
                        math.sqrt((atom.two_s * (atom.two_s + 2) - two_ms[jdx] * (two_ms[jdx] - 2)) / 4.0) *
                        math.sqrt((atom.two_i * (atom.two_i + 2) - two_mi[jdx] * (two_mi[jdx] + 2)) / 4.0))
                h_mat[idx, jdx] = sm_p

    evals, evecs = np.linalg.eigh(h_mat)
    return evals, evecs


@dataclass
class FieldChannel:
    idx: int
    two_ms1: int = 0
    two_mi1: int = 0
    two_ms2: int = 0
    two_mi2: int = 0
    two_f1: int = 0
    two_mf1: int = 0
    two_f2: int = 0
    two_mf2: int = 0
    two_S: int = 0
    two_I: int = 0
    two_F: int = 0
    two_MF: int = 0
    l_orb: int = 0
    m_l: int = 0
    two_Mtot: int = 0
    thresh_energy: float = 0.0


def build_field_collision_channels(atom1: ColdAtom, atom2: ColdAtom,
                                   basis_type: int, two_Mtot: int,
                                   l_max: int = 0) -> list[FieldChannel]:
    """给定守恒量 2*M_tot 与分波 l_max，构建碰撞通道列表."""
    channels = []
    ch_idx = 0

    for l_val in range(0, l_max + 1, 2):
        for ml_val in range(-l_val, l_val + 1):
            if basis_type == BASIS_UNCOUPLED:
                for ms1 in range(atom1.two_s, -atom1.two_s - 1, -2):
                    for mi1 in range(atom1.two_i, -atom1.two_i - 1, -2):
                        for ms2 in range(atom2.two_s, -atom2.two_s - 1, -2):
                            for mi2 in range(atom2.two_i, -atom2.two_i - 1, -2):
                                if ms1 + mi1 + ms2 + mi2 + 2 * ml_val == two_Mtot:
                                    ch_idx += 1
                                    channels.append(FieldChannel(
                                        idx=ch_idx, two_ms1=ms1, two_mi1=mi1,
                                        two_ms2=ms2, two_mi2=mi2,
                                        l_orb=l_val, m_l=ml_val, two_Mtot=two_Mtot
                                    ))

            elif basis_type == BASIS_F_COUPLED:
                for f1 in range(abs(atom1.two_s - atom1.two_i), atom1.two_s + atom1.two_i + 1, 2):
                    for mf1 in range(-f1, f1 + 1, 2):
                        for f2 in range(abs(atom2.two_s - atom2.two_i), atom2.two_s + atom2.two_i + 1, 2):
                            for mf2 in range(-f2, f2 + 1, 2):
                                if mf1 + mf2 + 2 * ml_val == two_Mtot:
                                    ch_idx += 1
                                    channels.append(FieldChannel(
                                        idx=ch_idx, two_f1=f1, two_mf1=mf1,
                                        two_f2=f2, two_mf2=mf2,
                                        l_orb=l_val, m_l=ml_val, two_Mtot=two_Mtot
                                    ))

            elif basis_type == BASIS_TOTAL_SPIN:
                for s_tot in range(abs(atom1.two_s - atom2.two_s), atom1.two_s + atom2.two_s + 1, 2):
                    for i_tot in range(abs(atom1.two_i - atom2.two_i), atom1.two_i + atom2.two_i + 1, 2):
                        for f_tot in range(abs(s_tot - i_tot), s_tot + i_tot + 1, 2):
                            for mf_tot in range(-f_tot, f_tot + 1, 2):
                                if mf_tot + 2 * ml_val == two_Mtot:
                                    ch_idx += 1
                                    channels.append(FieldChannel(
                                        idx=ch_idx, two_S=s_tot, two_I=i_tot,
                                        two_F=f_tot, two_MF=mf_tot,
                                        l_orb=l_val, m_l=ml_val, two_Mtot=two_Mtot
                                    ))

    return channels


def calc_basis_transform_matrix(atom1: ColdAtom, atom2: ColdAtom,
                                channels_src: list[FieldChannel],
                                channels_tgt: list[FieldChannel],
                                from_basis: int, to_basis: int,
                                b_gauss: float = 0.0) -> np.ndarray:
    """计算四大基组之间的严格幺正变换矩阵 U_{tgt, src}: |tgt> = sum_src U_{tgt, src} |src>."""
    n_ch = len(channels_src)
    if from_basis == to_basis:
        return np.eye(n_ch, dtype=float)

    # 1. Uncoupled -> f-coupled
    if from_basis == BASIS_UNCOUPLED and to_basis == BASIS_F_COUPLED:
        U = np.zeros((n_ch, n_ch), dtype=float)
        for i, tgt in enumerate(channels_tgt):
            for j, src in enumerate(channels_src):
                if tgt.l_orb != src.l_orb or tgt.m_l != src.m_l:
                    continue
                c1 = clebsch_gordan_half(atom1.two_s, src.two_ms1, atom1.two_i, src.two_mi1,
                                         tgt.two_f1, tgt.two_mf1)
                c2 = clebsch_gordan_half(atom2.two_s, src.two_ms2, atom2.two_i, src.two_mi2,
                                         tgt.two_f2, tgt.two_mf2)
                U[i, j] = c1 * c2
        return U

    elif from_basis == BASIS_F_COUPLED and to_basis == BASIS_UNCOUPLED:
        U = calc_basis_transform_matrix(atom1, atom2, channels_tgt, channels_src,
                                       BASIS_UNCOUPLED, BASIS_F_COUPLED, b_gauss)
        return U.T

    # 2. Uncoupled -> Total Spin
    elif from_basis == BASIS_UNCOUPLED and to_basis == BASIS_TOTAL_SPIN:
        U = np.zeros((n_ch, n_ch), dtype=float)
        for i, tgt in enumerate(channels_tgt):
            for j, src in enumerate(channels_src):
                if tgt.l_orb != src.l_orb or tgt.m_l != src.m_l:
                    continue
                ms = src.two_ms1 + src.two_ms2
                mi = src.two_mi1 + src.two_mi2
                cg_s = clebsch_gordan_half(atom1.two_s, src.two_ms1, atom2.two_s, src.two_ms2,
                                           tgt.two_S, ms)
                cg_i = clebsch_gordan_half(atom1.two_i, src.two_mi1, atom2.two_i, src.two_mi2,
                                           tgt.two_I, mi)
                cg_f = clebsch_gordan_half(tgt.two_S, ms, tgt.two_I, mi,
                                           tgt.two_F, tgt.two_MF)
                U[i, j] = cg_s * cg_i * cg_f
        return U

    elif from_basis == BASIS_TOTAL_SPIN and to_basis == BASIS_UNCOUPLED:
        U = calc_basis_transform_matrix(atom1, atom2, channels_tgt, channels_src,
                                       BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, b_gauss)
        return U.T

    # 3. Field-Dressed Channel Basis
    elif to_basis == BASIS_FIELD_DRESSED:
        h_asymp = build_asymptotic_hamiltonian(atom1, atom2, b_gauss, 0.0, from_basis, channels_src)
        evals, evecs = np.linalg.eigh(h_asymp)
        return evecs.T

    # 通用链式变换: from -> uncoupled -> to
    else:
        ch_unc = build_field_collision_channels(atom1, atom2, BASIS_UNCOUPLED,
                                               channels_src[0].two_Mtot, channels_src[0].l_orb)
        u1 = calc_basis_transform_matrix(atom1, atom2, channels_src, ch_unc, from_basis, BASIS_UNCOUPLED, b_gauss)
        u2 = calc_basis_transform_matrix(atom1, atom2, ch_unc, channels_tgt, BASIS_UNCOUPLED, to_basis, b_gauss)
        return u2 @ u1


def build_asymptotic_hamiltonian(atom1: ColdAtom, atom2: ColdAtom,
                                 b_gauss: float, e_field: float,
                                 basis_type: int, channels: list[FieldChannel]) -> np.ndarray:
    """构建渐近外场双原子哈密顿量 H_asymp(B, E) = h_1 + h_2."""
    ch_unc = build_field_collision_channels(atom1, atom2, BASIS_UNCOUPLED,
                                           channels[0].two_Mtot, channels[0].l_orb)
    n_ch = len(ch_unc)
    h_unc = np.zeros((n_ch, n_ch), dtype=float)
    b_au = b_gauss * GAUSS2AU

    for i, ch_i in enumerate(ch_unc):
        ms1 = ch_i.two_ms1 / 2.0
        mi1 = ch_i.two_mi1 / 2.0
        ms2 = ch_i.two_ms2 / 2.0
        mi2 = ch_i.two_mi2 / 2.0
        z1 = (atom1.g_s * MU_B_AU * ms1 - atom1.g_i * MU_N_AU * mi1) * b_au
        z2 = (atom2.g_s * MU_B_AU * ms2 - atom2.g_i * MU_N_AU * mi2) * b_au
        h_unc[i, i] = z1 + z2 + atom1.a_hf_au * (ms1 * mi1) + atom2.a_hf_au * (ms2 * mi2)

        for j, ch_j in enumerate(ch_unc):
            if i == j:
                continue
            # 原子 1 翻转
            if ch_i.two_ms2 == ch_j.two_ms2 and ch_i.two_mi2 == ch_j.two_mi2:
                if ch_i.two_ms1 == ch_j.two_ms1 + 2 and ch_i.two_mi1 == ch_j.two_mi1 - 2:
                    h_unc[i, j] = 0.5 * atom1.a_hf_au * \
                        math.sqrt((atom1.two_s * (atom1.two_s + 2) - ch_j.two_ms1 * (ch_j.two_ms1 + 2)) / 4.0) * \
                        math.sqrt((atom1.two_i * (atom1.two_i + 2) - ch_j.two_mi1 * (ch_j.two_mi1 - 2)) / 4.0)
                elif ch_i.two_ms1 == ch_j.two_ms1 - 2 and ch_i.two_mi1 == ch_j.two_mi1 + 2:
                    h_unc[i, j] = 0.5 * atom1.a_hf_au * \
                        math.sqrt((atom1.two_s * (atom1.two_s + 2) - ch_j.two_ms1 * (ch_j.two_ms1 - 2)) / 4.0) * \
                        math.sqrt((atom1.two_i * (atom1.two_i + 2) - ch_j.two_mi1 * (ch_j.two_mi1 + 2)) / 4.0)

            # 原子 2 翻转
            if ch_i.two_ms1 == ch_j.two_ms1 and ch_i.two_mi1 == ch_j.two_mi1:
                if ch_i.two_ms2 == ch_j.two_ms2 + 2 and ch_i.two_mi2 == ch_j.two_mi2 - 2:
                    h_unc[i, j] = 0.5 * atom2.a_hf_au * \
                        math.sqrt((atom2.two_s * (atom2.two_s + 2) - ch_j.two_ms2 * (ch_j.two_ms2 + 2)) / 4.0) * \
                        math.sqrt((atom2.two_i * (atom2.two_i + 2) - ch_j.two_mi2 * (ch_j.two_mi2 - 2)) / 4.0)
                elif ch_i.two_ms2 == ch_j.two_ms2 - 2 and ch_i.two_mi2 == ch_j.two_mi2 + 2:
                    h_unc[i, j] = 0.5 * atom2.a_hf_au * \
                        math.sqrt((atom2.two_s * (atom2.two_s + 2) - ch_j.two_ms2 * (ch_j.two_ms2 - 2)) / 4.0) * \
                        math.sqrt((atom2.two_i * (atom2.two_i + 2) - ch_j.two_mi2 * (ch_j.two_mi2 + 2)) / 4.0)

    if basis_type == BASIS_UNCOUPLED:
        return h_unc
    else:
        U = calc_basis_transform_matrix(atom1, atom2, ch_unc, channels, BASIS_UNCOUPLED, basis_type, b_gauss)
        return U @ h_unc @ U.T


def build_spin_exchange_matrix(channels: list[FieldChannel], basis_type: int) -> np.ndarray:
    """计算电子自旋交换算符 s1 . s2 矩阵元."""
    n_ch = len(channels)
    p_exc = np.zeros((n_ch, n_ch), dtype=float)

    if basis_type == BASIS_TOTAL_SPIN:
        for i, ch in enumerate(channels):
            s_val = ch.two_S / 2.0
            p_exc[i, i] = 0.5 * (s_val * (s_val + 1.0) - 0.75 - 0.75)
    else:
        for i, ch_i in enumerate(channels):
            ms1 = ch_i.two_ms1 / 2.0
            ms2 = ch_i.two_ms2 / 2.0
            p_exc[i, i] = ms1 * ms2
            for j, ch_j in enumerate(channels):
                if i == j:
                    continue
                if (ch_i.two_mi1 == ch_j.two_mi1 and ch_i.two_mi2 == ch_j.two_mi2 and
                    ch_i.l_orb == ch_j.l_orb and ch_i.m_l == ch_j.m_l):
                    if ch_i.two_ms1 == ch_j.two_ms1 + 2 and ch_i.two_ms2 == ch_j.two_ms2 - 2:
                        p_exc[i, j] = 0.5
                    elif ch_i.two_ms1 == ch_j.two_ms1 - 2 and ch_i.two_ms2 == ch_j.two_ms2 + 2:
                        p_exc[i, j] = 0.5

    return p_exc


def fit_feshbach_resonance_parameters(b_grid: np.ndarray, a_s_grid: np.ndarray) -> tuple[float, float, float]:
    """
    拟合磁 Feshbach 共振散射长度色散公式: a_s(B) = a_bg * [1 - Delta_B / (B - B_0)]
    返回 (B_0, Delta_B, a_bg)
    """
    n_b = len(b_grid)
    if n_b < 5:
        return 0.0, 0.0, a_s_grid[0]

    # 寻找正负跳跃最大处作为共振极点 B_0
    jumps = np.abs(np.diff(a_s_grid))
    max_idx = int(np.argmax(jumps))
    b_res = float(0.5 * (b_grid[max_idx] + b_grid[max_idx + 1]))

    # 背景散射长度: 两端平均
    a_bg = float(0.5 * (a_s_grid[0] + a_s_grid[-1]))
    if abs(a_bg) < 1e-10:
        a_bg = 1e-5

    # 寻找 a_s(B) = 0 处对应零点 B_zero = B_0 + Delta_B
    delta_b = float((b_grid[1] - b_grid[0]) * 2.0)
    for i in range(n_b - 1):
        if a_s_grid[i] * a_s_grid[i + 1] <= 0.0 and i != max_idx:
            b_zero = 0.5 * (b_grid[i] + b_grid[i + 1])
            delta_b = float(b_zero - b_res)
            break

    return b_res, delta_b, a_bg


def plot_breit_rabi_diagram(atom: ColdAtom, b_max_gauss: float = 200.0,
                            num_b: int = 150, save_path: Optional[str] = None):
    """绘制单原子 Zeeman-超精细 Breit-Rabi 能级图."""
    b_vals = np.linspace(0.0, b_max_gauss, num_b)
    n_states = (atom.two_s + 1) * (atom.two_i + 1)
    energies = np.zeros((num_b, n_states))

    for idx, b in enumerate(b_vals):
        evals, _ = calc_breit_rabi_energies(atom, b)
        energies[idx, :] = evals * AU2GHZ

    plt.figure(figsize=(8, 5))
    for s in range(n_states):
        plt.plot(b_vals, energies[:, s], lw=1.5)

    plt.xlabel("Magnetic Field $B$ (Gauss)", fontsize=12)
    plt.ylabel("Energy $E / h$ (GHz)", fontsize=12)
    plt.title(f"Breit-Rabi Diagram for $^{{{atom.name}}}$ ($s={atom.two_s}/2, i={atom.two_i}/2$)", fontsize=13)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()


def plot_magnetic_feshbach_resonance(b_grid: np.ndarray, a_s_grid: np.ndarray,
                                   b_res: Optional[float] = None,
                                   delta_b: Optional[float] = None,
                                   a_bg: Optional[float] = None,
                                   save_path: Optional[str] = None):
    """绘制磁 Feshbach 共振色散曲线 a_s(B) 及拟合对比."""
    if b_res is None or delta_b is None or a_bg is None:
        b_res, delta_b, a_bg = fit_feshbach_resonance_parameters(b_grid, a_s_grid)

    plt.figure(figsize=(8, 5))
    plt.scatter(b_grid, a_s_grid, color="black", s=25, label="Scattering Length $a_s(B)$")

    b_dense = np.linspace(b_grid[0], b_grid[-1], 500)
    a_fit = a_bg * (1.0 - delta_b / (b_dense - b_res))
    mask = np.abs(a_fit) < 1000.0
    plt.plot(b_dense[mask], a_fit[mask], "r--", lw=2,
             label=f"Fit: $B_0={b_res:.1f}$ G, $\\Delta B={delta_b:.1f}$ G, $a_{{bg}}={a_bg:.1f} a_0$")

    plt.axvline(b_res, color="blue", linestyle=":", alpha=0.7, label="Resonance Pole $B_0$")
    plt.axhline(0.0, color="gray", linestyle="-", alpha=0.5)
    plt.xlabel("Magnetic Field $B$ (Gauss)", fontsize=12)
    plt.ylabel("s-wave Scattering Length $a_s$ ($a_0$)", fontsize=12)
    plt.title("Magnetic Feshbach Resonance Profile", fontsize=13)
    plt.ylim(-500.0, 500.0)
    plt.legend(loc="upper right", frameon=True)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()
