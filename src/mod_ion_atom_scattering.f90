!> \file mod_ion_atom_scattering.f90
!> \brief 冷离子-中性原子杂化碰撞与极化相互作用物理模块
!> \details 涵盖 1/r^4 长程极化势、Langevin 经典碰撞与反应捕获截面、
!>          极化势量子散射相移与微分散射、修正有效力程展开 (MERE)、
!>          以及 Paul 射频阱中微运动导致的碰撞致热动力学。
!> \author LiHao
module mod_ion_atom_scattering
    use mod_constants, only: dp, PI, TWOPI, AMU2AU, KB
    implicit none
    private

    public :: ion_atom_system_t
    public :: init_ion_atom_system
    public :: calc_ion_atom_potential
    public :: calc_langevin_critical_impact_parameter
    public :: calc_langevin_cross_section
    public :: calc_langevin_rate_coefficient
    public :: calc_ion_atom_phase_shift
    public :: calc_mere_phase_shift_s_wave
    public :: calc_rf_micromotion_heating

    !> 冷离子-中性原子碰撞体系特征参数
    type :: ion_atom_system_t
        real(dp) :: mass_ion_amu      !< 离子质量 (amu)
        real(dp) :: mass_atom_amu     !< 中性原子质量 (amu)
        real(dp) :: mu_amu            !< 约化质量 (amu)
        real(dp) :: mu_au             !< 约化质量 (a.u.)
        real(dp) :: charge_ion_e      !< 离子电荷 (通常为 +1.0)
        real(dp) :: alpha_atom_au     !< 中性原子静态偶极极化率 (a.u.)
        real(dp) :: c4_au             !< 长程极化相互作用系数 C4 (a.u.)
        real(dp) :: r_star_bohr       !< 极化特征长度标度 R* = sqrt(2*mu*C4/hbar^2) (Bohr)
        real(dp) :: e_star_kelvin     !< 极化特征能量标度 E* = hbar^2 / (2*mu*(R*)^2) (Kelvin)
        real(dp) :: k_langevin_cm3_s  !< 能量无关经典 Langevin 反应速率常数 (cm^3/s)
    end type ion_atom_system_t

contains

    !> \brief 初始化冷离子-中性原子杂化体系
    !> \param[out] sys 体系结构体
    !> \param[in] m_ion_amu 离子质量 (amu)
    !> \param[in] m_atom_amu 原子质量 (amu)
    !> \param[in] charge_e 离子电荷数 (通常为 +1.0)
    !> \param[in] alpha_au 中性原子静态极化率 (a.u.，如 Rb ~ 318.8 a.u., Li ~ 164.1 a.u., Yb ~ 141 a.u.)
    !> \param[out] stat 状态码 (0 成功, -1 参数非法)
    subroutine init_ion_atom_system(sys, m_ion_amu, m_atom_amu, charge_e, alpha_au, stat)
        type(ion_atom_system_t), intent(out) :: sys
        real(dp), intent(in)                 :: m_ion_amu, m_atom_amu
        real(dp), intent(in)                 :: charge_e, alpha_au
        integer, optional, intent(out)       :: stat

        real(dp) :: e_star_au

        if (present(stat)) stat = 0
        if (m_ion_amu <= 0.0_dp .or. m_atom_amu <= 0.0_dp .or. alpha_au <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        sys%mass_ion_amu = m_ion_amu
        sys%mass_atom_amu = m_atom_amu
        sys%charge_ion_e = charge_e
        sys%alpha_atom_au = alpha_au

        sys%mu_amu = (m_ion_amu * m_atom_amu) / (m_ion_amu + m_atom_amu)
        sys%mu_au = sys%mu_amu * AMU2AU

        ! 极化相互作用能量系数 V(r) = - C4 / (2 * r^4)
        ! 在原子单位下 C4 = alpha * q^2
        sys%c4_au = alpha_au * (charge_e**2)

        ! 特征散射长度尺度: R* = sqrt(2 * mu * C4) (a.u. 下 hbar = 1)
        sys%r_star_bohr = sqrt(2.0_dp * sys%mu_au * sys%c4_au)

        ! 特征散射能量尺度: E* = 1 / (2 * mu * (R*)^2)
        e_star_au = 1.0_dp / (2.0_dp * sys%mu_au * (sys%r_star_bohr**2))
        ! 转换到开尔文温度标度 (1 a.u. energy = 3.157750248e5 K)
        sys%e_star_kelvin = e_star_au * 3.157750248e5_dp

        ! 经典 Langevin 碰撞速率: K_L = 2*pi*sqrt(C4 / mu)
        ! a.u. 到 cm^3/s: 1 a.u. velocity = 2.18769126364e8 cm/s,
        !                 1 a.u. area = (0.529177210903e-8 cm)^2 = 2.8002852e-17 cm^2
        ! 1 a.u. rate = 2.18769126364e8 * 2.8002852e-17 = 6.126049e-9 cm^3/s
        sys%k_langevin_cm3_s = TWOPI * sqrt(sys%c4_au / sys%mu_au) * 6.126049e-9_dp
    end subroutine init_ion_atom_system

    !> \brief 计算离子-原子径向相互作用势 (12-4 模型势)
    !> \details V(r) = C12 / r^12 - C4 / (2 * r^4)
    pure function calc_ion_atom_potential(sys, r_bohr, r_core_bohr) result(v_pot)
        type(ion_atom_system_t), intent(in) :: sys
        real(dp), intent(in)                :: r_bohr
        real(dp), intent(in)                :: r_core_bohr
        real(dp)                            :: v_pot

        real(dp) :: c12_au, r_eff

        r_eff = max(0.5_dp, r_bohr)
        ! 取平衡位置 r_min ~ r_core_bohr: dV/dr = 0 => -12*C12/r^13 + 2*C4/r^5 = 0 => C12 = (C4/6) * r_core^8
        c12_au = (sys%c4_au / 6.0_dp) * (r_core_bohr**8)

        v_pot = c12_au / (r_eff**12) - 0.5_dp * sys%c4_au / (r_eff**4)
    end function calc_ion_atom_potential

    !> \brief 计算经典 Langevin 临界碰撞参数 b_c(E)
    !> \details 当碰撞参数 b <= b_c 时，径向有效势无离心势垒阻碍，粒子螺旋落入中心发生反应
    pure function calc_langevin_critical_impact_parameter(sys, energy_au) result(b_c)
        type(ion_atom_system_t), intent(in) :: sys
        real(dp), intent(in)                :: energy_au
        real(dp)                            :: b_c

        if (energy_au <= 1.0e-30_dp) then
            b_c = 1.0e10_dp
            return
        end if

        ! b_c = (2 * C4 / E)^(1/4)
        b_c = (2.0_dp * sys%c4_au / energy_au)**0.25_dp
    end function calc_langevin_critical_impact_parameter

    !> \brief 计算经典 Langevin 俘获截面 sigma_L(E) = pi * b_c^2
    pure function calc_langevin_cross_section(sys, energy_au) result(sigma_l)
        type(ion_atom_system_t), intent(in) :: sys
        real(dp), intent(in)                :: energy_au
        real(dp)                            :: sigma_l

        real(dp) :: b_c

        if (energy_au <= 1.0e-30_dp) then
            sigma_l = 0.0_dp
            return
        end if

        b_c = calc_langevin_critical_impact_parameter(sys, energy_au)
        sigma_l = PI * (b_c**2)
    end function calc_langevin_cross_section

    !> \brief 计算能量无关的 Langevin 反应速率常数 (cm^3/s)
    pure function calc_langevin_rate_coefficient(sys) result(k_l)
        type(ion_atom_system_t), intent(in) :: sys
        real(dp)                            :: k_l

        k_l = sys%k_langevin_cm3_s
    end function calc_langevin_rate_coefficient

    !> \brief 利用修正有效力程展开 (MERE) 计算冷离子-原子 s 波散射相移
    !> \details 对于 1/r^4 长程极化势，有效力程展开包含特征奇数次幂修正:
    !>          k * cot(delta_0) = -1/a_s + (pi / (3 * a_s^2)) * R* * k +
    !>                             (4 / (3 * a_s)) * (R*)^2 * k^2 * ln(R* * k / 4) + ...
    subroutine calc_mere_phase_shift_s_wave(sys, energy_au, a_s_bohr, delta0, stat)
        type(ion_atom_system_t), intent(in) :: sys
        real(dp), intent(in)                :: energy_au
        real(dp), intent(in)                :: a_s_bohr
        real(dp), intent(out)               :: delta0
        integer, optional, intent(out)      :: stat

        real(dp) :: k_wave, cot_delta, term1, term2, term3, arg_ln

        if (present(stat)) stat = 0
        if (energy_au <= 0.0_dp .or. abs(a_s_bohr) < 1.0e-12_dp) then
            if (present(stat)) stat = -1
            delta0 = 0.0_dp
            return
        end if

        k_wave = sqrt(2.0_dp * sys%mu_au * energy_au)
        term1 = - 1.0_dp / a_s_bohr
        term2 = (PI / (3.0_dp * (a_s_bohr**2))) * sys%r_star_bohr * k_wave

        arg_ln = max(1.0e-15_dp, sys%r_star_bohr * k_wave / 4.0_dp)
        term3 = (4.0_dp / (3.0_dp * a_s_bohr)) * (sys%r_star_bohr**2) * (k_wave**2) * log(arg_ln)

        cot_delta = (term1 + term2 + term3) / k_wave
        delta0 = atan2(1.0_dp, cot_delta)
        if (delta0 < 0.0_dp) delta0 = delta0 + PI
    end subroutine calc_mere_phase_shift_s_wave

    !> \brief 对指定分波求解离子-原子量子相移 delta_l(E)
    !> \param[in] sys 体系参数
    !> \param[in] energy_au 碰撞动能 (a.u.)
    !> \param[in] l 角动量量子数
    !> \param[in] r_match 渐近匹配外边界 (Bohr, 建议 >= 2.5 * R*)
    !> \param[in] r_core 短程排斥核半径 (Bohr)
    !> \param[out] phase_shift 散射相移 (rad)
    !> \param[out] stat 状态码
    subroutine calc_ion_atom_phase_shift(sys, energy_au, l, r_match, r_core, phase_shift, stat)
        type(ion_atom_system_t), intent(in) :: sys
        real(dp), intent(in)                :: energy_au
        integer, intent(in)                 :: l
        real(dp), intent(in)                :: r_match, r_core
        real(dp), intent(out)               :: phase_shift
        integer, optional, intent(out)      :: stat

        integer, parameter :: N_STEPS = 2000
        integer  :: i
        real(dp) :: dr, r_start, r_curr, k_wave
        real(dp) :: u_prev, u_curr, u_next
        real(dp) :: f_prev, f_curr, f_next
        real(dp) :: v_eff, kr_match, log_der
        real(dp) :: jl, nl, djl, dnl

        if (present(stat)) stat = 0
        if (energy_au <= 0.0_dp .or. r_match <= r_core) then
            if (present(stat)) stat = -1
            phase_shift = 0.0_dp
            return
        end if

        k_wave = sqrt(2.0_dp * sys%mu_au * energy_au)
        r_start = max(0.2_dp, r_core * 0.7_dp)
        dr = (r_match - r_start) / real(N_STEPS, dp)

        ! 初始化 Numerov 推进步
        u_prev = 0.0_dp
        u_curr = 1.0e-8_dp

        r_curr = r_start
        v_eff = calc_ion_atom_potential(sys, r_curr, r_core) + real(l * (l + 1), dp) / (2.0_dp * sys%mu_au * r_curr**2)
        f_prev = 2.0_dp * sys%mu_au * (energy_au - v_eff)

        r_curr = r_start + dr
        v_eff = calc_ion_atom_potential(sys, r_curr, r_core) + real(l * (l + 1), dp) / (2.0_dp * sys%mu_au * r_curr**2)
        f_curr = 2.0_dp * sys%mu_au * (energy_au - v_eff)

        do i = 2, N_STEPS
            r_curr = r_start + real(i, dp) * dr
            v_eff = calc_ion_atom_potential(sys, r_curr, r_core) + &
                    real(l * (l + 1), dp) / (2.0_dp * sys%mu_au * r_curr**2)
            f_next = 2.0_dp * sys%mu_au * (energy_au - v_eff)

            u_next = (2.0_dp * (1.0_dp - (5.0_dp / 12.0_dp) * (dr**2) * f_curr) * u_curr - &
                      (1.0_dp + (dr**2 / 12.0_dp) * f_prev) * u_prev) / &
                     (1.0_dp + (dr**2 / 12.0_dp) * f_next)

            u_prev = u_curr
            u_curr = u_next
            f_prev = f_curr
            f_curr = f_next
        end do

        ! 在渐近区 r_match 计算波函数的对数导数 log_der = u' / u
        log_der = (u_curr - u_prev) / (dr * max(1.0e-30_dp, abs(u_curr)))
        kr_match = k_wave * r_match

        ! 0 分波与 1 分波的解析 Riccati-Bessel 函数
        if (l == 0) then
            jl  = sin(kr_match)
            nl  = -cos(kr_match)
            djl = k_wave * cos(kr_match)
            dnl = k_wave * sin(kr_match)
        else
            jl  = sin(kr_match) / kr_match - cos(kr_match)
            nl  = -cos(kr_match) / kr_match - sin(kr_match)
            djl = k_wave * (sin(kr_match) + (2.0_dp / kr_match) * (cos(kr_match) - sin(kr_match) / kr_match))
            dnl = k_wave * (-cos(kr_match) + (2.0_dp / kr_match) * (sin(kr_match) + cos(kr_match) / kr_match))
        end if

        ! K 矩阵元 tan(delta) = (k * j'_l - log_der * j_l) / (k * n'_l - log_der * n_l)
        phase_shift = atan2(djl - log_der * jl, dnl - log_der * nl)
        if (phase_shift < 0.0_dp) phase_shift = phase_shift + PI
    end subroutine calc_ion_atom_phase_shift

    !> \brief 计算 Paul 射频阱中微运动驱动的碰撞致热动力学参数
    !> \details 离子在微波/射频交变电场中具有无耗散快微运动 (Micromotion)。
    !>          与超冷原子碰撞会随机打乱微运动相位，将射频场能量泵入体系产生持续发热。
    !> \param[in] trap_q Paul 阱 Mathieu 稳定性参数 q (通常 q ~ 0.1 - 0.4)
    !> \param[in] m_ion 离子质量 (amu)
    !> \param[in] m_atom 中性原子质量 (amu)
    !> \param[in] temp_atom_k 中性原子气体温度 (K)
    !> \param[in] coll_rate_hz 离子与中性原子碰撞频率 (s^-1)
    !> \param[out] t_limit_k 微运动致热平衡极限温度 (K)
    !> \param[out] heating_rate_k_s 初始致热升温速率 (K/s)
    pure subroutine calc_rf_micromotion_heating(trap_q, m_ion, m_atom, temp_atom_k, coll_rate_hz, &
                                               t_limit_k, heating_rate_k_s)
        real(dp), intent(in)  :: trap_q
        real(dp), intent(in)  :: m_ion, m_atom
        real(dp), intent(in)  :: temp_atom_k
        real(dp), intent(in)  :: coll_rate_hz
        real(dp), intent(out) :: t_limit_k
        real(dp), intent(out) :: heating_rate_k_s

        real(dp) :: mass_ratio, heat_factor

        mass_ratio = m_atom / max(1.0e-6_dp, m_ion)

        ! Cetina et al. (Phys. Rev. Lett. 109, 253201):
        ! 能量增益因子取决于原子与离子质量比以及 Mathieu q 参数:
        ! epsilon = (q^2 / 4) * (m_atom / m_ion) / (1 - (q^2 / 8) * (m_atom / m_ion))
        heat_factor = (trap_q**2 / 4.0_dp) * mass_ratio / max(0.01_dp, 1.0_dp - (trap_q**2 / 8.0_dp) * mass_ratio)

        ! 微运动加热使离子能量达到远高于环境原子的平衡平台
        t_limit_k = temp_atom_k * (1.0_dp + heat_factor / max(1.0e-5_dp, 1.0_dp - heat_factor))
        heating_rate_k_s = heat_factor * temp_atom_k * coll_rate_hz
    end subroutine calc_rf_micromotion_heating

end module mod_ion_atom_scattering
