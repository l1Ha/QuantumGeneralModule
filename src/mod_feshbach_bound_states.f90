!> \brief 磁与光 Feshbach 共振：分子束缚态、闭通道成分与光致非弹性损耗模块
!> \details 基于多通道量子数亏损理论 (MQDT) 与两通道有效场论模型，
!>          涵盖磁 Feshbach 共振 (MFR) 的场强依赖散射长度 a(B)、共振特征长度 R*、
!>          弱束缚晕二聚体 (Halo Dimer) 普适及非普适结合能 E_b(B)、闭通道几率 Z(B)
!>          以及光 Feshbach 共振 (OFR) 的复散射长度与双体损耗率 K_2。
!> \author LiHao
!> \date 2026-09-13
module mod_feshbach_bound_states
    use mod_constants, only: dp, PI, AMU2AU, GAUSS2AU
    implicit none
    private

    public :: mfr_param_t, ofr_param_t
    public :: init_mfr_preset, init_ofr_param
    public :: calc_mfr_scattering_length, calc_mfr_r_star, calc_mfr_s_res
    public :: calc_mfr_bound_energy_universal, calc_mfr_bound_energy_coupled
    public :: calc_mfr_closed_channel_fraction
    public :: calc_ofr_complex_scattering_length, calc_ofr_inelastic_loss_rate

    !> \brief 磁 Feshbach 共振 (MFR) 参数类型
    type :: mfr_param_t
        character(len=16) :: name          = "6Li"         !< 共振系统标识
        real(dp)          :: mass_amu      = 6.015122_dp   !< 原子质量 (amu)
        real(dp)          :: mu_red_au     = 0.0_dp        !< 双原子约化质量 (a.u.)
        real(dp)          :: b0_gauss      = 832.2_dp      !< 共振极点磁场 B0 (Gauss)
        real(dp)          :: delta_b_gauss = -262.3_dp     !< 共振磁场宽度 Delta B (Gauss)
        real(dp)          :: a_bg_au       = -1405.0_dp    !< 背景散射长度 a_bg (Bohr)
        real(dp)          :: delta_mu_au   = 0.0_dp        !< 开闭通道差分磁矩 delta_mu (a.u.)
        real(dp)          :: c6_au         = 1393.0_dp     !< 范德瓦尔斯 C6 色散系数 (a.u.)
        real(dp)          :: r_vdw_au      = 0.0_dp        !< 范德瓦尔斯特征长度 R_vdW (Bohr)
        real(dp)          :: r_star_au     = 0.0_dp        !< 共振特征相互作用长度 R* (Bohr)
        real(dp)          :: s_res         = 0.0_dp        !< 无量纲共振强度参数 s_res
    end type mfr_param_t

    !> \brief 光 Feshbach 共振 (OFR) 参数类型
    type :: ofr_param_t
        character(len=16) :: name          = "87Rb"        !< 原子系统标识
        real(dp)          :: mass_amu      = 86.90918_dp   !< 原子质量 (amu)
        real(dp)          :: mu_red_au     = 0.0_dp        !< 约化质量 (a.u.)
        real(dp)          :: a_bg_au       = 100.0_dp      !< 背景散射长度 (Bohr)
        real(dp)          :: gamma_mol_hz  = 1.0e7_dp      !< 分子自发辐射线宽 (Hz)
        real(dp)          :: opt_length_au = 50.0_dp       !< 光学特征长度 l_opt (Bohr), 正比于光强
    end type ofr_param_t

contains

    ! ==========================================================================
    ! init_mfr_preset: 初始化著名实验磁 Feshbach 共振预设
    ! 包含 6Li (832 G 极宽共振), 40K (202 G 中等共振), 87Rb (1007 G 窄共振)
    ! ==========================================================================
    subroutine init_mfr_preset(preset_name, mfr, stat)
        character(len=*), intent(in)   :: preset_name
        type(mfr_param_t), intent(out) :: mfr
        integer, intent(out)           :: stat

        real(dp) :: delta_mu_mub

        stat = 0
        mfr%name = trim(preset_name)

        if (trim(preset_name) == "6Li" .or. trim(preset_name) == "Li6") then
            mfr%mass_amu      = 6.015122_dp
            mfr%b0_gauss      = 832.18_dp
            mfr%delta_b_gauss = -262.3_dp
            mfr%a_bg_au       = -1405.0_dp
            delta_mu_mub      = 2.0_dp           ! ~ 2.0 mu_B
            mfr%c6_au         = 1393.39_dp
        else if (trim(preset_name) == "40K" .or. trim(preset_name) == "K40") then
            mfr%mass_amu      = 39.963998_dp
            mfr%b0_gauss      = 202.10_dp
            mfr%delta_b_gauss = 7.04_dp
            mfr%a_bg_au       = 174.0_dp
            delta_mu_mub      = 1.68_dp          ! ~ 1.68 mu_B
            mfr%c6_au         = 3927.0_dp
        else if (trim(preset_name) == "87Rb" .or. trim(preset_name) == "Rb87") then
            mfr%mass_amu      = 86.90918_dp
            mfr%b0_gauss      = 1007.4_dp
            mfr%delta_b_gauss = 0.17_dp
            mfr%a_bg_au       = 100.5_dp
            delta_mu_mub      = 2.79_dp          ! ~ 2.79 mu_B
            mfr%c6_au         = 4707.0_dp
        else
            stat = -1
            mfr%mass_amu      = 6.015122_dp
            mfr%b0_gauss      = 832.18_dp
            mfr%delta_b_gauss = -262.3_dp
            mfr%a_bg_au       = -1405.0_dp
            delta_mu_mub      = 2.0_dp
            mfr%c6_au         = 1393.39_dp
        end if

        ! 约化质量 mu_red = m / 2 (同核双原子体系)
        mfr%mu_red_au = 0.5_dp * mfr%mass_amu * AMU2AU

        ! 玻尔磁子在原子单位下为 0.5 a.u. (e*hbar / (2*m_e) = 0.5)
        mfr%delta_mu_au = delta_mu_mub * 0.5_dp

        ! 范德瓦尔斯特征长度 R_vdW = 1/2 * (2 * mu_red * C6 / hbar^2)^(1/4)
        mfr%r_vdw_au = 0.5_dp * ((2.0_dp * mfr%mu_red_au * mfr%c6_au)**0.25_dp)

        ! 计算 R* 与无量纲共振强度 s_res
        mfr%r_star_au = calc_mfr_r_star(mfr)
        mfr%s_res     = calc_mfr_s_res(mfr)
    end subroutine init_mfr_preset

    ! ==========================================================================
    ! init_ofr_param: 初始化光 Feshbach 共振参数
    ! ==========================================================================
    subroutine init_ofr_param(name, mass_amu, a_bg_bohr, gamma_hz, l_opt_bohr, ofr, stat)
        character(len=*), intent(in)   :: name
        real(dp), intent(in)           :: mass_amu
        real(dp), intent(in)           :: a_bg_bohr
        real(dp), intent(in)           :: gamma_hz
        real(dp), intent(in)           :: l_opt_bohr
        type(ofr_param_t), intent(out) :: ofr
        integer, intent(out)           :: stat

        stat = 0
        if (mass_amu <= 0.0_dp .or. gamma_hz <= 0.0_dp) then
            stat = 1
            return
        end if

        ofr%name          = trim(name)
        ofr%mass_amu      = mass_amu
        ofr%mu_red_au     = 0.5_dp * mass_amu * AMU2AU
        ofr%a_bg_au       = a_bg_bohr
        ofr%gamma_mol_hz  = gamma_hz
        ofr%opt_length_au = l_opt_bohr
    end subroutine init_ofr_param

    ! ==========================================================================
    ! calc_mfr_scattering_length:
    ! 计算磁场 B (Gauss) 下的 s 波有效散射长度:
    !   a(B) = a_bg * [ 1 - Delta_B / (B - B_0) ]
    ! ==========================================================================
    pure function calc_mfr_scattering_length(mfr, b_gauss) result(a_au)
        type(mfr_param_t), intent(in) :: mfr
        real(dp), intent(in)          :: b_gauss
        real(dp)                      :: a_au

        real(dp) :: b_diff

        b_diff = b_gauss - mfr%b0_gauss
        if (abs(b_diff) < 1.0e-12_dp) then
            ! 处于极点共振位置
            if (mfr%delta_b_gauss * mfr%a_bg_au > 0.0_dp) then
                a_au = -1.0e12_dp
            else
                a_au = 1.0e12_dp
            end if
        else
            a_au = mfr%a_bg_au * (1.0_dp - mfr%delta_b_gauss / b_diff)
        end if
    end function calc_mfr_scattering_length

    ! ==========================================================================
    ! calc_mfr_r_star:
    ! 计算 Petrov-Chin 共振特征相互作用长度 R*:
    !   R* = hbar^2 / (2 * mu_red * |a_bg * delta_mu * Delta_B|)
    ! 在原子单位下 hbar = 1
    ! ==========================================================================
    pure function calc_mfr_r_star(mfr) result(r_star_au)
        type(mfr_param_t), intent(in) :: mfr
        real(dp)                      :: r_star_au

        real(dp) :: delta_b_au, denom

        delta_b_au = abs(mfr%delta_b_gauss) * GAUSS2AU
        denom = 2.0_dp * mfr%mu_red_au * abs(mfr%a_bg_au) * mfr%delta_mu_au * delta_b_au

        if (denom > 1.0e-30_dp) then
            r_star_au = 1.0_dp / denom
        else
            r_star_au = 1.0e15_dp
        end if
    end function calc_mfr_r_star

    ! ==========================================================================
    ! calc_mfr_s_res:
    ! 计算无量纲共振强度参数 s_res (Chin et al. RMP 2010):
    !   s_res = (a_bg / \bar{a}) * (delta_mu * Delta_B / \bar{E})
    ! 其中 \bar{a} = 0.955978 * R_vdW, \bar{E} = hbar^2 / (2 * mu_red * \bar{a}^2)
    ! s_res >> 1 为宽共振 (开通道主导), s_res << 1 为窄共振 (闭通道主导)
    ! ==========================================================================
    pure function calc_mfr_s_res(mfr) result(s_res)
        type(mfr_param_t), intent(in) :: mfr
        real(dp)                      :: s_res

        real(dp) :: a_bar, e_bar, delta_b_au

        a_bar = 0.9559784_dp * mfr%r_vdw_au
        if (a_bar <= 0.0_dp) then
            s_res = 0.0_dp
            return
        end if

        e_bar = 1.0_dp / (2.0_dp * mfr%mu_red_au * (a_bar**2))
        delta_b_au = abs(mfr%delta_b_gauss) * GAUSS2AU

        s_res = (abs(mfr%a_bg_au) / a_bar) * ((mfr%delta_mu_au * delta_b_au) / e_bar)
    end function calc_mfr_s_res

    ! ==========================================================================
    ! calc_mfr_bound_energy_universal:
    ! 普适弱束缚晕二聚体结合能 (零程近似 a(B) >> R_vdW):
    !   E_b = hbar^2 / (2 * mu_red * a(B)^2)
    ! 仅当 a(B) > 0 时存在二聚体分子束缚态
    ! ==========================================================================
    pure function calc_mfr_bound_energy_universal(mfr, b_gauss) result(e_b_au)
        type(mfr_param_t), intent(in) :: mfr
        real(dp), intent(in)          :: b_gauss
        real(dp)                      :: e_b_au

        real(dp) :: a_au

        a_au = calc_mfr_scattering_length(mfr, b_gauss)
        if (a_au <= 0.0_dp) then
            e_b_au = 0.0_dp
        else
            e_b_au = 1.0_dp / (2.0_dp * mfr%mu_red_au * (a_au**2))
        end if
    end function calc_mfr_bound_energy_universal

    ! ==========================================================================
    ! calc_mfr_bound_energy_coupled:
    ! 考虑有限共振特征长度 R* 的耦合通道精确分子束缚能 (Chin et al. RMP 2010):
    !   E_b = [ hbar^2 / (2 * mu_red * (R*)^2) ] * [ sqrt(1 + 2*R* / a(B)) - 1 ]^2
    ! 当 R* -> 0 时严格趋向普适解 E_b -> hbar^2 / (2*mu_red*a^2);
    ! 当 a >> R* 偏离普适区时自动过渡到准经典线性态 E_b ~ delta_mu * (B_0 - B).
    ! ==========================================================================
    pure function calc_mfr_bound_energy_coupled(mfr, b_gauss) result(e_b_au)
        type(mfr_param_t), intent(in) :: mfr
        real(dp), intent(in)          :: b_gauss
        real(dp)                      :: e_b_au

        real(dp) :: a_au, r_star, term, prefactor

        a_au = calc_mfr_scattering_length(mfr, b_gauss)
        if (a_au <= 0.0_dp) then
            e_b_au = 0.0_dp
            return
        end if

        r_star = mfr%r_star_au
        if (r_star <= 1.0e-12_dp) then
            ! 宽共振普适极限
            e_b_au = 1.0_dp / (2.0_dp * mfr%mu_red_au * (a_au**2))
            return
        end if

        prefactor = 1.0_dp / (2.0_dp * mfr%mu_red_au * (r_star**2))
        term = sqrt(1.0_dp + 2.0_dp * r_star / a_au) - 1.0_dp
        e_b_au = prefactor * (term**2)
    end function calc_mfr_bound_energy_coupled

    ! ==========================================================================
    ! calc_mfr_closed_channel_fraction:
    ! 分子态中闭通道分量占比 Z(B) = <psi|Q|psi>:
    !   由 Hellmann-Feynman 定理: Z(B) = (dE_b / dB) / delta_mu
    !   解析解: Z(B) = 1 / sqrt(1 + 2 * R* / a(B))
    ! 在共振极点 a -> inf 时 Z -> 0 (纯开通道分子); 在远共振区 a << R* 时 Z -> 1 (纯闭通道分子)
    ! ==========================================================================
    pure function calc_mfr_closed_channel_fraction(mfr, b_gauss) result(z_frac)
        type(mfr_param_t), intent(in) :: mfr
        real(dp), intent(in)          :: b_gauss
        real(dp)                      :: z_frac

        real(dp) :: a_au, r_star

        a_au = calc_mfr_scattering_length(mfr, b_gauss)
        if (a_au <= 0.0_dp) then
            z_frac = 0.0_dp
            return
        end if

        r_star = mfr%r_star_au
        z_frac = 1.0_dp - 1.0_dp / sqrt(1.0_dp + 2.0_dp * r_star / a_au)
    end function calc_mfr_closed_channel_fraction

    ! ==========================================================================
    ! calc_ofr_complex_scattering_length:
    ! 计算光 Feshbach 共振 (OFR) 下的有效复散射长度:
    !   \tilde{a}(Delta_L) = a_bg + \delta a(Delta_L) - i * b(Delta_L) / 2
    ! 其中:
    !   \delta a = l_opt * (gamma_mol * Delta_L) / [ Delta_L^2 + (gamma_mol / 2)^2 ]
    !   b        = 2 * l_opt * (gamma_mol / 2)^2 / [ Delta_L^2 + (gamma_mol / 2)^2 ]
    ! ==========================================================================
    pure subroutine calc_ofr_complex_scattering_length(ofr, delta_hz, a_real_au, a_imag_au)
        type(ofr_param_t), intent(in) :: ofr
        real(dp), intent(in)          :: delta_hz    !< 激光失谐 Delta_L (Hz)
        real(dp), intent(out)         :: a_real_au   !< 实部 Re(a) (Bohr)
        real(dp), intent(out)         :: a_imag_au   !< 虚部 Im(a) (Bohr), 负值对应损耗

        real(dp) :: gamma, denom, delta_a, b_loss

        gamma = ofr%gamma_mol_hz
        denom = delta_hz**2 + 0.25_dp * (gamma**2)

        if (denom > 1.0e-30_dp) then
            delta_a = ofr%opt_length_au * (gamma * delta_hz) / denom
            b_loss  = 2.0_dp * ofr%opt_length_au * (0.25_dp * gamma**2) / denom
        else
            delta_a = 0.0_dp
            b_loss  = 0.0_dp
        end if

        a_real_au = ofr%a_bg_au + delta_a
        a_imag_au = -0.5_dp * b_loss
    end subroutine calc_ofr_complex_scattering_length

    ! ==========================================================================
    ! calc_ofr_inelastic_loss_rate:
    ! 计算光致双体非弹性损失速率常数 K_2 (cm^3 / s):
    !   K_2 = (4 * pi * hbar / mu_red) * (-Im(\tilde{a}))
    ! ==========================================================================
    pure function calc_ofr_inelastic_loss_rate(ofr, delta_hz) result(k2_cm3_s)
        type(ofr_param_t), intent(in) :: ofr
        real(dp), intent(in)          :: delta_hz
        real(dp)                      :: k2_cm3_s

        real(dp) :: a_re, a_im, k2_au

        call calc_ofr_complex_scattering_length(ofr, delta_hz, a_re, a_im)

        ! 原子单位下的损耗速率常数: K_2(a.u.) = (4 * pi / mu_red) * (-a_imag)
        k2_au = (4.0_dp * PI / ofr%mu_red_au) * (-a_im)

        ! 单位转换: a.u. 体积/时间 -> cm^3 / s
        ! 1 a.u. length = 5.29177210903e-9 cm
        ! 1 a.u. time = 2.4188843265857e-17 s
        k2_cm3_s = k2_au * ((5.29177210903e-9_dp)**3) / 2.4188843265857e-17_dp
    end function calc_ofr_inelastic_loss_rate

end module mod_feshbach_bound_states
