!> \brief 阿秒瞬态吸收光谱 (ATAS) 与强场受驱自电离动力学模块
!> \details 基于相位微扰模型 (Phase Perturbation Model, PPM) 与光诱导态 (LIS) 理论，
!>          涵盖极紫外 (XUV) 阿秒脉冲激发自电离态、近红外 (NIR) 激光强场 Stark 频移、
!>          动态 Fano 线型 q-参数调控、光诱导态劈裂以及时延相关的双电子量子拍频。
!> \author LiHao
!> \date 2026-09-15
module mod_attosecond_transient_absorption
    use mod_constants, only: dp, PI, TWOPI, C_LIGHT, HBAR, EV2AU, AU2EV, FS2AU, AU2FS, &
                             W_CM2AU, NM2AU
    implicit none
    private

    public :: atas_state_t, atas_config_t
    public :: init_atas_helium_benchmark
    public :: calc_laser_dressed_fano_q, calc_light_induced_state_energy
    public :: calc_quantum_beat_period_fs, calc_atas_spectrum

    !> \brief 吸收共振能级参数类型
    type :: atas_state_t
        character(len=16) :: name         = "He_2s2p"   !< 能级名称 (如 He 2s2p 1P)
        real(dp)          :: energy_ev     = 60.15_dp    !< 共振本征能量 (eV)
        real(dp)          :: energy_au     = 0.0_dp      !< 共振本征能量 (a.u.)
        real(dp)          :: gamma_ev      = 0.037_dp    !< 自电离线宽 Gamma (eV), 对应 ~18 fs 寿命
        real(dp)          :: gamma_au      = 0.0_dp      !< 线宽 (a.u.)
        real(dp)          :: q_fano        = -2.80_dp    !< 无场稳态 Fano 不对称参数 q_0
        real(dp)          :: polarizability= 10.0_dp     !< 激发态极化率 alpha (a.u.)
        real(dp)          :: dipole_trans  = 0.25_dp     !< 跃迁偶极矩 (a.u.)
    end type atas_state_t

    !> \brief ATAS 激光场与时延实验配置
    type :: atas_config_t
        real(dp) :: nir_wavelength_nm   = 800.0_dp      !< NIR 调制激光波长 (nm)
        real(dp) :: nir_intensity_w_cm2 = 1.0e12_dp     !< NIR 激光光强 (W/cm^2)
        real(dp) :: nir_duration_fs     = 30.0_dp       !< NIR 脉冲持续时间 FWHM (fs)
        real(dp) :: xuv_duration_as     = 200.0_dp      !< XUV 阿秒单脉冲宽度 (as)
        real(dp) :: omega_nir_au        = 0.05696_dp    !< NIR 载波频率 (a.u.)
        real(dp) :: field_nir_au        = 0.0_dp        !< NIR 峰值场强 (a.u.)
    end type atas_config_t

contains

    ! ==========================================================================
    ! init_atas_helium_benchmark: 初始化氦原子 2s2p(^1P) 经典自电离态基准
    ! E_0 = 60.15 eV, Gamma = 0.037 eV (tau ~ 17.8 fs), q0 = -2.80
    ! ==========================================================================
    subroutine init_atas_helium_benchmark(bright_state, stat)
        type(atas_state_t), intent(out) :: bright_state
        integer, intent(out)            :: stat

        stat = 0
        bright_state%name          = "He_2s2p_1P"
        bright_state%energy_ev     = 60.150_dp
        bright_state%energy_au     = bright_state%energy_ev * EV2AU
        bright_state%gamma_ev      = 0.0370_dp
        bright_state%gamma_au      = bright_state%gamma_ev * EV2AU
        bright_state%q_fano        = -2.80_dp
        bright_state%polarizability= 12.50_dp
        bright_state%dipole_trans  = 0.225_dp
    end subroutine init_atas_helium_benchmark

    ! ==========================================================================
    ! calc_laser_dressed_fano_q:
    ! 根据强场相位微扰模型 (PPM) 计算受激红外光场调制后的有效 Fano 参数:
    !   q(tau) = [ q_0 + tan(Delta_Phi / 2) ] / [ 1 - q_0 * tan(Delta_Phi / 2) ]
    ! 展现激光驱动下的 Lorentzian 与 Fano 之间的连续相干变形 (Science 340, 716)
    ! ==========================================================================
    pure function calc_laser_dressed_fano_q(q0, delta_phi) result(q_eff)
        real(dp), intent(in) :: q0
        real(dp), intent(in) :: delta_phi
        real(dp)             :: q_eff

        real(dp) :: tan_half_phi, denom

        tan_half_phi = tan(0.5_dp * delta_phi)
        denom = 1.0_dp - q0 * tan_half_phi

        if (abs(denom) < 1.0e-7_dp) then
            q_eff = 1.0e7_dp * sign(1.0_dp, q0 + tan_half_phi)
        else
            q_eff = (q0 + tan_half_phi) / denom
        end if
    end function calc_laser_dressed_fano_q

    ! ==========================================================================
    ! calc_light_induced_state_energy:
    ! 计算光诱导态 (Light-Induced States, LIS) 能量:
    !   E_LIS = E_dark + hbar * omega_NIR 或避免交叉本征值
    ! ==========================================================================
    pure function calc_light_induced_state_energy(e_bright, e_dark, omega_nir, rabi_freq) &
        result(e_lis)
        real(dp), intent(in) :: e_bright
        real(dp), intent(in) :: e_dark
        real(dp), intent(in) :: omega_nir
        real(dp), intent(in) :: rabi_freq
        real(dp)             :: e_lis

        real(dp) :: delta_detuning, discr

        delta_detuning = (e_dark + omega_nir) - e_bright
        discr = sqrt(delta_detuning**2 + 4.0_dp * (rabi_freq**2))

        ! 返回避免交叉分裂下支或光修饰伴峰能量
        e_lis = 0.5_dp * (e_bright + (e_dark + omega_nir) - discr)
    end function calc_light_induced_state_energy

    ! ==========================================================================
    ! calc_quantum_beat_period_fs:
    ! 计算双态相干重叠引起的量子拍频周期:
    !   T_beat = 2 * pi * hbar / Delta_E
    ! ==========================================================================
    pure function calc_quantum_beat_period_fs(delta_e_ev) result(period_fs)
        real(dp), intent(in) :: delta_e_ev
        real(dp)             :: period_fs

        real(dp) :: delta_e_au

        if (delta_e_ev <= 0.0_dp) then
            period_fs = 0.0_dp
            return
        end if

        delta_e_au = delta_e_ev * EV2AU
        period_fs  = (TWOPI / delta_e_au) * AU2FS
    end function calc_quantum_beat_period_fs

    ! ==========================================================================
    ! calc_atas_spectrum:
    ! 模拟生成阿秒瞬态吸收谱矩阵 Delta OD(omega, tau)
    ! 包含 XUV-NIR 交叉时延从 tau_min 到 tau_max 的吸收线型动态演化
    ! ==========================================================================
    subroutine calc_atas_spectrum(bright_state, nir_intensity_w_cm2, nir_wavelength_nm, &
                                  n_energy, e_min_ev, e_max_ev, n_delay, &
                                  tau_min_fs, tau_max_fs, e_grid_ev, tau_grid_fs, spec_2d)
        type(atas_state_t), intent(in) :: bright_state
        real(dp), intent(in)           :: nir_intensity_w_cm2
        real(dp), intent(in)           :: nir_wavelength_nm
        integer, intent(in)            :: n_energy
        real(dp), intent(in)           :: e_min_ev, e_max_ev
        integer, intent(in)            :: n_delay
        real(dp), intent(in)           :: tau_min_fs, tau_max_fs
        real(dp), intent(out)          :: e_grid_ev(n_energy)
        real(dp), intent(out)          :: tau_grid_fs(n_delay)
        real(dp), intent(out)          :: spec_2d(n_energy, n_delay)

        integer  :: i_e, i_tau
        real(dp) :: de_step, dtau_step, tau, e_ev
        real(dp) :: i_au, f0_au, lambda_au, omega_au
        real(dp) :: ac_stark_au, delta_phi, q_eff, e_shifted_ev
        real(dp) :: epsilon, sigma_fano, sigma_unperturbed

        ! 构建能量与时延网格
        de_step = (e_max_ev - e_min_ev) / real(n_energy - 1, dp)
        do i_e = 1, n_energy
            e_grid_ev(i_e) = e_min_ev + real(i_e - 1, dp) * de_step
        end do

        dtau_step = (tau_max_fs - tau_min_fs) / real(n_delay - 1, dp)
        do i_tau = 1, n_delay
            tau_grid_fs(i_tau) = tau_min_fs + real(i_tau - 1, dp) * dtau_step
        end do

        ! 激光场参数
        i_au  = nir_intensity_w_cm2 * W_CM2AU
        f0_au = sqrt(i_au)
        lambda_au = nir_wavelength_nm * NM2AU
        omega_au  = TWOPI * 137.035999084_dp / lambda_au

        ! NIR 极化交流 Stark 位移: Delta E_Stark = -1/4 * alpha * F0^2
        ac_stark_au = 0.25_dp * bright_state%polarizability * (f0_au**2)

        do i_tau = 1, n_delay
            tau = tau_grid_fs(i_tau)

            ! 相位累积模型: 当 XUV 先于 NIR 到达 (tau > 0), 激发态经历 NIR 完整波包作用
            ! 当 tau < 0 (NIR 先于 XUV), 作用衰减
            if (tau >= 0.0_dp) then
                delta_phi = ac_stark_au * exp(- (tau / 25.0_dp)**2) * (25.0_dp * FS2AU)
            else
                delta_phi = ac_stark_au * exp(- (abs(tau) / 10.0_dp)**2) * (10.0_dp * FS2AU)
            end if

            q_eff = calc_laser_dressed_fano_q(bright_state%q_fano, delta_phi)
            e_shifted_ev = bright_state%energy_ev + ac_stark_au * AU2EV * exp(- (tau / 20.0_dp)**2)

            do i_e = 1, n_energy
                e_ev = e_grid_ev(i_e)

                ! 无量纲失谐能量 epsilon = 2 * (E - E0) / Gamma
                epsilon = 2.0_dp * (e_ev - e_shifted_ev) / max(bright_state%gamma_ev, 1.0e-6_dp)
                sigma_fano = ((q_eff + epsilon)**2) / (1.0_dp + epsilon**2)

                ! 无场参考线型
                epsilon = 2.0_dp * (e_ev - bright_state%energy_ev) / max(bright_state%gamma_ev, 1.0e-6_dp)
                sigma_unperturbed = ((bright_state%q_fano + epsilon)**2) / (1.0_dp + epsilon**2)

                ! 瞬态吸收差分吸光度: Delta OD proportional to (sigma - sigma_0)
                spec_2d(i_e, i_tau) = (sigma_fano - sigma_unperturbed) * 0.1_dp
            end do
        end do
    end subroutine calc_atas_spectrum

end module mod_attosecond_transient_absorption
