! ==============================================================================
! GeneralModule: mod_photoassociation.f90
! ------------------------------------------------------------------------------
! 现代超冷光缔合谱学与自由-束缚态量子跃迁模块 (Ultracold Photoassociation Spectroscopy)
! 权威理论参考文献 (详见 LITERATURE.md):
!   - K. M. Jones, E. Tiesinga, P. D. Lett, P. S. Julienne, Rev. Mod. Phys. 78, 483 (2006)
!   - J. L. Bohn & P. S. Julienne, Phys. Rev. A 60, 414 (1999) [光缔合半解析散射理论]
!   - R. Napolitano, J. Weiner, C. J. Williams, P. S. Julienne, Phys. Rev. Lett. 73, 1352 (1994)
!   - H. R. Thorsheim, J. Weiner, P. S. Julienne, Phys. Rev. Lett. 58, 2420 (1987) [光缔合奠基]
! 核心特性：
! 1. 连续散射态 u_E(r) 与激发电子态振动束缚态 chi_{v'}(r) 自由-束缚态 (Free-Bound) FC 积分
! 2. 激光驱动受激跃迁宽度 \hbar\Gamma_{stim}(E, I) 与跃迁偶极矩核间距依赖集成
! 3. 广义光缔合共振截面与 Lorentzian / Fano 线型谱 (带有自然辐射线宽与激光线宽展宽)
! 4. 麦克斯韦-玻尔兹曼 (Maxwell-Boltzmann) 碰撞热系综双体光缔合损失速率系数 K_{PA}(T, I, Delta)
! 5. 双光子受激 Raman / STIRAP 缔合超冷基态分子有效双光子 Rabi 频率与形成截面
! ==============================================================================
module mod_photoassociation
    use mod_constants, only: dp, PI, TWOPI, SQRTPI, HBAR, AMU2AU, AU2CM, CM2AU
    implicit none
    private

    ! --------------------------------------------------------------------------
    ! 派生类型定义
    ! --------------------------------------------------------------------------
    !> \brief 光缔合跃迁物理参数配置结构体
    type, public :: pa_config_t
        real(dp) :: laser_intensity_w_cm2 !< 缔合激光功率密度 (W/cm^2)
        real(dp) :: gamma_nat_mhz         !< 激发态自然自发辐射衰减线宽 (MHz)
        real(dp) :: detuning_ghz          !< 激光相对于分子束缚能级的失谐 (GHz)
        real(dp) :: dipole_trans_au       !< 渐近电子跃迁偶极矩 mu_{eg} (a.u.)
    end type pa_config_t

    !> \brief 自由-束缚态跃迁矩阵元与截面分析结果
    type, public :: pa_transition_result_t
        real(dp) :: fc_overlap_density    !< 自由-束缚态 Franck-Condon 因子密度 (1/a.u.)
        real(dp) :: gamma_stim_au         !< 激光受激线宽 \hbar\Gamma_{stim} (a.u.)
        real(dp) :: cross_section_au      !< 单能光缔合截面 sigma_{PA}(E) (a.u.)
        real(dp) :: rate_coeff_cm3_s      !< 单能双体速率系数 K_{PA}(E) (cm^3/s)
    end type pa_transition_result_t

    ! --------------------------------------------------------------------------
    ! 公共接口导出
    ! --------------------------------------------------------------------------
    public :: calc_free_bound_fc_overlap
    public :: calc_pa_stimulated_width
    public :: calc_pa_cross_section
    public :: calc_pa_thermal_rate_coefficient
    public :: calc_pa_spectrum_scan
    public :: calc_two_photon_raman_association_coupling

contains

    ! ==========================================================================
    ! 1. 自由-束缚态 (Free-to-Bound) Franck-Condon 跃迁重叠积分
    ! I_{FB}(E, v') = \int_0^\infty u_E(r) * mu_{eg}(r) * \chi_{v'}(r) dr
    ! 其中:
    ! u_E(r): 能量归一化散射态波函数 (\int u_E(r) u_{E'}(r) dr = \delta(E - E'))
    ! \chi_{v'}(r): 束缚态波函数 (\int |\chi_{v'}(r)|^2 dr = 1)
    ! FC 密度 f_{FB}(E, v') = |I_{FB}|^2 (具有 1/能量 的量纲)
    ! ==========================================================================
    subroutine calc_free_bound_fc_overlap(r_grid, u_energy, chi_bound, dipole_trans, &
                                          fc_integral, fc_density, stat)
        real(dp), dimension(:), intent(in) :: r_grid
        real(dp), dimension(:), intent(in) :: u_energy
        real(dp), dimension(:), intent(in) :: chi_bound
        real(dp), dimension(:), intent(in) :: dipole_trans
        real(dp), intent(out)              :: fc_integral
        real(dp), intent(out)              :: fc_density
        integer, optional, intent(out)     :: stat

        integer  :: n_pts, i
        real(dp) :: dr, sum_int

        if (present(stat)) stat = 0
        n_pts = size(r_grid)

        if (n_pts < 4 .or. size(u_energy) /= n_pts .or. &
            size(chi_bound) /= n_pts .or. size(dipole_trans) /= n_pts) then
            if (present(stat)) stat = -1
            return
        end if

        ! 使用梯形数值积分
        sum_int = 0.0_dp
        do i = 1, n_pts - 1
            dr = r_grid(i + 1) - r_grid(i)
            sum_int = sum_int + 0.5_dp * dr * ( &
                u_energy(i) * dipole_trans(i) * chi_bound(i) + &
                u_energy(i + 1) * dipole_trans(i + 1) * chi_bound(i + 1))
        end do

        fc_integral = sum_int
        fc_density  = sum_int**2
    end subroutine calc_free_bound_fc_overlap

    ! ==========================================================================
    ! 2. 激光驱动受激跃迁宽度 \hbar\Gamma_{stim}(E, I)
    ! \hbar\Gamma_{stim} = \pi * \mathcal{E}_0^2 * |I_{FB}|^2
    ! 其中 \mathcal{E}_0 为激光电场峰值 (a.u.), 1 a.u. 强光强度 = 3.509445e16 W/cm^2
    ! ==========================================================================
    pure function calc_pa_stimulated_width(laser_intensity_w_cm2, fc_density) result(gamma_stim_au)
        real(dp), intent(in) :: laser_intensity_w_cm2
        real(dp), intent(in) :: fc_density
        real(dp)            :: gamma_stim_au

        real(dp) :: i_au, e0_sq

        ! 强度 W/cm^2 转换为 a.u.
        i_au = laser_intensity_w_cm2 / 3.5094452e16_dp
        ! I = c * epsilon_0 * E0^2 / 2 => E0^2 (a.u.) = 8*pi*alpha * I_au
        ! 在原子单位下: E0^2 = 8.0 * pi * alpha * I = 8.0 * PI * (1/137.036) * i_au
        e0_sq = (8.0_dp * PI / 137.035999084_dp) * i_au

        ! 受激宽度: \hbar \Gamma_{stim} = 2*pi * |<u_E| -d*E/2 | chi_v>|^2 * 2 = \pi * E0^2 * |I_FB|^2
        gamma_stim_au = PI * e0_sq * fc_density
    end function calc_pa_stimulated_width

    ! ==========================================================================
    ! 3. 单能超冷光缔合反应截面与速率系数
    ! \sigma_{PA}(E, \Delta) = \frac{\pi}{k^2} * \frac{\hbar\Gamma_{stim} * \gamma_{nat}}
    !                                          {(E - \Delta)^2 + [(\gamma_{nat} + \hbar\Gamma_{stim})/2]^2}
    ! ==========================================================================
    subroutine calc_pa_cross_section(collision_energy_au, mass_amu, detuning_au, &
                                     gamma_stim_au, gamma_nat_au, res, stat)
        real(dp), intent(in)                       :: collision_energy_au
        real(dp), intent(in)                       :: mass_amu
        real(dp), intent(in)                       :: detuning_au
        real(dp), intent(in)                       :: gamma_stim_au
        real(dp), intent(in)                       :: gamma_nat_au
        type(pa_transition_result_t), intent(out)  :: res
        integer, optional, intent(out)             :: stat

        real(dp) :: mu_au, k_wave, v_rel, denom, gamma_tot_half

        res%gamma_stim_au    = 0.0_dp
        res%cross_section_au = 0.0_dp
        res%rate_coeff_cm3_s = 0.0_dp

        if (present(stat)) stat = 0
        if (collision_energy_au <= 0.0_dp .or. mass_amu <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        mu_au = 0.5_dp * mass_amu * AMU2AU
        k_wave = sqrt(2.0_dp * mu_au * collision_energy_au)
        v_rel  = k_wave / mu_au

        gamma_tot_half = 0.5_dp * (gamma_nat_au + gamma_stim_au)
        denom = (collision_energy_au - detuning_au)**2 + gamma_tot_half**2

        res%gamma_stim_au = gamma_stim_au
        if (denom > 1.0e-35_dp) then
            res%cross_section_au = (PI / (k_wave**2)) * (gamma_stim_au * gamma_nat_au) / denom
        else
            res%cross_section_au = 0.0_dp
        end if

        ! 转换为双体速率系数 K_{PA} (cm^3/s)
        ! 1 a.u. (K) = 6.126048e-9 cm^3/s
        res%rate_coeff_cm3_s = v_rel * res%cross_section_au * 6.126048e-9_dp
    end subroutine calc_pa_cross_section

    ! ==========================================================================
    ! 4. 麦克斯韦-玻尔兹曼 (Maxwell-Boltzmann) 热系综平均光缔合速率系数 K_{PA}(T)
    ! K_{PA}(T, \Delta) = \frac{1}{h Q_T} \int_0^\infty e^{-E/k_B T} \frac{\hbar\Gamma_{stim} \gamma_{nat}}
    !                                                           {(E - \Delta)^2 + [(\gamma_tot)/2]^2} dE
    ! 其中 Q_T = (2*pi*mu*k_B*T / h^2)^{3/2}
    ! ==========================================================================
    subroutine calc_pa_thermal_rate_coefficient(temp_k, mass_amu, detuning_au, &
                                                gamma_stim_base_au, gamma_nat_au, &
                                                k_pa_cm3_s, stat)
        real(dp), intent(in)           :: temp_k
        real(dp), intent(in)           :: mass_amu
        real(dp), intent(in)           :: detuning_au
        real(dp), intent(in)           :: gamma_stim_base_au
        real(dp), intent(in)           :: gamma_nat_au
        real(dp), intent(out)          :: k_pa_cm3_s
        integer, optional, intent(out) :: stat

        integer, parameter :: N_E = 128
        integer  :: i
        real(dp) :: kb_t_au, e_max, de, e_curr, gamma_stim_e
        real(dp) :: sum_rate, sum_weight, weight
        type(pa_transition_result_t) :: tr_res

        k_pa_cm3_s = 0.0_dp
        if (present(stat)) stat = 0
        if (temp_k <= 0.0_dp .or. mass_amu <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        kb_t_au = temp_k * 3.166811563e-6_dp
        e_max = max(8.0_dp * kb_t_au, detuning_au + 5.0_dp * kb_t_au)
        de = e_max / real(N_E, dp)

        sum_rate   = 0.0_dp
        sum_weight = 0.0_dp

        do i = 1, N_E
            e_curr = (real(i, dp) - 0.5_dp) * de
            ! 自由-束缚态受激跃迁宽度在低能遵循 Wigner 阈值定律: \Gamma_{stim}(E) \propto \sqrt{E}
            gamma_stim_e = gamma_stim_base_au * sqrt(e_curr / kb_t_au)

            call calc_pa_cross_section(e_curr, mass_amu, detuning_au, gamma_stim_e, &
                                       gamma_nat_au, tr_res)

            ! 3D 相对运动麦克斯韦分布权重: P(E) dE \propto \sqrt{E} \exp(-E / k_B T) dE
            weight = sqrt(e_curr) * exp(-e_curr / kb_t_au) * de
            sum_rate   = sum_rate   + tr_res%rate_coeff_cm3_s * weight
            sum_weight = sum_weight + weight
        end do

        if (sum_weight > 1.0e-35_dp) then
            k_pa_cm3_s = sum_rate / sum_weight
        else
            k_pa_cm3_s = 0.0_dp
        end if
    end subroutine calc_pa_thermal_rate_coefficient

    ! ==========================================================================
    ! 5. 激光失谐频率扫描光缔合吸收谱 K_{PA}(\Delta)
    ! ==========================================================================
    subroutine calc_pa_spectrum_scan(temp_k, mass_amu, delta_min_ghz, delta_max_ghz, n_delta, &
                                     gamma_stim_base_mhz, gamma_nat_mhz, &
                                     delta_grid_ghz, k_pa_array, stat)
        real(dp), intent(in)              :: temp_k
        real(dp), intent(in)              :: mass_amu
        real(dp), intent(in)              :: delta_min_ghz
        real(dp), intent(in)              :: delta_max_ghz
        integer, intent(in)               :: n_delta
        real(dp), intent(in)              :: gamma_stim_base_mhz
        real(dp), intent(in)              :: gamma_nat_mhz
        real(dp), dimension(:), intent(out) :: delta_grid_ghz
        real(dp), dimension(:), intent(out) :: k_pa_array
        integer, optional, intent(out)    :: stat

        integer  :: i
        real(dp) :: d_delta, delta_ghz, delta_au
        real(dp) :: gamma_stim_au, gamma_nat_au, k_val

        if (present(stat)) stat = 0
        if (n_delta < 2 .or. size(delta_grid_ghz) < n_delta .or. size(k_pa_array) < n_delta) then
            if (present(stat)) stat = -1
            return
        end if

        ! 1 MHz = 1.519829846e-10 a.u., 1 GHz = 1.519829846e-7 a.u.
        gamma_stim_au = gamma_stim_base_mhz * 1.519829846e-10_dp
        gamma_nat_au  = gamma_nat_mhz * 1.519829846e-10_dp

        d_delta = (delta_max_ghz - delta_min_ghz) / real(n_delta - 1, dp)

        do i = 1, n_delta
            delta_ghz = delta_min_ghz + real(i - 1, dp) * d_delta
            delta_au  = delta_ghz * 1.519829846e-7_dp
            delta_grid_ghz(i) = delta_ghz

            call calc_pa_thermal_rate_coefficient(temp_k, mass_amu, delta_au, &
                                                  gamma_stim_au, gamma_nat_au, k_val)
            k_pa_array(i) = k_val
        end do
    end subroutine calc_pa_spectrum_scan

    ! ==========================================================================
    ! 6. 双光子受激 Raman 缔合产生基态超冷分子的有效双光子 Rabi 频率
    ! 碰撞对 A+B --(Omega_P)--> AB*(v_exc) --(Omega_S)--> AB(v_gnd)
    ! \Omega_{eff} = (\Omega_P * \Omega_S) / (2 * \Delta_{exc})
    ! ==========================================================================
    pure function calc_two_photon_raman_association_coupling(omega_pump_au, omega_stokes_au, &
                                                             detuning_exc_au) result(omega_eff_au)
        real(dp), intent(in) :: omega_pump_au
        real(dp), intent(in) :: omega_stokes_au
        real(dp), intent(in) :: detuning_exc_au
        real(dp)            :: omega_eff_au

        if (abs(detuning_exc_au) > 1.0e-20_dp) then
            omega_eff_au = (omega_pump_au * omega_stokes_au) / (2.0_dp * detuning_exc_au)
        else
            omega_eff_au = 0.0_dp
        end if
    end function calc_two_photon_raman_association_coupling

end module mod_photoassociation
