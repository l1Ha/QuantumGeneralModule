! ==============================================================================
! GeneralModule Example 10: 超冷光缔合谱学与自由-束缚态分子生成模拟
! 理论参考文献:
!   - K. M. Jones, E. Tiesinga, P. D. Lett, P. S. Julienne, Rev. Mod. Phys. 78, 483 (2006)
!   - J. L. Bohn & P. S. Julienne, Phys. Rev. A 60, 414 (1999)
!   - R. Napolitano et al., Phys. Rev. Lett. 73, 1352 (1994)
! 演示内容：
! 1. 超冷 87Rb 原子 (T = 50 uK) 碰撞态与激发电子态长程分子振动态自由-束缚态 FC 积分
! 2. 缔合激光失谐 Delta 扫描 (-60 ~ +60 MHz) 模拟高分辨率光缔合线型谱
! 3. 激光强度依赖 (I = 20, 100, 500 W/cm^2) 的受激线宽展宽与损耗速率系数 K_PA
! 4. 双光子受激 Raman 缔合产生 X^1\Sigma_g^+ (v=0) 极性分子的有效拉比频率
! ==============================================================================
program ex10_photoassociation_spectroscopy
    use mod_constants, only: dp, PI, TO_AU, FROM_AU
    use mod_photoassociation
    use mod_io_utils, only: print_banner
    implicit none

    integer, parameter :: N_PTS = 300, N_DELTA = 31
    real(dp) :: r_grid(N_PTS), u_scat(N_PTS), chi_v(N_PTS), dip_mat(N_PTS)
    real(dp) :: delta_scan(N_DELTA), k_pa_low(N_DELTA), k_pa_mid(N_DELTA), k_pa_high(N_DELTA)
    real(dp) :: fc_int, fc_dens, gamma_stim_mid, omega_eff
    real(dp) :: r, dr, norm_sq
    integer  :: i, u, st

    call print_banner("GeneralModule Example 10: Photoassociation Spectroscopy", 70)

    ! --------------------------------------------------------------------------
    ! 1. 构造 87Rb 自由态与长程激发态分子振动态模型波函数
    ! --------------------------------------------------------------------------
    print *, ">> [1/3] Preparing Free Continuum & Bound Molecular Wavefunctions..."
    dr = (30.0_dp - 3.0_dp) / real(N_PTS - 1, dp)
    do i = 1, N_PTS
        r = 3.0_dp + real(i - 1, dp) * dr
        r_grid(i) = r
        ! 激发分子振动态 (长程 Condon 点位于 R_C = 12.0 a_0)
        chi_v(i)  = exp(-0.5_dp * ((r - 12.0_dp) / 1.5_dp)**2)
        ! 低能 s-波连续态散射波函数 (E = 50 uK)
        u_scat(i) = sin(0.05_dp * (r - 4.5_dp))
        ! 渐近跃迁偶极矩 mu_{eg} ~ 1.2 a.u.
        dip_mat(i) = 1.2_dp
    end do

    ! 激发束缚态归一化
    norm_sq = 0.0_dp
    do i = 1, N_PTS - 1
        norm_sq = norm_sq + 0.5_dp * (r_grid(i + 1) - r_grid(i)) * (chi_v(i)**2 + chi_v(i + 1)**2)
    end do
    chi_v = chi_v / sqrt(norm_sq)

    ! 计算自由-束缚态 Franck-Condon 跃迁密度
    call calc_free_bound_fc_overlap(r_grid, u_scat, chi_v, dip_mat, fc_int, fc_dens, st)
    print '(A, ES14.6, A)', "   Franck-Condon Overlap Integral I_FB : ", fc_int, " a.u."
    print '(A, ES14.6, A)', "   Franck-Condon Density f_FB          : ", fc_dens, " a.u.^-1"

    ! --------------------------------------------------------------------------
    ! 2. 模拟三组激光光强下光缔合线型谱随失谐 Delta (-60 ~ +60 MHz) 扫描
    ! --------------------------------------------------------------------------
    print *
    print *, ">> [2/3] Scanning Laser Detuning Delta for Photoassociation Line Shapes..."
    print *, "   Collision Temperature: T = 50 uK"
    print *, "   Natural Line Width   : Gamma_nat = 6.0 MHz"
    print *, "   Laser Intensities    : 20 W/cm^2 (Low), 100 W/cm^2 (Mid), 500 W/cm^2 (High)"
    print *, "----------------------------------------------------------------------"
    print '(A15, A18, A18, A18)', "Delta (MHz)", "K_PA [20 W/cm^2]", "K_PA [100 W/cm^2]", "K_PA [500 W/cm^2]"
    print *, "----------------------------------------------------------------------"

    ! I = 20 W/cm^2
    call calc_pa_spectrum_scan(50.0e-6_dp, 86.909_dp, -0.06_dp, 0.06_dp, N_DELTA, &
                               0.2_dp, 6.0_dp, delta_scan, k_pa_low, st)
    ! I = 100 W/cm^2
    call calc_pa_spectrum_scan(50.0e-6_dp, 86.909_dp, -0.06_dp, 0.06_dp, N_DELTA, &
                               1.0_dp, 6.0_dp, delta_scan, k_pa_mid, st)
    ! I = 500 W/cm^2
    call calc_pa_spectrum_scan(50.0e-6_dp, 86.909_dp, -0.06_dp, 0.06_dp, N_DELTA, &
                               5.0_dp, 6.0_dp, delta_scan, k_pa_high, st)

    open(newunit=u, file="photoassociation_spectrum.dat", status="replace", action="write")
    write(u, '(A)') "# Ultracold Photoassociation Spectrum vs Laser Detuning"
    write(u, '(A)') "# Detuning(MHz)  K_PA_20W(cm^3/s)  K_PA_100W(cm^3/s)  K_PA_500W(cm^3/s)"

    do i = 1, N_DELTA
        write(u, '(4ES16.7)') delta_scan(i) * 1000.0_dp, k_pa_low(i), k_pa_mid(i), k_pa_high(i)
        if (mod(i, 5) == 1 .or. i == N_DELTA) then
            print '(F12.2, 3ES18.6)', delta_scan(i) * 1000.0_dp, k_pa_low(i), k_pa_mid(i), k_pa_high(i)
        end if
    end do
    close(u)
    print *, ">> Photoassociation spectrum saved to: photoassociation_spectrum.dat"

    ! --------------------------------------------------------------------------
    ! 3. 双光子受激 Raman 缔合产生超冷基态分子
    ! --------------------------------------------------------------------------
    print *
    print *, ">> [3/3] Evaluating Two-Photon Raman STIRAP Association to Ground State..."
    ! Pump 光场 Rabi 频率 2*pi * 5 MHz, Stokes 光场 2*pi * 20 MHz, 单光子失谐 2*pi * 500 MHz
    omega_eff = calc_two_photon_raman_association_coupling(5.0e-9_dp, 2.0e-8_dp, 5.0e-7_dp)
    print '(A, ES14.6, A)', "   Effective Two-Photon Rabi Frequency Omega_eff : ", omega_eff, " a.u."
    print '(A, F10.3, A)',  "   Effective Ground Molecule Rabi Period T_2pi   : ", &
                            (2.0_dp * PI / omega_eff) * 2.4188843265857e-17 * 1.0e6, " us"

    print *
    print *, "======================================================================"
    print *, "   SUCCESS: Example 10 (Photoassociation) Completed Successfully!     "
    print *, "======================================================================"

end program ex10_photoassociation_spectroscopy
