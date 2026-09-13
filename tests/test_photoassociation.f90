! ==============================================================================
! GeneralModule Unit Tests: Ultracold Photoassociation Spectroscopy
! Testing:
!   1. Free-to-Bound Franck-Condon overlap integral & density f_FB(E, v')
!   2. Stimulated transition width \hbar\Gamma_{stim}(E, I) laser intensity scaling
!   3. Lorentzian resonance line shape & cross section peak at E = Delta
!   4. Maxwell-Boltzmann thermal average rate coefficient K_{PA}(T)
!   5. Detuning spectrum scan and resonance peak identification
!   6. Two-photon Raman association effective coupling \Omega_{eff}
! ==============================================================================
program test_photoassociation
    use mod_constants, only: dp, PI
    use mod_photoassociation
    implicit none

    integer :: n_pass = 0, n_total = 0
    real(dp) :: tol = 1.0e-5_dp

    ! 临时测试变量
    integer, parameter :: N_PTS = 200
    real(dp) :: r_grid(N_PTS), u_scat(N_PTS), chi_bound(N_PTS), dip_mat(N_PTS)
    real(dp) :: r0, sigma, k_wave, fc_int, fc_dens
    real(dp) :: gamma_stim_1, gamma_stim_2
    real(dp) :: detuning_au, gamma_nat_au, k_pa_therm
    real(dp) :: delta_scan(21), k_scan(21), omega_eff
    type(pa_transition_result_t) :: pa_res
    integer  :: i, st, peak_idx
    real(dp) :: max_k

    print *, "=================================================="
    print *, "   GeneralModule Unit Tests: Photoassociation     "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! 构造测试用模型波函数 (径向网格 r in [1.0, 20.0] a_0)
    ! --------------------------------------------------------------------------
    do i = 1, N_PTS
        r_grid(i) = 1.0_dp + real(i - 1, dp) * (19.0_dp / real(N_PTS - 1, dp))
        ! 激发态束缚态高斯波包 (模拟长程激发分子振动态, 居中于 r=8.0 a_0)
        chi_bound(i) = exp(-0.5_dp * ((r_grid(i) - 8.0_dp) / 1.0_dp)**2)
        ! 低能连续散射波函数 (模拟 s-波正弦振荡)
        k_wave = 0.1_dp
        u_scat(i) = sin(k_wave * (r_grid(i) - 2.0_dp))
        ! 跃迁偶极矩 mu_{eg}(r) 渐近趋于 1.0 a.u.
        dip_mat(i) = 1.0_dp
    end do

    ! 束缚态波函数归一化
    call normalize_wf(r_grid, chi_bound)

    ! --------------------------------------------------------------------------
    ! 1. 自由-束缚态 Franck-Condon 跃迁重叠积分检验
    ! --------------------------------------------------------------------------
    call calc_free_bound_fc_overlap(r_grid, u_scat, chi_bound, dip_mat, fc_int, fc_dens, st)
    call assert_true(st == 0, "Free-bound FC overlap status = 0", n_pass, n_total)
    call assert_true(fc_dens > 0.0_dp, "Free-bound FC overlap density > 0", n_pass, n_total)
    call assert_real_equal(fc_dens, fc_int**2, "FC density is square of overlap integral", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 2. 受激线宽 \hbar\Gamma_{stim}(E, I) 与激光光强正比线性度检验
    ! --------------------------------------------------------------------------
    gamma_stim_1 = calc_pa_stimulated_width(100.0_dp, fc_dens)  ! 100 W/cm^2
    gamma_stim_2 = calc_pa_stimulated_width(200.0_dp, fc_dens)  ! 200 W/cm^2

    call assert_true(gamma_stim_1 > 0.0_dp, "Stimulated width > 0", n_pass, n_total)
    call assert_real_equal(gamma_stim_2 / gamma_stim_1, 2.0_dp, &
                           "Stimulated width linearly scales with intensity (x2)", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 3. 单能光缔合截面与共振峰检验
    ! --------------------------------------------------------------------------
    detuning_au  = 1.0e-8_dp     ! 设定共振能量
    gamma_nat_au = 1.0e-9_dp     ! 自然线宽 (~6 MHz)

    ! 共振处 E = Delta
    call calc_pa_cross_section(detuning_au, 86.909_dp, detuning_au, gamma_stim_1, &
                               gamma_nat_au, pa_res, st)
    call assert_true(st == 0, "Resonance PA cross section status = 0", n_pass, n_total)
    call assert_true(pa_res%cross_section_au > 1.0e-3_dp, "Resonant PA cross section > 1e-3 a_0^2", n_pass, n_total)
    call assert_true(pa_res%rate_coeff_cm3_s > 0.0_dp, "Resonant rate coefficient K_PA > 0", n_pass, n_total)

    ! 远失谐处 E = 10 * Delta, 截面应当显著下降 (> 10 倍衰减)
    call calc_pa_cross_section(10.0_dp * detuning_au, 86.909_dp, detuning_au, gamma_stim_1, &
                               gamma_nat_au, pa_res, st)
    call assert_true(pa_res%cross_section_au < 1.0e-4_dp, &
                     "Off-resonant cross section strongly attenuated (< 1e-4)", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 4. 麦克斯韦-玻尔兹曼热系综平均光缔合速率 K_{PA}(T)
    ! --------------------------------------------------------------------------
    call calc_pa_thermal_rate_coefficient(1.0e-4_dp, 86.909_dp, detuning_au, gamma_stim_1, &
                                          gamma_nat_au, k_pa_therm, st)
    call assert_true(st == 0, "Thermal PA rate coefficient status = 0", n_pass, n_total)
    call assert_true(k_pa_therm > 0.0_dp, "Thermal PA rate coefficient > 0", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 5. 激光失谐频率扫描谱 K_{PA}(\Delta) 与共振峰定位
    ! --------------------------------------------------------------------------
    call calc_pa_spectrum_scan(1.0e-4_dp, 86.909_dp, -0.05_dp, 0.05_dp, 21, &
                               1.0_dp, 6.0_dp, delta_scan, k_scan, st)
    call assert_true(st == 0, "PA spectrum scan status = 0", n_pass, n_total)

    ! 寻找谱线最大值位置
    max_k = -1.0_dp
    peak_idx = 1
    do i = 1, 21
        if (k_scan(i) > max_k) then
            max_k = k_scan(i)
            peak_idx = i
        end if
    end do
    ! 共振峰应该出现在中心附近 (|Delta| < 0.02 GHz)
    call assert_true(abs(delta_scan(peak_idx)) < 0.02_dp, &
                     "PA resonance peak centered near zero detuning", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 6. 双光子受激 Raman 缔合耦合强度
    ! --------------------------------------------------------------------------
    ! \Omega_P = 1.0e-5, \Omega_S = 1.0e-5, \Delta_{exc} = 1.0e-3
    ! \Omega_{eff} = (1e-5 * 1e-5) / (2 * 1e-3) = 5.0e-8 a.u.
    omega_eff = calc_two_photon_raman_association_coupling(1.0e-5_dp, 1.0e-5_dp, 1.0e-3_dp)
    call assert_real_equal(omega_eff, 5.0e-8_dp, "Two-photon Raman Rabi frequency", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 汇总报告
    ! --------------------------------------------------------------------------
    print *, "--------------------------------------------------"
    print '(A, I3, A, I3, A)', "Photoassociation Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All photoassociation tests passed."
    else
        stop 1
    end if

contains

    subroutine normalize_wf(r, wf)
        real(dp), dimension(:), intent(in) :: r
        real(dp), dimension(:), intent(inout) :: wf
        integer  :: j, n
        real(dp) :: dr, norm_sq

        n = size(r)
        norm_sq = 0.0_dp
        do j = 1, n - 1
            dr = r(j + 1) - r(j)
            norm_sq = norm_sq + 0.5_dp * dr * (wf(j)**2 + wf(j + 1)**2)
        end do
        if (norm_sq > 1.0e-30_dp) then
            wf = wf / sqrt(norm_sq)
        end if
    end subroutine normalize_wf

    subroutine assert_real_equal(actual, expected, desc, p_count, t_count)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: desc
        integer, intent(inout) :: p_count, t_count
        t_count = t_count + 1
        if (abs(actual - expected) <= tol) then
            print '(A, A)', " [PASS] ", desc
            p_count = p_count + 1
        else
            print '(A, A, A, ES14.6, A, ES14.6)', " [FAIL] ", desc, &
                  " (Actual: ", actual, ", Expected: ", expected, ")"
        end if
    end subroutine assert_real_equal

    subroutine assert_true(cond, desc, p_count, t_count)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: desc
        integer, intent(inout) :: p_count, t_count
        t_count = t_count + 1
        if (cond) then
            print '(A, A)', " [PASS] ", desc
            p_count = p_count + 1
        else
            print '(A, A)', " [FAIL] ", desc
        end if
    end subroutine assert_true

end program test_photoassociation
