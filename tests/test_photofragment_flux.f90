!> \brief 光解离碎片动能释放谱 (KER) 与自相关吸收谱单元测试
program test_photofragment_flux
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer, parameter :: nt = 500, nw = 100, ne = 80
    real(dp) :: dt, t_arr(nt), omega_arr(nw), e_kin_grid(ne)
    complex(dp) :: psi0(128), psit(128), c_hist(nt), psi_rd(nt), amp_e(ne)
    real(dp) :: abs_spec(nw), ker_spec(ne), ker_channels(ne, 2), yields(2), ratios(2)
    real(dp) :: omega_0, gamma_decay, t_val, sum_ratios, peak_w
    integer :: it, iw, ie, max_idx

    print '(A)', "=================================================="
    print '(A)', " GeneralModule Unit Tests: Photofragment & Flux   "
    print '(A)', "=================================================="

    ! 1. 测试波包自相关函数 C(0) = 1.0
    psi0 = (0.0_dp, 0.0_dp)
    psi0(64) = cmplx(1.0_dp / sqrt(0.1_dp), 0.0_dp, kind=dp)
    call assert_close("Autocorrelation C(0) = 1.0", 1.0_dp, abs(calculate_autocorrelation(psi0, psi0, 0.1_dp)), 1.0e-12_dp)

    ! 2. 模拟高斯阻尼自相关函数 C(t) = exp(-i * omega_0 * t - 0.5 * (gamma * t)^2)
    ! 对应吸收截面应在 omega = omega_0 处呈现高斯吸收峰
    omega_0 = 0.15_dp    ! 目标共振跃迁光子能量 (a.u.)
    gamma_decay = 0.02_dp
    dt = 0.5_dp * FS2AU

    do it = 1, nt
        t_val = real(it - 1, dp) * dt
        t_arr(it) = t_val
        c_hist(it) = exp(-EYE * omega_0 * t_val - 0.5_dp * (gamma_decay * t_val)**2)
    end do

    do iw = 1, nw
        omega_arr(iw) = 0.05_dp + real(iw - 1, dp) * (0.25_dp - 0.05_dp) / real(nw - 1, dp)
    end do

    call calculate_absorption_spectrum(t_arr, c_hist, 0.0_dp, omega_arr, abs_spec, dt)

    ! 查找吸收谱峰值对应频率
    max_idx = maxloc(abs_spec, dim=1)
    peak_w = omega_arr(max_idx)
    call assert_close("Absorption spectrum resonance peak at omega_0", omega_0, peak_w, 0.005_dp)

    ! 3. 测试渐近解离面能谱通量与 KER 谱
    ! 构造穿过渐近面 R=R_d 的波包时间包络: psi(R_d, t) ~ exp(-i * E_frag * t) * exp(-((t - t0)/tau)^2)
    do it = 1, nt
        t_val = t_arr(it)
        psi_rd(it) = exp(-EYE * 0.08_dp * t_val) * exp(-((t_val - 100.0_dp * FS2AU) / (25.0_dp * FS2AU))**2)
    end do

    do ie = 1, ne
        e_kin_grid(ie) = 0.01_dp + real(ie - 1, dp) * (0.15_dp - 0.01_dp) / real(ne - 1, dp)
    end do

    call energy_resolved_flux_amplitude(t_arr, psi_rd, e_kin_grid, amp_e, dt)
    call calculate_ker_spectrum(amp_e, e_kin_grid, 2000.0_dp, ker_spec)

    ! KER 峰值应位于 E_frag = 0.08 a.u.
    max_idx = maxloc(ker_spec, dim=1)
    call assert_close("KER kinetic energy peak at 0.08 a.u.", 0.08_dp, e_kin_grid(max_idx), 0.005_dp)

    ! 4. 测试双通道解离分支比统计: 假设通道 1 强度为通道 2 的 3 倍 (BR1 = 0.75, BR2 = 0.25)
    ker_channels(:, 1) = 3.0_dp * ker_spec(:)
    ker_channels(:, 2) = 1.0_dp * ker_spec(:)

    call calculate_branching_ratios(ker_channels, e_kin_grid, yields, ratios)
    sum_ratios = ratios(1) + ratios(2)
    call assert_close("Branching ratios sum = 1.0", 1.0_dp, sum_ratios, 1.0e-12_dp)
    call assert_close("Channel 1 branching ratio = 0.75", 0.75_dp, ratios(1), 1.0e-10_dp)
    call assert_close("Channel 2 branching ratio = 0.25", 0.25_dp, ratios(2), 1.0e-10_dp)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Photofragment & Flux Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some flux tests failed!"
        stop 1
    else
        print '(A)', "SUCCESS: All photofragment & flux tests passed."
    end if

contains

    subroutine assert_close(name, expected, actual, eps)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, eps
        real(dp) :: diff

        n_tests = n_tests + 1
        diff = abs(expected - actual)
        if (diff <= eps) then
            n_passed = n_passed + 1
            print '(A, A40, A)', " [PASS] ", name, ""
        else
            print '(A, A40, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_photofragment_flux
