!> \brief 激光脉冲时域合成单元测试
program test_laser_pulse
    use general_module
    implicit none

    integer :: n_tests = 0
    integer :: n_passed = 0
    real(dp) :: tol = 1.0e-8_dp

    type(pulse_config_t) :: p_cfg
    real(dp) :: env, ef, stark
    real(dp) :: t_arr(101), e_arr(101), env_arr(101)
    integer :: i

    print '(A)', "=================================================="
    print '(A)', "      GeneralModule Unit Tests: Laser Pulse       "
    print '(A)', "=================================================="

    ! 1. 高斯包络峰值与 FWHM 半高点测试
    p_cfg%shape_type = PULSE_GAUSSIAN
    p_cfg%field_peak = 0.05_dp
    p_cfg%freq_central = 0.057_dp ! ~800 nm (Ti:Sapphire)
    p_cfg%duration = 50.0_dp * FS2AU
    p_cfg%t_center = 0.0_dp
    p_cfg%cep_phase = 0.0_dp

    ! 峰值处 (t=0) 包络应为 1.0
    env = pulse_envelope(0.0_dp, p_cfg)
    call assert_close("Gaussian envelope peak (t=0)", 1.0_dp, env, tol)

    ! 半高全宽处 (t = +/- tau/2) 包络应为 0.5
    env = pulse_envelope(0.5_dp * p_cfg%duration, p_cfg)
    call assert_close("Gaussian envelope at FWHM edge", 0.5_dp, env, tol)

    ! 2. 瞬时电场峰值
    ef = pulse_electric_field(0.0_dp, p_cfg)
    call assert_close("Peak electric field at t=0", p_cfg%field_peak, ef, tol)

    ! 3. Sin^2 包络边界测试
    p_cfg%shape_type = PULSE_SIN2
    env = pulse_envelope(0.0_dp, p_cfg)
    call assert_close("Sin^2 envelope peak (t=0)", 1.0_dp, env, tol)

    env = pulse_envelope(p_cfg%duration, p_cfg)
    call assert_close("Sin^2 envelope at t=tau", 0.0_dp, env, tol)

    env = pulse_envelope(1.5_dp * p_cfg%duration, p_cfg)
    call assert_close("Sin^2 envelope outside [-tau, tau]", 0.0_dp, env, tol)

    ! 4. AC Stark 位移测试
    ! theta = 0 -> shift = -0.25 * alpha_par * E^2
    p_cfg%shape_type = PULSE_GAUSSIAN
    stark = pulse_stark_shift(0.0_dp, p_cfg, 10.0_dp, 5.0_dp, 0.0_dp)
    call assert_close("Stark shift at theta=0", -0.25_dp * 10.0_dp * (p_cfg%field_peak**2), stark, tol)

    ! theta = pi/2 -> shift = -0.25 * alpha_perp * E^2
    stark = pulse_stark_shift(0.0_dp, p_cfg, 10.0_dp, 5.0_dp, HALFPI)
    call assert_close("Stark shift at theta=pi/2", -0.25_dp * 5.0_dp * (p_cfg%field_peak**2), stark, tol)

    ! 5. 序列生成测试
    do i = 1, 101
        t_arr(i) = -100.0_dp * FS2AU + real(i - 1, dp) * (200.0_dp * FS2AU / 100.0_dp)
    end do
    call pulse_generate_timeseries(t_arr, p_cfg, e_arr, env_arr)
    call assert_close("Timeseries center envelope", 1.0_dp, env_arr(51), tol)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Laser Pulse Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some tests did not pass!"
        stop 1
    else
        print '(A)', "SUCCESS: All laser pulse tests passed."
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
            print '(A, A35, A)', " [PASS] ", name, ""
        else
            print '(A, A35, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_laser_pulse
