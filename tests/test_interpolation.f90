!> \brief 样条插值与势能面外推单元测试
program test_interpolation
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    type(spline_1d_t) :: spl
    integer, parameter :: n_pts = 21
    real(dp) :: x(n_pts), y(n_pts)
    real(dp) :: x_test, y_test, dy_test, y_exact, dy_exact, err
    integer :: i
    real(dp) :: r_raw(5), v_raw(5), r_grid(20), v_grid(20)

    print '(A)', "=================================================="
    print '(A)', "     GeneralModule Unit Tests: Interpolation      "
    print '(A)', "=================================================="

    ! 1. 测试标准三次样条插值: y(x) = sin(x), x ∈ [0, pi]
    do i = 1, n_pts
        x(i) = real(i - 1, dp) * (PI / real(n_pts - 1, dp))
        y(i) = sin(x(i))
    end do

    call spline_init(x, y, spl, BC_NATURAL)

    ! 节点精确还原校验
    err = 0.0_dp
    do i = 1, n_pts
        err = max(err, abs(spline_eval(spl, x(i)) - y(i)))
    end do
    call assert_close("Spline node value reconstruction", 0.0_dp, err, 1.0e-12_dp)

    ! 区间中间点插值精度校验: x = pi/4 -> sin(pi/4) = 1/sqrt(2)
    x_test = 0.25_dp * PI
    y_test = spline_eval(spl, x_test)
    y_exact = sin(x_test)
    call assert_close("Spline interior interpolation error < 1e-3", y_exact, y_test, 1.0e-3_dp)

    ! 一阶导数插值精度校验: d/dx sin(x) = cos(x)
    dy_test = spline_eval_deriv(spl, x_test)
    dy_exact = cos(x_test)
    call assert_close("Spline first derivative error < 5e-3", dy_exact, dy_test, 5.0e-3_dp)

    ! 二阶导数自然边界条件校验: y''(0) = 0, y''(pi) = 0
    call assert_close("Natural BC y''(0) = 0", 0.0_dp, spl%y2(1), 1.0e-12_dp)
    call assert_close("Natural BC y''(pi) = 0", 0.0_dp, spl%y2(n_pts), 1.0e-12_dp)

    call spline_clean(spl)

    ! 2. 从头算势能曲面格点投影与外推测试 (短程斥力 + 长程范德华)
    r_raw = [1.2_dp, 1.6_dp, 2.0_dp, 2.5_dp, 3.2_dp]
    v_raw = [0.15_dp, 0.02_dp, 0.00_dp, 0.015_dp, 0.035_dp]

    do i = 1, 20
        r_grid(i) = 0.8_dp + real(i - 1, dp) * (4.0_dp - 0.8_dp) / 19.0_dp
    end do

    call interpolate_pes_to_grid(r_raw, v_raw, r_grid, v_grid, v_inf=0.040_dp)

    ! 短程排斥区 V(0.8) > V(1.2)
    call assert_true("Short-range repulsive extrapolation V(0.8) > V(1.2)", v_grid(1) > v_raw(1))
    ! 长程渐近区 V(4.0) 趋近解离极限 V_inf (0.040 a.u.)
    call assert_close("Long-range asymptotic V(4.0) -> V_inf", 0.040_dp, v_grid(20), 0.005_dp)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Interpolation Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some interpolation tests failed!"
        stop 1
    else
        print '(A)', "SUCCESS: All interpolation tests passed."
    end if

contains

    subroutine assert_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        n_tests = n_tests + 1
        if (condition) then
            n_passed = n_passed + 1
            print '(A, A40, A)', " [PASS] ", name, ""
        else
            print '(A, A40, A)', " [FAIL] ", name, " (Condition was FALSE)"
        end if
    end subroutine assert_true

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

end program test_interpolation
