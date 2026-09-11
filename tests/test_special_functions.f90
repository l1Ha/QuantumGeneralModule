!> \brief 特殊函数与角动量耦合算法单元测试
program test_special_functions
    use general_module
    implicit none

    integer :: n_tests = 0
    integer :: n_passed = 0
    real(dp) :: tol = 1.0e-10_dp
    real(dp) :: res, exp_val

    print '(A)', "=================================================="
    print '(A)', "   GeneralModule Unit Tests: Special Functions    "
    print '(A)', "=================================================="

    ! 1. 普通勒让德多项式
    call assert_close("Legendre P0(0.5)", 1.0_dp, legendre_poly(0, 0.5_dp), tol)
    call assert_close("Legendre P1(0.5)", 0.5_dp, legendre_poly(1, 0.5_dp), tol)
    call assert_close("Legendre P2(0.5)", -0.125_dp, legendre_poly(2, 0.5_dp), tol)
    call assert_close("Legendre P3(0.5)", -0.4375_dp, legendre_poly(3, 0.5_dp), tol)

    ! 2. 缔合勒让德多项式 (含 Condon-Shortley 相位)
    ! P_1^1(0.6) = -sqrt(1 - 0.6^2) = -0.8
    call assert_close("Assoc Legendre P1^1(0.6)", -0.8_dp, assoc_legendre_poly(1, 1, 0.6_dp), tol)
    ! P_2^1(0.6) = -3 * 0.6 * 0.8 = -1.44
    call assert_close("Assoc Legendre P2^1(0.6)", -1.44_dp, assoc_legendre_poly(2, 1, 0.6_dp), tol)

    ! 3. Wigner 3j 符号
    ! (1 1 0; 0 0 0) = -1 / sqrt(3)
    exp_val = -1.0_dp / sqrt(3.0_dp)
    call assert_close("Wigner 3j (1 1 0; 0 0 0)", exp_val, wigner_3j(1, 1, 0, 0, 0, 0), tol)

    ! (1 1 2; 0 0 0) = sqrt(2/15)
    exp_val = sqrt(2.0_dp / 15.0_dp)
    call assert_close("Wigner 3j (1 1 2; 0 0 0)", exp_val, wigner_3j(1, 1, 2, 0, 0, 0), tol)

    ! 选择定则校验: m1+m2+m3 /= 0 应严格为 0
    call assert_close("Wigner 3j selection rule", 0.0_dp, wigner_3j(1, 1, 1, 1, 0, 0), tol)

    ! 4. Clebsch-Gordan 系数
    ! <1,0, 1,0 | 0,0> = -1 / sqrt(3)
    exp_val = -1.0_dp / sqrt(3.0_dp)
    call assert_close("CG <1,0,1,0|0,0>", exp_val, clebsch_gordan(1, 0, 1, 0, 0, 0), tol)

    ! 5. 转动偶极与取向矩阵元
    ! <0,0|cos(theta)|1,0> = 1/sqrt(3)
    exp_val = 1.0_dp / sqrt(3.0_dp)
    call assert_close("rot_matrix_cos_theta <0|cos|1>", exp_val, rot_matrix_cos_theta(0, 1, 0), tol)

    ! <0,0|cos^2(theta)|0,0> = 1/3
    exp_val = 1.0_dp / 3.0_dp
    call assert_close("rot_matrix_cos2_theta <0|cos^2|0>", exp_val, rot_matrix_cos2_theta(0, 0, 0), tol)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Special Functions Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some tests did not pass!"
        stop 1
    else
        print '(A)', "SUCCESS: All special functions tests passed."
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

end program test_special_functions
