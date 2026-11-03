!> \brief 多原子 N 模振动模块自动化单元测试 (多原子拓展第三步)
!> \details 验证层次: (1) 分离谐振极限 (lam = 0) 的解析谱 sum w(v+1/2);
!>          (2) 双模双线性耦合的解析简正频率 ZPE 与第一激发对照;
!>          (3) 简并旋转不变性 (w1 = w2 时任意 lam 只旋转简正坐标);
!>          (4) 网格收敛 (n_pts 翻倍, ZPE 漂移 < 1e-6 a.u.)。
!> \author LiHao
program test_polyatomic_vib
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_polyatomic_vib
    implicit none

    integer :: total_tests = 0
    integer :: passed_tests = 0

    print *, "================================================================"
    print *, "   TEST SUITE: Polyatomic N-Mode Vibrations (Multi-Atomic Step 3)"
    print *, "================================================================"

    call test_separable_harmonic()
    call test_coupled2_analytic()
    call test_degenerate_rotation()
    call test_grid_convergence()

    print *, ""
    print *, "================================================================"
    write(*, '(A, I3, A, I3, A)') " Polyatomic Vibration Test Results: ", &
        passed_tests, "/", total_tests, " PASSED"
    print *, "================================================================"

    if (passed_tests /= total_tests) then
        error stop "Some tests failed in test_polyatomic_vib!"
    end if

contains

    subroutine assert_near(val, target, tol, name)
        real(dp), intent(in) :: val, target, tol
        character(len=*), intent(in) :: name
        call count_test(abs(val - target) <= tol)
        if (abs(val - target) <= tol) then
            write(*, '(A, A, A)') "  [PASS] ", trim(name), ""
        else
            write(*, '(A, A, A, ES12.4, A, ES12.4)') "  [FAIL] ", trim(name), &
                "  Got: ", val, " Expected: ", target
        end if
    end subroutine assert_near

    subroutine assert_true(cond, name)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: name
        call count_test(cond)
        if (cond) then
            write(*, '(A, A, A)') "  [PASS] ", trim(name), ""
        else
            write(*, '(A, A, A)') "  [FAIL] ", trim(name), ""
        end if
    end subroutine assert_true

    subroutine count_test(cond)
        logical, intent(in) :: cond
        total_tests = total_tests + 1
        if (cond) passed_tests = passed_tests + 1
    end subroutine count_test

    !> 分离谐振极限 (lam = 0): E = sum_k w_k (v_k + 1/2), 三模 ZPE 与单模激发
    subroutine test_separable_harmonic()
        real(dp) :: om(3), lam(3, 3), x_half(3)
        real(dp), allocatable :: e(:)
        integer :: stat

        print *, ""
        print *, "--- 1. Separable harmonic limit (lam = 0) ---"
        om = (/ 0.30_dp, 0.20_dp, 0.12_dp /)
        x_half = (/ 9.0_dp, 11.0_dp, 14.0_dp /)
        lam = 0.0_dp

        call poly_vib_product_levels(3, om, lam, x_half, 10, e, stat)
        call assert_true(stat == 0, "3-mode DVR solve succeeded")
        ! ZPE = (w1 + w2 + w3)/2 (sinc 变分, n=10 容差 1e-4)
        call assert_near(e(1), 0.5_dp * (om(1) + om(2) + om(3)), 1.0e-4_dp, &
                         "3-mode ZPE = sum w/2 (analytic)")
        ! 能级按升序排列: 第一激发 = 最软模量子 min(w) = w3
        call assert_near(e(2) - e(1), minval(om), 1.0e-4_dp, &
                         "First excited quantum = min(w) (sorted spectrum)")
        ! 第二低能级 = 次软单模量子 (w2 = 0.20 < 2 w3 = 0.24)
        call assert_near(e(3) - e(1), 0.20_dp, 1.0e-3_dp, &
                         "Second excited = mode-2 quantum (sorted)")
        ! 网格收敛: 加密后 ZPE 漂移
        block
            real(dp), allocatable :: e2(:)
            call poly_vib_product_levels(3, om, lam, x_half, 12, e2, stat)
            call assert_near(e2(1), e(1), 1.0e-4_dp, "ZPE grid-converged (10 -> 12 pts)")
        end block
    end subroutine test_separable_harmonic

    !> 双模双线性耦合: ZPE = (Om_lo + Om_hi)/2 (解析简正频率)
    subroutine test_coupled2_analytic()
        real(dp) :: om(2), lam(2, 2), x_half(2)
        real(dp), allocatable :: e(:)
        integer :: stat
        real(dp) :: om_lo, om_hi, zpe_exact

        print *, ""
        print *, "--- 2. Two-mode bilinear coupling vs analytic normal modes ---"
        om = (/ 0.28_dp, 0.19_dp /)
        x_half = (/ 10.0_dp, 13.0_dp /)
        lam(1, 1) = 0.0_dp
        lam(1, 2) = 0.045_dp
        lam(2, 2) = 0.0_dp

        call coupled2_normal_freqs(om(1), om(2), lam(1, 2), om_lo, om_hi)
        zpe_exact = 0.5_dp * (om_lo + om_hi)
        write(*, '(A, 2F10.6)') "  analytic normal freqs:", om_lo, om_hi

        call poly_vib_product_levels(2, om, lam, x_half, 20, e, stat)
        call assert_true(stat == 0, "2-mode coupled DVR solve succeeded")
        call assert_near(e(1), zpe_exact, 1.0e-6_dp, "Coupled ZPE = (Om_lo + Om_hi)/2")
        ! 第一激发 (低频简正模的一个量子): E - ZPE = Om_lo
        call assert_near(e(2) - e(1), om_lo, 1.0e-4_dp, &
                         "First excited = Om_lo (normal mode)")
        ! 高频简正模激发位置: 排序后需含 Om_hi 量子 (变分容差放宽至 1e-4)
        block
            integer :: k
            logical :: found_hi
            found_hi = .false.
            do k = 2, size(e)
                if (abs((e(k) - e(1)) - om_hi) < 1.0e-4_dp) found_hi = .true.
            end do
            call assert_true(found_hi, "Level at Om_hi present in spectrum")
        end block
    end subroutine test_coupled2_analytic

    !> 简并旋转不变性: w1 = w2 时动力学矩阵的本征矢旋转任意角度不变,
    !> 故耦合后的简正频率只依赖 lam 的正交不变量: Om^2 = w^2 -/+ lam
    subroutine test_degenerate_rotation()
        real(dp) :: om(2), lam(2, 2), x_half(2)
        real(dp), allocatable :: e(:)
        integer :: stat
        real(dp) :: w0, om_lo, om_hi

        print *, ""
        print *, "--- 3. Degenerate rotation invariance (w1 = w2) ---"
        w0 = 0.22_dp
        om = (/ w0, w0 /)
        x_half = (/ 12.0_dp, 12.0_dp /)
        lam(1, 1) = 0.0_dp
        lam(1, 2) = 0.03_dp
        lam(2, 2) = 0.0_dp

        call coupled2_normal_freqs(w0, w0, lam(1, 2), om_lo, om_hi)
        call poly_vib_product_levels(2, om, lam, x_half, 18, e, stat)
        call assert_near(e(1), 0.5_dp * (om_lo + om_hi), 1.0e-6_dp, &
                         "Degenerate coupled ZPE")
        call assert_near(om_hi - om_lo, &
                         sqrt(w0 * w0 + lam(1, 2)) - sqrt(w0 * w0 - lam(1, 2)), &
                         1.0e-8_dp, &
                         "Normal-frequency splitting (analytic)")
    end subroutine test_degenerate_rotation

    !> 网格收敛: n_pts 14 -> 20, 耦合 ZPE 漂移 < 1e-6 a.u.
    subroutine test_grid_convergence()
        real(dp) :: om(2), lam(2, 2), x_half(2)
        real(dp), allocatable :: e1(:), e2(:)
        integer :: stat

        print *, ""
        print *, "--- 4. Grid convergence of the coupled ZPE ---"
        om = (/ 0.26_dp, 0.17_dp /)
        x_half = (/ 10.0_dp, 13.0_dp /)
        lam = 0.0_dp
        lam(1, 2) = 0.02_dp

        call poly_vib_product_levels(2, om, lam, x_half, 14, e1, stat)
        call poly_vib_product_levels(2, om, lam, x_half, 20, e2, stat)
        call assert_near(e2(1), e1(1), 1.0e-6_dp, "Coupled ZPE grid-converged")
    end subroutine test_grid_convergence

end program test_polyatomic_vib
