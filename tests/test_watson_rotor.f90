!> \brief 非对称陀螺 Watson 刚转子模块自动化单元测试 (多原子拓展第一步)
!> \details 四组验证: (1) 长对称陀螺极限 B=C 的解析谱; (2) 扁对称陀螺极限 A=B;
!>          (3) 球陀螺 A=B=C 的 (2J+1) 重简并; (4) 迹规则 sum E = (2J+1)J(J+1)(A+B+C)/3;
!>          (5) 由几何计算 H2O 转动常数 (垂直轴定理 + 与实验值 1% 内对照)。
!> \author LiHao
program test_watson_rotor
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_watson_rotor
    implicit none

    integer :: total_tests = 0
    integer :: passed_tests = 0

    print *, "================================================================"
    print *, "    TEST SUITE: Watson Rotor (Asymmetric Top, Polyatomic Step 1)"
    print *, "================================================================"

    call test_prolate_limit()
    call test_oblate_limit()
    call test_spherical_limit()
    call test_trace_rule()
    call test_water_constants()
    call test_centrifugal_distortion()
    call test_watson_d1_d2()

    print *, ""
    print *, "================================================================"
    write(*, '(A, I3, A, I3, A)') " Watson Rotor Test Results: ", passed_tests, "/", total_tests, " PASSED"
    print *, "================================================================"

    if (passed_tests /= total_tests) then
        error stop "Some tests failed in test_watson_rotor!"
    end if

contains

    subroutine assert_near(val, target, tol, name)
        real(dp), intent(in) :: val, target, tol
        character(len=*), intent(in) :: name
        call count_test(abs(val - target) <= tol)
        if (abs(val - target) > tol) then
            write(*, '(A, A, A, ES14.6, A, ES14.6)') "  [FAIL] ", trim(name), &
                "  Got: ", val, " Expected: ", target
        else
            write(*, '(A, A, A)') "  [PASS] ", trim(name), ""
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

    !> 长对称陀螺极限 B = C: E = B J(J+1) + (A-B) K^2。
    !> 升序槽位: K=0 占第 1 位, |K|=k>=1 简并对占第 2k 与 2k+1 位。
    subroutine test_prolate_limit()
        real(dp), parameter :: a_c = 27.42_dp, b_c = 14.59_dp
        real(dp), allocatable :: e(:), e_exp(:)
        integer :: stat, j_qnum, k, n_dim

        print *, ""
        print *, "--- 1. Prolate limit (B = C): E = B J(J+1) + (A-B) K^2 ---"
        do j_qnum = 0, 4
            n_dim = 2 * j_qnum + 1
            call asymmetric_top_levels(a_c, b_c, b_c, j_qnum, e, stat)
            allocate(e_exp(n_dim))
            e_exp(1) = b_c * real(j_qnum * (j_qnum + 1), dp)
            do k = 1, j_qnum
                e_exp(2 * k) = b_c * real(j_qnum * (j_qnum + 1), dp) &
                               + (a_c - b_c) * real(k * k, dp)
                e_exp(2 * k + 1) = e_exp(2 * k)
            end do
            write(*, '(A, I2, A, ES10.2, A, L1)') "  J =", j_qnum, &
                  "  max deviation = ", maxval(abs(e - e_exp)), &
                  "  exact=", maxval(abs(e - e_exp)) < 1.0e-10_dp
            call assert_true(stat == 0 .and. maxval(abs(e - e_exp)) < 1.0e-10_dp, &
                             "Prolate limit spectrum, J")
            deallocate(e, e_exp)
        end do
    end subroutine test_prolate_limit

    !> 扁对称陀螺极限 A = B: E = A J(J+1) - (A-C) K_c^2。
    !> 升序槽位: |K_c|=J 简并对占第 1,2 位, |K_c|=J-1 占第 3,4 位, ..., K_c=0 占末位。
    subroutine test_oblate_limit()
        real(dp), parameter :: a_c = 27.42_dp, c_c = 9.53_dp
        real(dp), allocatable :: e(:), e_exp(:)
        integer :: stat, j_qnum, kc, pos, n_dim

        print *, ""
        print *, "--- 2. Oblate limit (A = B): E = A J(J+1) - (A-C) K_c^2 ---"
        do j_qnum = 0, 4
            n_dim = 2 * j_qnum + 1
            call asymmetric_top_levels(a_c, a_c, c_c, j_qnum, e, stat)
            allocate(e_exp(n_dim))
            pos = 1
            do kc = j_qnum, 1, -1
                e_exp(pos) = a_c * real(j_qnum * (j_qnum + 1), dp) &
                             - (a_c - c_c) * real(kc * kc, dp)
                e_exp(pos + 1) = e_exp(pos)
                pos = pos + 2
            end do
            e_exp(n_dim) = a_c * real(j_qnum * (j_qnum + 1), dp)
            write(*, '(A, I2, A, ES10.2, A, L1)') "  J =", j_qnum, &
                  "  max deviation = ", maxval(abs(e - e_exp)), &
                  "  exact=", maxval(abs(e - e_exp)) < 1.0e-10_dp
            call assert_true(stat == 0 .and. maxval(abs(e - e_exp)) < 1.0e-10_dp, &
                             "Oblate limit spectrum, J")
            deallocate(e, e_exp)
        end do
    end subroutine test_oblate_limit

    !> 球陀螺极限 A = B = C: E = A J(J+1), (2J+1) 重简并
    subroutine test_spherical_limit()
        real(dp), parameter :: a_c = 10.0_dp
        real(dp), allocatable :: e(:)
        integer :: stat, j_qnum, k
        logical :: all_ok

        print *, ""
        print *, "--- 3. Spherical limit (A = B = C): E = A J(J+1) ---"
        do j_qnum = 1, 4
            call asymmetric_top_levels(a_c, a_c, a_c, j_qnum, e, stat)
            all_ok = (stat == 0) .and. (size(e) == 2 * j_qnum + 1)
            do k = 1, 2 * j_qnum + 1
                all_ok = all_ok .and. (abs(e(k) - a_c * real(j_qnum * (j_qnum + 1), dp)) < 1.0e-9_dp)
            end do
            write(*, '(A, I2, A, L1)') "  J =", j_qnum, "  spherical degeneracy exact: ", all_ok
            call assert_true(all_ok, "Spherical limit spectrum, J")
        end do
    end subroutine test_spherical_limit

    !> 迹规则: sum E = (2J+1) J(J+1) (A+B+C)/3 (三轴迹相等)
    subroutine test_trace_rule()
        real(dp), parameter :: a_c = 27.42_dp, b_c = 14.59_dp, c_c = 9.53_dp
        real(dp), allocatable :: e(:)
        integer :: stat, j_qnum, k
        real(dp) :: s_sum, target

        print *, ""
        print *, "--- 4. Trace rule: sum E = (2J+1) J(J+1) (A+B+C)/3 ---"
        do j_qnum = 1, 8
            call asymmetric_top_levels(a_c, b_c, c_c, j_qnum, e, stat)
            s_sum = sum(e)
            target = real(2 * j_qnum + 1, dp) * real(j_qnum * (j_qnum + 1), dp) &
                     * (a_c + b_c + c_c) / 3.0_dp
            write(*, '(A, I2, A, F14.8, A, F14.8, A, L1)') "  J =", j_qnum, &
                  "  sum(E)=", s_sum, "  target=", target, &
                  "  ok=", abs(s_sum - target) < 1.0e-8_dp
            call assert_near(s_sum, target, 1.0e-8_dp, "Trace rule, J")
        end do
    end subroutine test_trace_rule

    !> 水分子: 由几何 (r_OH=0.9572 A, 角 104.52 deg) 计算 A,B,C,
    !> 垂直轴定理自检 + 与实验值 27.33/14.57/9.49 cm^-1 在 1% 内对照
    subroutine test_water_constants()
        real(dp) :: pos(3, 3), masses(3)
        real(dp) :: a_c, b_c, c_c, i_p(3), com(3), r_oh, half_ang, x_h, y_h
        integer :: stat

        print *, ""
        print *, "--- 5. H2O rotational constants from geometry ---"
        r_oh = 0.9572_dp          ! angstrom
        half_ang = 0.5_dp * 104.52_dp * acos(-1.0_dp) / 180.0_dp
        x_h = r_oh * sin(half_ang)
        y_h = r_oh * cos(half_ang)
        ! O 在原点, 两个 H 关于 y 轴 (C2 平分轴) 对称
        pos(1, 1) = 0.0_dp;      pos(2, 1) = 0.0_dp;      pos(3, 1) = 0.0_dp
        pos(1, 2) = x_h;         pos(2, 2) = y_h;         pos(3, 2) = 0.0_dp
        pos(1, 3) = -x_h;        pos(2, 3) = y_h;         pos(3, 3) = 0.0_dp
        masses = (/ 15.9949_dp, 1.00783_dp, 1.00783_dp /)

        call inertia_constants_from_geometry(pos, masses, a_c, b_c, c_c, i_p, com, stat)
        print '(A, 3F10.4)', "  principal moments (amu A^2):", i_p
        print '(A, 3F10.4)', "  rotational constants (cm^-1):", a_c, b_c, c_c
        print '(A, 3F10.6)', "  center of mass:", com

        call assert_true(stat == 0, "H2O geometry conversion succeeded")
        ! 平面分子垂直轴定理: I_c = I_a + I_b
        call assert_near(i_p(3), i_p(1) + i_p(2), 1.0e-10_dp, "Perpendicular axis theorem")
        ! 与实验基态值 (A=27.33, B=14.57, C=9.49 cm^-1) 在 1% 内对照
        call assert_near(a_c, 27.33_dp, 0.30_dp, "H2O A constant vs experiment (1%)")
        call assert_near(b_c, 14.57_dp, 0.15_dp, "H2O B constant vs experiment (1%)")
        call assert_near(c_c, 9.49_dp, 0.10_dp, "H2O C constant vs experiment (1%)")
        ! 非对称性: A > B > C, J=1 的三条能级分裂
        block
            real(dp), allocatable :: e(:)
            integer :: jj
            call asymmetric_top_levels(a_c, b_c, c_c, 1, e, stat)
            call assert_true(stat == 0 .and. abs(e(1) - e(2)) > 1.0e-3_dp &
                             .and. abs(e(2) - e(3)) > 1.0e-3_dp, &
                             "J=1 triplet is split (asymmetric)")
            jj = 0
        end block
    end subroutine test_water_constants

    !> 离心畸变 (A reduction 对角四次项): (a) 刚体回归; (b) 球陀螺 D_J 解析;
    !> (c) 含四次项的迹规则 (d1 = d2 = 0); (d) 长对称陀螺 K=0 槽位的 H2 WDC 对照
    subroutine test_centrifugal_distortion()
        real(dp), parameter :: a_c = 27.42_dp, b_c = 14.59_dp, c_c = 9.53_dp
        real(dp), parameter :: d_j = 2.76e-3_dp, d_jk = 1.0e-2_dp, d_k = 4.0e-2_dp
        real(dp), allocatable :: e_rig(:), e_cd(:)
        integer :: stat, j_qnum
        real(dp) :: s_cd, s_rig, s4_trace, target4

        print *, ""
        print *, "--- 6. Centrifugal distortion (A reduction diagonal quartics) ---"

        ! (a) 刚体回归: 全零畸变 -> 与 rigid 例程逐位一致
        do j_qnum = 0, 5
            call asymmetric_top_levels(a_c, b_c, c_c, j_qnum, e_rig, stat)
            call asymmetric_top_levels_cd(a_c, b_c, c_c, 0.0_dp, 0.0_dp, 0.0_dp, j_qnum, e_cd, stat)
            call assert_true(size(e_cd) == size(e_rig) .and. &
                             maxval(abs(e_cd - e_rig)) < 1.0e-12_dp, &
                             "CD zero-regression vs rigid, J")
            write(*, '(A, I2, A, L1)') "  J =", j_qnum, "  zero-distortion regression: ", &
                  size(e_cd) == size(e_rig) .and. maxval(abs(e_cd - e_rig)) < 1.0e-12_dp
        end do

        ! (b) 球陀螺 + D_J: E = [A - D_J J(J+1)] J(J+1) (解析)
        j_qnum = 6
        call asymmetric_top_levels_cd(a_c, a_c, a_c, d_j, 0.0_dp, 0.0_dp, j_qnum, e_cd, stat)
        block
            real(dp) :: jj
            jj = real(j_qnum * (j_qnum + 1), dp)
            call assert_near(e_cd(1), (a_c - d_j * jj) * jj, 1.0e-9_dp, &
                             "Spherical top + D_J analytic (J=6 lowest)")
        end block

        ! (c) 迹规则扩展: sum E = (2J+1)J(J+1)(A+B+C)/3 + D_J (2J+1)[J(J+1)]^2
        !     + D_JK (2J+1)[J(J+1)]^2/3 + D_K (2J+1) J(J+1)(3J^2+3J-1)/15
        j_qnum = 7
        call asymmetric_top_levels_cd(a_c, b_c, c_c, d_j, d_jk, d_k, j_qnum, e_cd, stat)
        block
            real(dp) :: jj, s_sum, t_trace
            jj = real(j_qnum * (j_qnum + 1), dp)
            s_sum = sum(e_cd)
            t_trace = real(2 * j_qnum + 1, dp) * jj * (a_c + b_c + c_c) / 3.0_dp &
                      - d_j * real(2 * j_qnum + 1, dp) * jj * jj &
                      - d_jk * real(2 * j_qnum + 1, dp) * jj * jj / 3.0_dp &
                      - d_k * real(2 * j_qnum + 1, dp) * jj &
                        * (3.0_dp * real(j_qnum * j_qnum, dp) + 3.0_dp * real(j_qnum, dp) - 1.0_dp) / 15.0_dp
            call assert_near(s_sum, t_trace, 1.0e-8_dp, "Trace rule with quartics, J")
        end block

        ! (d) 长对称陀螺 (B=C) 的 K=0 槽位: B J(J+1) - D_J [J(J+1)]^2
        !     以书中 H2 数值 (B0 = 60.18 cm^-1, D_J = 4B0^3/w0^2 = 0.04895 cm^-1) 对照
        block
            real(dp), parameter :: b_h2 = 60.18_dp, d_h2 = 4.895e-2_dp
            real(dp) :: jj
            j_qnum = 5
            call asymmetric_top_levels_cd(b_h2, b_h2, b_h2, d_h2, 0.0_dp, 0.0_dp, j_qnum, e_cd, stat)
            jj = real(j_qnum * (j_qnum + 1), dp)
            call assert_near(e_cd(1), b_h2 * jj - d_h2 * jj * jj, 1.0e-8_dp, &
                             "H2 WDC diatomic limit (K=0 slot, J=5)")
        end block
    end subroutine test_centrifugal_distortion

    !> Watson d1/d2 非对角畸变 (多原子拓展第八步): Delta K = +/-2/4 梯级项。
    !> 迹定理: d1/d2 项纯非对角 (Tr = 0), sum(E) 与 d1/d2 无关 (精确, 一切阶)
    !> —— A/S reduction 在求和规则下解耦, 使 D 常数拟合对 reduction 选择稳健。
    !> 近简并 ±K 对的 d1 位移是一阶的 (简并微扰), 非近简并区为二阶。
    subroutine test_watson_d1_d2()
        real(dp), allocatable :: e0(:), e1(:), e2(:)
        real(dp) :: a_c, b_c, c_c, dj, djk, dk, s0, s1, s2, ratio
        integer :: stat, jj

        print *, ""
        print *, "--- 7. Watson d1/d2 off-diagonal quartics (trace theorem) ---"
        a_c = 2.0_dp; b_c = 0.9_dp; c_c = 0.6_dp
        dj = 2.76e-4_dp; djk = 1.0e-3_dp; dk = 4.0e-3_dp

        ! (1) 迹不变: sum(E) 与 d1/d2 无关 (精确, J = 3..8)
        do jj = 3, 8
            call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, jj, e0, stat)
            call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, jj, e1, stat, &
                                          d1_in=2.0e-6_dp, d2_in=5.0e-7_dp)
            call assert_near(sum(e1) - sum(e0), 0.0_dp, 1.0e-8_dp, &
                             "trace theorem, sum(E) d1/d2-independent, J")
        end do

        ! (2) d2 单独: 迹不变
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 6, e0, stat)
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 6, e1, stat, &
                                      d1_in=0.0_dp, d2_in=5.0e-7_dp)
        call assert_near(sum(e1) - sum(e0), 0.0_dp, 1.0e-8_dp, &
                         "trace theorem, d2-only")

        ! (3) 零回归: d1 = d2 = 0 与既有对角版本逐位一致
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 8, e0, stat)
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 8, e1, stat, &
                                      d1_in=0.0_dp, d2_in=0.0_dp)
        call assert_near(maxval(abs(e0 - e1)), 0.0_dp, 1.0e-12_dp, &
                         "d1=d2=0 regression vs diagonal-only version")

        ! (4) 近简并 ±K 对的 d1 位移是一阶的 (简并微扰): d1 -> 2d1 时
        !     max|E shift| 比值 = 2 (非近简并区应为 4 = 二阶)
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 8, e0, stat)
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 8, e1, stat, &
                                      d1_in=1.0e-8_dp)
        call asymmetric_top_levels_cd(a_c, b_c, c_c, dj, djk, dk, 8, e2, stat, &
                                      d1_in=2.0e-8_dp)
        s1 = maxval(abs(e1 - e0))
        s2 = maxval(abs(e2 - e0))
        ratio = s2 / s1
        call assert_near(ratio, 2.0_dp, 0.1_dp, &
                         "near-degenerate +/-K pair shifts are first order in d1")
    end subroutine test_watson_d1_d2

end program test_watson_rotor
