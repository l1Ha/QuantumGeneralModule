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
    call test_morse_anharmonic()
    call test_coriolis_vibration()

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

    !> Morse 非谐 (多原子拓展第四步): 单模精确谱对照 (分层容差: 顶部 Morse
    !> 能级尾部衰减率 kappa = sqrt(2(D-E)) 饱和, 盒截断误差只随盒宽收敛,
    !> 不随格点数收敛), 束缚能级计数, 非谐红移与间距收缩, 多模可分离谱。
    subroutine test_morse_anharmonic()
        real(dp) :: lam(3, 3), d_e(2), a_m(2), xh1(1), xh2(2)
        real(dp), allocatable :: e(:), e2(:)
        integer :: stat, i, n_bound_dvr
        real(dp) :: w0_1

        print *, ""
        print *, "--- 5. Morse anharmonic force field (exact spectrum) ---"
        lam = 0.0_dp

        ! (a) 单模: D = 0.1744, a = 0.2 -> w0 = 0.1181, 3 个束缚态
        xh1 = (/ 28.0_dp /)
        call poly_vib_morse_levels(1, (/ 0.1744_dp /), (/ 0.2_dp /), lam, xh1, 60, e, stat)
        w0_1 = 0.2_dp * sqrt(2.0_dp * 0.1744_dp)
        call assert_true(stat == 0, "1-mode Morse DVR solve succeeded")
        call assert_near(e(1), morse_bound_energy(0.0_dp, 0.1744_dp, 0.2_dp), 1.0e-10_dp, &
                         "Morse G(0) vs exact spectrum")
        call assert_near(e(2), morse_bound_energy(1.0_dp, 0.1744_dp, 0.2_dp), 1.0e-5_dp, &
                         "Morse G(1) vs exact spectrum")
        ! 顶部能级接近解离: 尾巴长 (见子程序 docstring), 容差按物理放宽
        call assert_near(e(3), morse_bound_energy(2.0_dp, 0.1744_dp, 0.2_dp), 5.0e-3_dp, &
                         "Morse G(2) vs exact (near-dissociation tail)")

        ! (b) 束缚能级计数: E < D_e 的 DVR 能级数 = morse_bound_count = 3
        n_bound_dvr = 0
        do i = 1, size(e)
            if (e(i) < 0.1744_dp) n_bound_dvr = n_bound_dvr + 1
        end do
        call assert_true(n_bound_dvr == morse_bound_count(0.1744_dp, 0.2_dp), &
                         "bound-level count = morse_bound_count")

        ! (c) 非谐红移 + 间距收缩: E1-E0 < w0 且 E2-E1 < E1-E0
        call assert_true(e(2) - e(1) < w0_1, "anharmonic red shift (E1-E0 < w0)")
        call assert_true(e(3) - e(2) < e(2) - e(1), "spacing contraction (E2-E1 < E1-E0)")

        ! (d) 双模可分离 (lam = 0): E(v1,v2) = G1(v1) + G2(v2) 精确。
        !     模 2 取深阱 (D=0.4) 使其 v=1 远离解离、第一激发由模 1 承担
        !     (模 1 的 v=1 尾巴在盒内的截断误差 ~2e-7, 决定该断言容差)。
        d_e = (/ 0.1744_dp, 0.4000_dp /)
        a_m = (/ 0.2_dp, 0.24_dp /)
        xh2 = (/ 28.0_dp, 20.0_dp /)
        call poly_vib_morse_levels(2, d_e, a_m, lam, xh2, 40, e, stat)
        call assert_near(e(1), morse_bound_energy(0.0_dp, d_e(1), a_m(1)) &
                             + morse_bound_energy(0.0_dp, d_e(2), a_m(2)), 1.0e-10_dp, &
                         "2-mode separable Morse ZPE = G1(0)+G2(0)")
        ! 第一激发 = 全谱排序后最低单量子 (排序陷阱的可分离 Morse 版)
        block
            real(dp) :: d1q, d2q
            d1q = morse_bound_energy(1.0_dp, d_e(1), a_m(1)) - morse_bound_energy(0.0_dp, d_e(1), a_m(1))
            d2q = morse_bound_energy(1.0_dp, d_e(2), a_m(2)) - morse_bound_energy(0.0_dp, d_e(2), a_m(2))
            call assert_near(e(2) - e(1), min(d1q, d2q), 5.0e-6_dp, &
                             "first excited = softest single quantum (sorted)")
            call assert_true(d2q > d1q, "mode-2 quantum above mode-1 (deep well)")
        end block

        ! (e) 可分离 ZPE 网格收敛 (盒误差对 n 不敏感, 见 docstring)
        call poly_vib_morse_levels(2, d_e, a_m, lam, xh2, 48, e2, stat)
        call assert_near(e2(1), e(1), 1.0e-10_dp, "separable Morse ZPE grid-converged")
    end subroutine test_morse_anharmonic

    !> Coriolis 振动角动量 (多原子拓展第五步): 简并模对 + L_z 耦合。
    !> 精确谱 E(v,l) = w(v+1) - g l; 验收: (1) g=0 与谐振乘积机制逐位一致;
    !> (2) 2 倍维实形式的每个能级精确二重 (H_c 与 H_c* 同谱);
    !> (3) ZPE (l=0) 与 g 无关且等于 w; (4) v=1 Pi 双重态 2w -/+ g;
    !> (5) v=2 Sigma/Delta 三重态 3w, 3w -/+ 2g。
    subroutine test_coriolis_vibration()
        real(dp), allocatable :: e0(:), eg(:), eh(:)
        real(dp) :: om2(2), lam2(2, 2), xh2(2), w, g
        integer :: stat, i, n_bad

        print *, ""
        print *, "--- 6. Coriolis vibrational angular momentum (degenerate pair) ---"
        w = 0.12_dp
        g = 0.03_dp

        ! (1) g = 0 回归: 与第三步谐振乘积机制逐位一致
        call poly_vib_coriolis_levels(w, 0.0_dp, 10.0_dp, 20, e0, stat)
        call assert_true(stat == 0, "coriolis g=0 solve succeeded")
        om2 = (/ w, w /); lam2 = 0.0_dp; xh2 = (/ 10.0_dp, 10.0_dp /)
        call poly_vib_product_levels(2, om2, lam2, xh2, 20, eh, stat)
        call assert_near(e0(1), eh(1), 1.0e-12_dp, "g=0 ZPE = harmonic product ZPE")
        call assert_near(e0(3), eh(2), 1.0e-12_dp, "g=0 first excited = harmonic product")

        ! (2) 2 倍维实形式: 每个物理能级精确出现两次
        n_bad = 0
        do i = 1, size(e0) - 1, 2
            if (abs(e0(i) - e0(i + 1)) > 1.0e-11_dp) n_bad = n_bad + 1
        end do
        call assert_true(n_bad == 0, "doubled-real spectrum exactly pairwise degenerate")

        ! (3) g = 0.03: ZPE (l=0) 与 g 无关且等于 w
        call poly_vib_coriolis_levels(w, g, 10.0_dp, 20, eg, stat)
        ! QL 特征解的绝对噪声随矩阵范数增长 (g != 0 引入 ~g*x*p 量级非对角元),
        ! 实测 ZPE 差 ~2e-8, 容差取 1e-6
        call assert_near(eg(1), e0(1), 1.0e-6_dp, "ZPE g-independent (l=0 unshifted)")
        call assert_near(eg(1), w, 2.0e-5_dp, "ZPE = omega (2-mode DVR accuracy)")

        ! (4) v=1 Pi 双重态: 2w -/+ g
        call assert_near(eg(3), 2.0_dp * w - g, 5.0e-4_dp, "v=1 doublet lower = 2w - g")
        call assert_near(eg(5), 2.0_dp * w + g, 5.0e-4_dp, "v=1 doublet upper = 2w + g")

        ! (5) v=2 Sigma/Delta: 3w 与 3w -/+ 2g 均在谱中 (排序后位置由 l 决定)
        call assert_near(eg(7), 3.0_dp * w - 2.0_dp * g, 2.0e-3_dp, "v=2 Delta lower = 3w - 2g")
        call assert_near(eg(9), 3.0_dp * w, 2.0e-3_dp, "v=2 Sigma = 3w")
    end subroutine test_coriolis_vibration

end program test_polyatomic_vib
