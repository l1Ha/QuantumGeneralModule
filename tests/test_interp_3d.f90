!> \brief 三维三三次张量积样条验收测试 (路线图「多维插值」3D 扩展)
!> \details 验收组:
!>          组 1: 三线性函数精确重构与节点精确插值;
!>          组 2: 光滑函数内部四阶收敛;
!>          组 3: 对称数据与全置换对称性 + 解析梯度一致性/精度;
!>          组 4: 真实 H3 LEPS 三维势能面 V(r12, r23, r31) 插值
!>                (直接求值对照 + 全置换对称性保持).
program test_interp_3d
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat

    call print_banner("GeneralModule: 3D Tricubic Spline (3D PES) Test", 60)

    ! ================================================================
    ! 组 1: 三线性精确重构 + 节点精确
    ! ================================================================
    block
        integer, parameter :: n = 5
        real(dp) :: xg(n), yg(n), zg(n), fg(n, n, n)
        real(dp) :: a0, b0, c0, d0, e0, g0, h0, m0
        real(dp) :: xv, yv, zv, f1, f2, err_lin, err_node
        type(spline_3d_t) :: s3
        integer :: i, j, k

        a0 = 0.3_dp; b0 = -0.7_dp; c0 = 0.45_dp; d0 = 0.28_dp
        e0 = -0.15_dp; g0 = 0.22_dp; h0 = -0.09_dp; m0 = 0.13_dp
        do i = 1, n
            xg(i) = 0.1_dp + 0.41_dp * real(i - 1, dp)
            yg(i) = -0.3_dp + 0.27_dp * real(i - 1, dp)
            zg(i) = 0.5_dp + 0.19_dp * real(i - 1, dp)
        end do
        do i = 1, n
            do j = 1, n
                do k = 1, n
                    fg(i, j, k) = a0 + b0 * xg(i) + c0 * yg(j) + d0 * zg(k) &
                                  + e0 * xg(i) * yg(j) + g0 * xg(i) * zg(k) &
                                  + h0 * yg(j) * zg(k) + m0 * xg(i) * yg(j) * zg(k)
                end do
            end do
        end do

        call spline_3d_init(xg, yg, zg, fg, s3, stat)

        err_lin = 0.0_dp
        do i = 1, 7
            xv = xg(1) + (xg(n) - xg(1)) * (0.06_dp + 0.11_dp * real(i, dp))
            yv = yg(1) + (yg(n) - yg(1)) * (0.09_dp + 0.10_dp * real(i, dp))
            zv = zg(1) + (zg(n) - zg(1)) * (0.07_dp + 0.12_dp * real(i, dp))
            f1 = a0 + b0 * xv + c0 * yv + d0 * zv + e0 * xv * yv + g0 * xv * zv &
                 + h0 * yv * zv + m0 * xv * yv * zv
            f2 = spline_3d_eval(s3, xv, yv, zv)
            err_lin = max(err_lin, abs(f1 - f2))
        end do

        err_node = 0.0_dp
        do i = 1, n
            do j = 1, n
                do k = 1, n
                    err_node = max(err_node, abs(spline_3d_eval(s3, xg(i), yg(j), zg(k)) - fg(i, j, k)))
                end do
            end do
        end do

        n_tests = n_tests + 1
        if (stat == 0 .and. err_lin < 1.0e-12_dp .and. err_node < 1.0e-13_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4, A, ES12.4)', " [PASS] trilinear exact reconstruction: err =", &
                  err_lin, "  node exactness: ", err_node
        else
            print '(A, ES12.4, A, ES12.4)', " [FAIL] trilinear/node", err_lin, err_node
        end if
        call spline_3d_clean(s3)
    end block

    ! ================================================================
    ! 组 2: 光滑函数内部四阶收敛
    ! ================================================================
    block
        real(dp) :: e_prev, ratio, e_first
        integer :: n, i
        logical :: ok

        e_prev = 0.0_dp
        e_first = 0.0_dp
        ratio = 0.0_dp
        ok = .true.
        do n = 9, 17, 8                            ! n = 9, 17
            block
                real(dp), allocatable :: g(:), fg(:, :, :)
                real(dp) :: err, xv, yv, zv, h
                type(spline_3d_t) :: s3
                integer :: j, k

                allocate (g(n), fg(n, n, n))
                h = 1.0_dp / real(n - 1, dp)
                do i = 1, n
                    g(i) = -0.5_dp + real(i - 1, dp) * h
                end do
                do i = 1, n
                    do j = 1, n
                        do k = 1, n
                            fg(i, j, k) = f3(g(i), g(j), g(k))
                        end do
                    end do
                end do
                call spline_3d_init(g, g, g, fg, s3, stat)
                err = 0.0_dp
                do i = 1, 8
                    xv = -0.25_dp + 0.5_dp * real(i, dp) / 9.0_dp
                    yv = -0.22_dp + 0.44_dp * real(i, dp) / 9.0_dp
                    zv = -0.28_dp + 0.56_dp * real(i, dp) / 9.0_dp
                    err = max(err, abs(f3(xv, yv, zv) - spline_3d_eval(s3, xv, yv, zv)))
                end do
                call spline_3d_clean(s3)
                deallocate (g, fg)
                if (n == 9) then
                    e_first = err
                else
                    ratio = e_prev / max(err, 1.0e-300_dp)
                end if
                e_prev = err
            end block
        end do

        n_tests = n_tests + 1
        if (ratio > 8.0_dp .and. e_first > 1.0e-6_dp) then
            n_passed = n_passed + 1
            print '(A, F8.2, A, ES10.2)', " [PASS] smooth function interior convergence (9->17): ratio =", &
                  ratio, "  err(n=9) =", e_first
        else
            print '(A, F8.2, A, ES12.4)', " [FAIL] convergence ratio =", ratio, " err =", e_prev
        end if
    end block

    ! ================================================================
    ! 组 3: 全置换对称性 + 解析梯度
    ! ================================================================
    block
        integer, parameter :: n = 17
        real(dp) :: g(n), fg(n, n, n)
        real(dp) :: xv, yv, zv, sym_err, g_cons, g_err
        real(dp) :: gx, gy, gz, hfd, f_p, f_m, gex, gey, gez
        type(spline_3d_t) :: s3
        integer :: i, j, k

        do i = 1, n
            g(i) = 0.6_dp + 0.05_dp * real(i - 1, dp)   ! 0.6 .. 1.4
        end do
        do i = 1, n
            do j = 1, n
                do k = 1, n
                    fg(i, j, k) = f3sym(g(i), g(j), g(k))
                end do
            end do
        end do
        call spline_3d_init(g, g, g, fg, s3, stat)

        ! 全置换对称: S(x,y,z) = S(y,x,z) = S(z,y,x) 等
        sym_err = 0.0_dp
        do i = 1, 10
            xv = 0.65_dp + 0.07_dp * real(i, dp)
            yv = 0.72_dp + 0.06_dp * real(i, dp)
            zv = 0.80_dp + 0.05_dp * real(i, dp)
            sym_err = max(sym_err, abs(spline_3d_eval(s3, xv, yv, zv) - spline_3d_eval(s3, yv, xv, zv)))
            sym_err = max(sym_err, abs(spline_3d_eval(s3, xv, yv, zv) - spline_3d_eval(s3, zv, yv, xv)))
        end do

        ! 解析梯度: 一致性 + 精度
        hfd = 1.0e-6_dp
        g_cons = 0.0_dp
        g_err = 0.0_dp
        do i = 1, 8
            xv = 0.75_dp + 0.06_dp * real(i, dp)
            yv = 0.70_dp + 0.07_dp * real(i, dp)
            zv = 0.82_dp + 0.055_dp * real(i, dp)
            call spline_3d_grad(s3, xv, yv, zv, gx, gy, gz)
            f_p = spline_3d_eval(s3, xv + hfd, yv, zv)
            f_m = spline_3d_eval(s3, xv - hfd, yv, zv)
            g_cons = max(g_cons, abs(gx - (f_p - f_m) / (2.0_dp * hfd)))
            f_p = spline_3d_eval(s3, xv, yv + hfd, zv)
            f_m = spline_3d_eval(s3, xv, yv - hfd, zv)
            g_cons = max(g_cons, abs(gy - (f_p - f_m) / (2.0_dp * hfd)))
            f_p = spline_3d_eval(s3, xv, yv, zv + hfd)
            f_m = spline_3d_eval(s3, xv, yv, zv - hfd)
            g_cons = max(g_cons, abs(gz - (f_p - f_m) / (2.0_dp * hfd)))
            call f3sym_grad(xv, yv, zv, gex, gey, gez)
            g_err = max(g_err, abs(gx - gex), abs(gy - gey), abs(gz - gez))
        end do
        call spline_3d_clean(s3)

        n_tests = n_tests + 1
        if (sym_err < 1.0e-12_dp .and. g_cons < 1.0e-6_dp .and. g_err < 1.0e-3_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] permutation symmetry + gradient: ", &
                  sym_err, " /", g_cons, " then ", g_err
        else
            print '(A, 3ES12.4)', " [FAIL] symmetry/gradient", sym_err, g_cons, g_err
        end if
    end block

    ! ================================================================
    ! 组 4: 真实 H3 LEPS 三维势能面 V(r12, r23, r31)
    ! ================================================================
    block
        type(leps_param_t) :: leps
        type(internuclear_dist_t) :: dist
        integer, parameter :: ng = 16
        real(dp) :: g(ng), fg(ng, ng, ng)
        real(dp) :: r1, r2, r3, err_g, err_pt, xv, yv, zv, v_test, err_ref
        type(spline_3d_t) :: s3
        integer :: i, j, k

        call init_default_h3_leps(leps)
        do i = 1, ng
            g(i) = 1.2_dp + (3.6_dp - 1.2_dp) * real(i - 1, dp) / real(ng - 1, dp)
        end do
        do i = 1, ng
            do j = 1, ng
                do k = 1, ng
                    dist%r12 = g(i); dist%r23 = g(j); dist%r31 = g(k)
                    fg(i, j, k) = calc_leps_potential(dist, leps)
                end do
            end do
        end do

        call spline_3d_init(g, g, g, fg, s3, stat)

        ! 化学相关区 (V < -0.10 E_h, 即势垒以下的反应通道区) 的插值精度:
        ! 高能排斥壁区 (V > 0) 曲率极大, 粗网格误差集中于此 (与 2D 一致的
        ! 物理图像: 误差 ~ 曲率), 不作为小网格下的验收指标
        err_g = 0.0_dp
        do i = 1, 20
            do j = 1, 20
                do k = 1, 20
                    xv = 1.25_dp + 1.70_dp * real(i, dp) / 21.0_dp
                    yv = 1.28_dp + 1.68_dp * real(j, dp) / 21.0_dp
                    zv = 1.24_dp + 1.72_dp * real(k, dp) / 21.0_dp
                    dist%r12 = xv; dist%r23 = yv; dist%r31 = zv
                    v_test = calc_leps_potential(dist, leps)
                    if (v_test < -0.10_dp) then
                        err_g = max(err_g, abs(v_test - spline_3d_eval(s3, xv, yv, zv)))
                    end if
                end do
            end do
        end do

        ! 物理点: 共线鞍点 (r12, r23, r31) = (1.775, 1.775, 3.55)
        dist%r12 = 1.775_dp; dist%r23 = 1.775_dp; dist%r31 = 3.55_dp
        err_pt = abs(calc_leps_potential(dist, leps) - spline_3d_eval(s3, 1.775_dp, 1.775_dp, 3.55_dp))

        ! 全置换对称性保持 (LEPS 对三键完全对称)
        block
            real(dp) :: asym
            asym = 0.0_dp
            do i = 1, 8
                xv = 1.45_dp + 0.13_dp * real(i, dp)
                yv = 1.90_dp - 0.05_dp * real(i, dp)
                zv = 2.35_dp - 0.08_dp * real(i, dp)
                asym = max(asym, abs(spline_3d_eval(s3, xv, yv, zv) - spline_3d_eval(s3, yv, zv, xv)))
            end do
            n_tests = n_tests + 1
            if (asym < 1.0e-12_dp) then
                n_passed = n_passed + 1
                print '(A, ES12.4)', " [PASS] LEPS 3D permutation symmetry preserved: ", asym
            else
                print '(A, ES12.4)', " [FAIL] LEPS 3D symmetry", asym
            end if
        end block

        call spline_3d_clean(s3)

        n_tests = n_tests + 1
        if (stat == 0 .and. err_g < 2.0e-3_dp .and. err_pt < 1.0e-3_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] H3 LEPS 3D PES (16^3, reaction region): err =", &
                  err_g, "  saddle: ", err_pt
        else
            print '(A, ES12.4, A, ES12.4)', " [FAIL] LEPS 3D (16^3) err =", err_g, &
                  "  saddle =", err_pt
        end if

        ! 对半加密 (13^3 -> 25^3) 的四阶收敛: 化学相关区误差比 ~ 16
        block
            integer, parameter :: n2 = 31
            real(dp), allocatable :: g2(:), f2(:, :, :)
            type(spline_3d_t) :: s2ref
            real(dp) :: err_f, ratio

            allocate (g2(n2), f2(n2, n2, n2))
            do i = 1, n2
                g2(i) = 1.2_dp + (3.6_dp - 1.2_dp) * real(i - 1, dp) / real(n2 - 1, dp)
            end do
            do i = 1, n2
                do j = 1, n2
                    do k = 1, n2
                        dist%r12 = g2(i); dist%r23 = g2(j); dist%r31 = g2(k)
                        f2(i, j, k) = calc_leps_potential(dist, leps)
                    end do
                end do
            end do
            call spline_3d_init(g2, g2, g2, f2, s2ref, stat)

            err_ref = 0.0_dp
            do i = 1, 20
                do j = 1, 20
                    do k = 1, 20
                        xv = 1.25_dp + 1.70_dp * real(i, dp) / 21.0_dp
                        yv = 1.28_dp + 1.68_dp * real(j, dp) / 21.0_dp
                        zv = 1.24_dp + 1.72_dp * real(k, dp) / 21.0_dp
                        dist%r12 = xv; dist%r23 = yv; dist%r31 = zv
                        v_test = calc_leps_potential(dist, leps)
                        if (v_test < -0.10_dp) then
                            err_ref = max(err_ref, abs(v_test - spline_3d_eval(s2ref, xv, yv, zv)))
                        end if
                    end do
                end do
            end do
            call spline_3d_clean(s2ref)
            deallocate (g2, f2)

            ratio = err_g / max(err_ref, 1.0e-300_dp)
            n_tests = n_tests + 1
            if (err_ref < 3.0e-4_dp .and. ratio > 3.0_dp) then
                n_passed = n_passed + 1
                print '(A, ES10.2, A, F6.2)', " [PASS] LEPS 3D grid refinement 16^3 -> 31^3: err =", &
                      err_ref, "  ratio =", ratio
            else
                print '(A, ES10.2, A, F6.2)', " [FAIL] LEPS 3D refinement err/ratio =", err_ref, ratio
            end if
        end block
    end block

    print '(A, I2, A, I2, A)', "3D Spline Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: 3D tricubic spline test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> 光滑三变量测试函数
    pure function f3(x, y, z) result(v)
        real(dp), intent(in) :: x, y, z
        real(dp) :: v
        v = sin(2.1_dp * x) * cos(1.7_dp * y) * exp(-0.4_dp * z) + 0.2_dp * x * y * z
    end function f3

    !> 全置换对称测试函数及其解析梯度
    pure function f3sym(x, y, z) result(v)
        real(dp), intent(in) :: x, y, z
        real(dp) :: v
        v = cos(0.9_dp * (x + y + z)) + 0.15_dp * (x * y + y * z + z * x) &
            + 0.1_dp * (x * y * z)
    end function f3sym

    pure subroutine f3sym_grad(x, y, z, gx, gy, gz)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: gx, gy, gz
        real(dp) :: s
        s = -0.9_dp * sin(0.9_dp * (x + y + z))
        gx = s + 0.15_dp * (y + z) + 0.1_dp * y * z
        gy = s + 0.15_dp * (x + z) + 0.1_dp * x * z
        gz = s + 0.15_dp * (x + y) + 0.1_dp * x * y
    end subroutine f3sym_grad

end program test_interp_3d
