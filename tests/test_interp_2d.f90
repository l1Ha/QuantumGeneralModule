!> \brief 二维双三次张量积样条验收测试 (§15.10 路线图项: 多维插值)
!> \details 闭式与收敛验收:
!>          组 1: 双线性函数精确重构 (自然边界下机器精度) 与节点精确插值;
!>          组 2: 光滑函数内部四阶收敛 (误差 h -> h/2 比 ~ 16);
!>          组 3: 对称数据 S(x,y) = S(y,x) 与解析梯度一致性;
!>          组 4: 真实 H3 LEPS 势能面插值 (直接求值对照 + 鞍点区域精度).
program test_interp_2d
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat

    call print_banner("GeneralModule: 2D Bicubic Spline (PES) Test", 60)

    ! ================================================================
    ! 组 1: 双线性精确重构 + 节点插值精确
    ! ================================================================
    block
        integer, parameter :: nx = 6, ny = 7
        real(dp) :: xg(nx), yg(ny), zg(nx, ny)
        real(dp) :: f1, f2, err_lin, err_node, xv, yv
        real(dp) :: a0, b0, c0, d0
        type(spline_2d_t) :: s2
        integer :: i, j

        a0 = 0.7_dp; b0 = -1.3_dp; c0 = 0.45_dp; d0 = 0.28_dp
        do i = 1, nx
            xg(i) = 0.2_dp + 0.37_dp * real(i - 1, dp)
        end do
        do j = 1, ny
            yg(j) = -0.4_dp + 0.29_dp * real(j - 1, dp)
        end do
        do i = 1, nx
            do j = 1, ny
                zg(i, j) = a0 + b0 * xg(i) + c0 * yg(j) + d0 * xg(i) * yg(j)
            end do
        end do

        call spline_2d_init(xg, yg, zg, s2, stat)

        ! 内部随机点: 双线性函数应被自然边界三次样条精确重构
        err_lin = 0.0_dp
        do i = 1, 9
            xv = xg(1) + (xg(nx) - xg(1)) * (0.05_dp + 0.09_dp * real(i, dp))
            yv = yg(1) + (yg(ny) - yg(1)) * (0.08_dp + 0.07_dp * real(i, dp))
            f1 = a0 + b0 * xv + c0 * yv + d0 * xv * yv
            f2 = spline_2d_eval(s2, xv, yv)
            err_lin = max(err_lin, abs(f1 - f2))
        end do

        ! 节点处逐点精确
        err_node = 0.0_dp
        do i = 1, nx
            do j = 1, ny
                err_node = max(err_node, abs(spline_2d_eval(s2, xg(i), yg(j)) - zg(i, j)))
            end do
        end do

        n_tests = n_tests + 1
        if (stat == 0 .and. err_lin < 1.0e-12_dp .and. err_node < 1.0e-13_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4, A, ES12.4)', " [PASS] bilinear exact reconstruction: err =", &
                  err_lin, "  node exactness: ", err_node
        else
            print '(A, ES12.4, A, ES12.4)', " [FAIL] bilinear/node", err_lin, err_node
        end if
        call spline_2d_clean(s2)
    end block

    ! ================================================================
    ! 组 2: 光滑函数内部四阶收敛
    ! ================================================================
    block
        real(dp) :: e_prev, ratio_min, e_first
        integer :: n, i, j, kstep

        ratio_min = 1.0e30_dp
        e_first = 0.0_dp
        e_prev = 0.0_dp
        kstep = 0
        do n = 9, 33, 8                            ! n = 9, 17, 25, 33
            kstep = kstep + 1
            block
                real(dp) :: xg(n), yg(n), zg(n, n), err, xv, yv, exact
                type(spline_2d_t) :: s2
                real(dp) :: h

                h = 1.0_dp / real(n - 1, dp)
                do i = 1, n
                    xg(i) = -0.5_dp + real(i - 1, dp) * h
                    yg(i) = xg(i)
                end do
                do i = 1, n
                    do j = 1, n
                        zg(i, j) = f_smooth(xg(i), yg(j))
                    end do
                end do
                call spline_2d_init(xg, yg, zg, s2, stat)

                ! 仅取深内部区域 (距边界 >= 3h) 测四阶收敛
                err = 0.0_dp
                do i = 1, 12
                    xv = -0.30_dp + 0.60_dp * real(i, dp) / 13.0_dp
                    yv = -0.28_dp + 0.58_dp * real(i, dp) / 13.0_dp
                    exact = f_smooth(xv, yv)
                    err = max(err, abs(exact - spline_2d_eval(s2, xv, yv)))
                end do
                call spline_2d_clean(s2)
                if (kstep == 1) e_first = err
                if (kstep == 2) ratio_min = e_prev / max(err, 1.0e-300_dp)   ! 9 -> 17 严格对半
                e_prev = err
            end block
        end do

        ! n = 9 -> 17 严格减半: 4 阶收敛期望比 ~ 16; 后续 (非严格减半) 略低
        n_tests = n_tests + 1
        if (ratio_min > 6.0_dp .and. e_first / e_prev > 40.0_dp) then
            n_passed = n_passed + 1
            print '(A, F8.2, A, ES10.2)', " [PASS] smooth function interior convergence: min ratio =", &
                  ratio_min, "  err(n=9) =", e_first
        else
            print '(A, F8.2, A, ES12.4)', " [FAIL] convergence ratio =", ratio_min, " err =", e_prev
        end if
    end block

    ! ================================================================
    ! 组 3: 对称性 + 解析梯度一致性/精度
    ! ================================================================
    block
        integer, parameter :: n = 41
        real(dp) :: xg(n), yg(n), zg(n, n)
        real(dp) :: xv, yv, sym_err, g_fd_x, g_fd_y, hfd, xp, xm
        real(dp) :: gx, gy, gex, gey, g_err, g_cons
        type(spline_2d_t) :: s2
        integer :: i, j

        do i = 1, n
            xg(i) = 0.5_dp + 0.025_dp * real(i - 1, dp)  ! 0.5 .. 1.5
            yg(i) = xg(i)
        end do
        do i = 1, n
            do j = 1, n
                ! 对称但非平凡: 距对角线的函数
                zg(i, j) = f_sym(xg(i), yg(j))
            end do
        end do
        call spline_2d_init(xg, yg, zg, s2, stat)

        ! 对称性 S(x,y) = S(y,x)
        sym_err = 0.0_dp
        do i = 1, 15
            xv = 0.55_dp + 0.06_dp * real(i, dp)
            yv = 0.60_dp + 0.055_dp * real(i, dp)
            sym_err = max(sym_err, abs(spline_2d_eval(s2, xv, yv) - spline_2d_eval(s2, yv, xv)))
        end do

        ! 梯度: 解析 vs 中心差分 (一致性) 与 vs 精确梯度 (精度)
        hfd = 1.0e-6_dp
        g_cons = 0.0_dp
        g_err = 0.0_dp
        do i = 1, 12
            xv = 0.70_dp + 0.055_dp * real(i, dp)
            yv = 0.75_dp + 0.05_dp * real(i, dp)
            call spline_2d_grad(s2, xv, yv, gx, gy)
            xp = xv + hfd; xm = xv - hfd
            g_fd_x = (spline_2d_eval(s2, xp, yv) - spline_2d_eval(s2, xm, yv)) / (2.0_dp * hfd)
            xp = yv + hfd; xm = yv - hfd
            g_fd_y = (spline_2d_eval(s2, xv, xp) - spline_2d_eval(s2, xv, xm)) / (2.0_dp * hfd)
            g_cons = max(g_cons, abs(gx - g_fd_x), abs(gy - g_fd_y))
            call f_sym_grad(xv, yv, gex, gey)
            g_err = max(g_err, abs(gx - gex), abs(gy - gey))
        end do
        call spline_2d_clean(s2)

        n_tests = n_tests + 1
        if (sym_err < 1.0e-12_dp .and. g_cons < 1.0e-6_dp .and. g_err < 3.0e-4_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] symmetry + analytic gradient: ", &
                  sym_err, " /", g_cons, " vs FD /", g_err
        else
            print '(A, 3ES12.4)', " [FAIL] symmetry/gradient", sym_err, g_cons, g_err
        end if
    end block

    ! ================================================================
    ! 组 4: 真实 H3 LEPS 势能面插值
    ! ================================================================
    block
        type(leps_param_t) :: leps
        type(internuclear_dist_t) :: dist
        integer, parameter :: ng = 16
        real(dp) :: xg(ng), yg(ng), zg(ng, ng)
        real(dp) :: r1, r2, v_direct, v_spline, err_grid, err_ref, xv, yv
        real(dp) :: saddle_direct, saddle_spl
        type(spline_2d_t) :: s2
        integer :: i, j

        call init_default_h3_leps(leps)
        do i = 1, ng
            xg(i) = 1.2_dp + (3.0_dp - 1.2_dp) * real(i - 1, dp) / real(ng - 1, dp)
            yg(i) = xg(i)
        end do
        do i = 1, ng
            do j = 1, ng
                dist%r12 = xg(i)
                dist%r23 = yg(j)
                dist%r31 = xg(i) + yg(j)             ! 共线
                zg(i, j) = calc_leps_potential(dist, leps)
            end do
        end do

        call spline_2d_init(xg, yg, zg, s2, stat)

        ! 粗网格 (16x16) 上任意内部点的插值误差
        err_grid = 0.0_dp
        do i = 1, 12
            xv = 1.35_dp + 1.45_dp * real(i, dp) / 13.0_dp
            yv = 1.30_dp + 1.50_dp * real(i, dp) / 13.0_dp
            dist%r12 = xv
            dist%r23 = yv
            dist%r31 = xv + yv
            v_direct = calc_leps_potential(dist, leps)
            v_spline = spline_2d_eval(s2, xv, yv)
            err_grid = max(err_grid, abs(v_direct - v_spline))
        end do

        ! 鞍点 (1.775, 1.775) 区域 (已知文献值附近) 的插值精度
        xv = 1.775_dp
        yv = 1.775_dp
        dist%r12 = xv
        dist%r23 = yv
        dist%r31 = xv + yv
        saddle_direct = calc_leps_potential(dist, leps)
        saddle_spl = spline_2d_eval(s2, xv, yv)
        err_ref = abs(saddle_direct - saddle_spl)

        call spline_2d_clean(s2)

        ! 16x16 网格 (h = 0.12 a0) 的四阶样条误差预期 ~ 1e-4 E_h 量级
        n_tests = n_tests + 1
        if (stat == 0 .and. err_grid < 1.0e-3_dp .and. err_ref < 1.0e-3_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] H3 LEPS PES interpolation (16x16): max err =", &
                  err_grid, "  saddle region: ", err_ref
        else
            print '(A, ES12.4, A, ES12.4)', " [FAIL] LEPS interpolation", err_grid, err_ref
        end if

        ! 对称性: LEPS H3 等质量时 V(r1, r2) = V(r2, r1)
        call spline_2d_init(xg, yg, zg, s2, stat)
        block
            real(dp) :: asym
            asym = 0.0_dp
            do i = 1, 10
                xv = 1.4_dp + 0.15_dp * real(i, dp)
                yv = 2.0_dp - 0.06_dp * real(i, dp)
                asym = max(asym, abs(spline_2d_eval(s2, xv, yv) - spline_2d_eval(s2, yv, xv)))
            end do
            call spline_2d_clean(s2)
            n_tests = n_tests + 1
            if (asym < 1.0e-12_dp) then
                n_passed = n_passed + 1
                print '(A, ES12.4)', " [PASS] LEPS symmetry V(r1,r2) = V(r2,r1) preserved: ", asym
            else
                print '(A, ES12.4)', " [FAIL] LEPS symmetry", asym
            end if
        end block
    end block

    print '(A, I2, A, I2, A)', "2D Spline Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: 2D bicubic spline test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> 光滑测试函数 (C^\infty, 内部四阶收敛用)
    pure function f_smooth(x, y) result(v)
        real(dp), intent(in) :: x, y
        real(dp) :: v
        v = sin(2.1_dp * x) * cos(1.7_dp * y) + 0.3_dp * exp(-(x**2 + y**2))
    end function f_smooth

    !> 对称测试函数 (关于 x = y 严格对称) 及其解析梯度
    pure function f_sym(x, y) result(v)
        real(dp), intent(in) :: x, y
        real(dp) :: v
        v = cos(1.3_dp * (x + y)) * cos(0.9_dp * (x - y)) + 0.2_dp * (x * y)**2
    end function f_sym

    pure subroutine f_sym_grad(x, y, gx, gy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: gx, gy
        real(dp) :: c1, s1, c2, s2
        c1 = cos(1.3_dp * (x + y)); s1 = sin(1.3_dp * (x + y))
        c2 = cos(0.9_dp * (x - y)); s2 = sin(0.9_dp * (x - y))
        ! d/dx: -1.3 s1 c2 + c1 * (-0.9 s2); d/dy: -1.3 s1 c2 + c1 * (+0.9 s2)
        gx = -1.3_dp * s1 * c2 - 0.9_dp * c1 * s2 + 0.4_dp * x * y**2
        gy = -1.3_dp * s1 * c2 + 0.9_dp * c1 * s2 + 0.4_dp * x**2 * y
    end subroutine f_sym_grad

end program test_interp_2d
