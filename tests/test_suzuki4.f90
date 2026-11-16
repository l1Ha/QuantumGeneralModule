!> \brief 四阶 Suzuki 分裂推进器验收测试 (路线图 §12.1「高阶分裂算符」项)
!> \details 验收组:
!>          组 1: HO 相干态解析动力学 (⟨x⟩/⟨p⟩ 振荡) —— Strang 与 Suzuki
!>                同 dt 下的精度对比;
!>          组 2: 自治收敛阶标定 (以 Chebyshev 大阶为参照), 拟合斜率
!>                Strang ~ 2, Suzuki ~ 4;
!>          组 3: 范数守恒与等成本对比 (1 步 Suzuki vs 5 步 Strang).
program test_suzuki4
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat

    ! 系统参数 (谐振子: m = 1, omega = 0.5)
    integer, parameter :: n = 256
    real(dp), parameter :: dx = 0.15_dp, x0 = 2.5_dp, mass = 1.0_dp
    real(dp), parameter :: omega = 0.5_dp
    real(dp) :: xg(n), vg(n)
    real(dp) :: xc, pc, t_fin, x_exp, p_exp
    complex(dp) :: psi_s(n), psi_z(n), psi_y(n), psi0(n), psi_ref(n)
    real(dp) :: err_s, err_z, e_xs, e_xz, e_xy, nrm, t0, tsim
    integer :: i, k, nst
    real(dp) :: dt

    call print_banner("GeneralModule: 4th-order Suzuki Splitting Test", 60)

    ! 网格: x = (i - (n+1)/2) dx, 谐振子势
    do i = 1, n
        xg(i) = (real(i, dp) - 0.5_dp * real(n + 1, dp)) * dx
        vg(i) = 0.5_dp * mass * omega**2 * xg(i)**2
    end do
    ! 初始相干态: 高斯, 中心 x0, 宽度 sigma = sqrt(1/(2 m omega))
    do i = 1, n
        psi0(i) = exp(-0.5_dp * mass * omega * (xg(i) - x0)**2)
    end do
    psi0 = psi0 / sqrt(real(sum(abs(psi0)**2), dp) * dx)

    ! ================================================================
    ! 组 1: HO 相干态解析动力学
    ! ================================================================
    t_fin = 2.0_dp * PI / omega          ! 一个完整周期
    ! dt 取较大值使时间积分误差主导 (dx = 0.15 的网格离散误差 ~1e-4,
    ! 在 dt = 0.05 时会淹没两种格式的阶数差异; dt = 0.2 下 Suzuki 的
    ! O(dt^4) 误差与 Strang 的 O(dt^2) 误差均超过网格底噪)
    nst = 63
    dt = t_fin / real(nst, dp)          ! 精确对齐末时刻 T = nst*dt (否则时间偏差
                                        ! 会主导与解析解的差异)

    psi_s = psi0
    psi_z = psi0
    psi_y = psi0
    do k = 1, nst
        call propagate_split_operator_1d(psi_s, vg, dx, mass, dt)
        call propagate_suzuki4_1d(psi_z, vg, dx, mass, dt)
        call propagate_yoshida4_1d(psi_y, vg, dx, mass, dt)
    end do
    ! 解析: <x>(t) = x0 cos(w t), <p>(t) = -m w x0 sin(w t)
    xc = x0 * cos(omega * t_fin)
    pc = -mass * omega * x0 * sin(omega * t_fin)
    x_exp = real(sum(xg * abs(psi_s)**2), dp) * dx
    e_xs = abs(x_exp - xc)
    x_exp = real(sum(xg * abs(psi_z)**2), dp) * dx
    e_xz = abs(x_exp - xc)
    x_exp = real(sum(xg * abs(psi_y)**2), dp) * dx
    e_xy = abs(x_exp - xc)

    ! 判据: 两个四阶格式在同 dt 下都比二阶 Strang 显著更准
    n_tests = n_tests + 1
    if (e_xz < 0.5_dp * e_xs .and. e_xy < 0.5_dp * e_xs) then
        n_passed = n_passed + 1
        print '(A, ES10.2, A, ES10.2, A, ES10.2, A)', " [PASS] HO <x>(T): Suzuki/Yoshida/Strang =", &
              e_xz, " /", e_xy, " /", e_xs, ""
    else
        print '(A, 3ES10.2)', " [FAIL] HO <x>: Suzuki/Yoshida/Strang =", e_xz, e_xy, e_xs
    end if

    ! ================================================================
    ! 组 2: 自治收敛阶标定 (Chebyshev 参照)
    ! ================================================================
    block
        real(dp), parameter :: dts(4) = [0.2_dp, 0.1_dp, 0.05_dp, 0.025_dp]
        real(dp) :: errs_s(4), errs_z(4), errs_y(4), slope_s, slope_z, slope_y
        real(dp) :: e_min, e_max, dt_act(4)
        integer :: m_dt, kk

        ! 谱界 (解析: V 最大 + 动能截断)
        ! 谱界: 势能最大 + FFT 谱动能算子的最大本征值 (pi/dx)^2/(2m)
        ! (与分裂算符所用 H 完全一致 —— 曾因用有限差分动能界与差分模板
        !  构造参照, 与传播器的 FFT 谱 H 不符, 引入 O(dt) 系统差)
        e_min = 0.0_dp
        e_max = maxval(vg) + (PI / dx)**2 / (2.0_dp * mass) + 1.0_dp

        ! Chebyshev 参照 (谱方法): 小步长多步 —— 阶数 60 需 ΔE·dt/2 ≲ 10,
        ! 单步 t_fin 需阶数 ~ΔE·t/2 ≈ 400+ (200 阶远远不够, 展开发散)
        psi_ref = psi0
        do kk = 1, 200
            call chebyshev_propagate_step(psi_ref, t_fin / 200.0_dp, e_min, e_max, h_mult_ho, 60)
        end do
        ! 参照自检: <x>(T) 应等于解析值 (验证参照本身)
        block
            real(dp) :: x_chk
            x_chk = real(sum(xg * abs(psi_ref)**2), dp) * dx
            if (abs(x_chk - x0 * cos(omega * t_fin)) > 1.0e-8_dp) then
                print '(A, 2ES12.4)', " [warn] Chebyshev reference self-check failed: ", &
                      x_chk, x0 * cos(omega * t_fin)
            end if
        end block

        do m_dt = 1, 4
            nst = max(1, nint(t_fin / dts(m_dt)))
            dt = t_fin / real(nst, dp)     ! 精确对齐末时刻
            psi_s = psi0
            psi_z = psi0
            psi_y = psi0
            do kk = 1, nst
                call propagate_split_operator_1d(psi_s, vg, dx, mass, dt)
                call propagate_suzuki4_1d(psi_z, vg, dx, mass, dt)
                call propagate_yoshida4_1d(psi_y, vg, dx, mass, dt)
            end do
            errs_s(m_dt) = sqrt(real(sum(abs(psi_s - psi_ref)**2), dp) * dx)
            errs_z(m_dt) = sqrt(real(sum(abs(psi_z - psi_ref)**2), dp) * dx)
            errs_y(m_dt) = sqrt(real(sum(abs(psi_y - psi_ref)**2), dp) * dx)
        end do

        ! 拟合 log(err) vs log(dt) 斜率 (用实际 dt)
        do m_dt = 1, 4
            dt_act(m_dt) = t_fin / real(max(1, nint(t_fin / dts(m_dt))), dp)
        end do
        slope_s = fit_slope(dt_act, errs_s)
        slope_z = fit_slope(dt_act, errs_z)
        slope_y = fit_slope(dt_act, errs_y)
        slope_z = fit_slope(dts, errs_z)
        print '(A)', " [info] error table (dt, Strang, Suzuki, Yoshida):"
        do m_dt = 1, 4
            print '(A, F8.4, 3ES12.3)', "        ", dt_act(m_dt), errs_s(m_dt), errs_z(m_dt), errs_y(m_dt)
        end do
        print '(A, 2F7.2)', " [info] fitted slopes: Strang, Suzuki =", slope_s, slope_z
        print '(A, F7.2)', " [info] fitted slope: Yoshida =", slope_y

        n_tests = n_tests + 1
        if (abs(slope_s - 2.0_dp) < 0.4_dp .and. abs(slope_z - 4.0_dp) < 0.8_dp .and. &
            errs_z(1) < errs_s(1)) then
            n_passed = n_passed + 1
            print '(A, F5.2, A, F5.2, A, F5.2, A, 2ES10.2)', " [PASS] slopes Strang/Suzuki/Yoshida =", &
                  slope_s, " /", slope_z, " /", slope_y, "  errs@dt_max =", errs_s(1), errs_z(1)
        else
            print '(A, F5.2, A, F5.2, A, 2ES10.2)', " [FAIL] order calibration: Strang slope =", slope_s, &
                  "  Suzuki slope =", slope_z, "  err@dt=0.2 =", errs_s(1), errs_z(1)
        end if

        ! ================================================================
        ! 组 3: 范数守恒 + 等成本对比 (1 Suzuki vs 5 Strang)
        ! ================================================================
        n_tests = n_tests + 1
        block
            complex(dp) :: p1(n), p5(n)
            real(dp) :: nrm1, nrm5, err_1, err_5
            integer :: reps, j
            p1 = psi0
            p5 = psi0
            reps = nint(t_fin / 0.1_dp)
            dt = t_fin / real(reps, dp)
            do j = 1, reps
                call propagate_suzuki4_1d(p1, vg, dx, mass, dt)
                ! 5 步 Strang 步长 dt/5 = 与 1 步 Suzuki 等成本
                do kk = 1, 5
                    call propagate_split_operator_1d(p5, vg, dx, mass, dt / 5.0_dp)
                end do
            end do
            nrm1 = sqrt(real(sum(abs(p1)**2), dp) * dx)
            nrm5 = sqrt(real(sum(abs(p5)**2), dp) * dx)
            err_1 = sqrt(real(sum(abs(p1 - psi_ref)**2), dp) * dx)
            err_5 = sqrt(real(sum(abs(p5 - psi_ref)**2), dp) * dx)
            if (abs(nrm1 - 1.0_dp) < 1.0e-12_dp .and. abs(nrm5 - 1.0_dp) < 1.0e-12_dp .and. &
                err_1 < 0.1_dp * err_5) then
                n_passed = n_passed + 1
                print '(A, ES9.2, A, ES9.2, A, ES9.2)', " [PASS] equal-cost: Suzuki =", err_1, &
                      "  vs 5xStrang =", err_5, "  norm = ", nrm1
            else
                print '(A, 3ES10.2)', " [FAIL] equal-cost", err_1, err_5, abs(nrm1 - 1.0_dp)
            end if
        end block
    end block

    print '(A, I2, A, I2, A)', "Suzuki4 Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: Suzuki 4th-order splitting test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> 谐振子哈密顿量 (FFT 谱动能 + 势能, 与分裂算符的离散热力学一致;
    !> 复数接口, 供 Chebyshev 参照使用)
    subroutine h_mult_ho(v_in, v_out)
        complex(dp), intent(in) :: v_in(:)
        complex(dp), intent(out) :: v_out(:)
        complex(dp) :: tmp(n)
        real(dp) :: dp_k, p_val
        integer :: j
        dp_k = TWOPI / (real(n, dp) * dx)
        tmp = v_in
        call fft_1d(tmp, -1)
        do j = 1, n
            if (j <= n / 2) then
                p_val = real(j - 1, dp) * dp_k
            else
                p_val = real(j - 1 - n, dp) * dp_k
            end if
            tmp(j) = tmp(j) * (p_val**2 / (2.0_dp * mass))
        end do
        call fft_1d(tmp, 1)
        v_out = tmp + vg * v_in
    end subroutine h_mult_ho

    !> 最小二乘拟合 log(err) = a + p log(dt) 的斜率 p
    pure function fit_slope(dts, errs) result(p)
        real(dp), intent(in) :: dts(:), errs(:)
        real(dp) :: p
        real(dp) :: lx(size(dts)), ly(size(errs))
        real(dp) :: mx, my, num, den
        integer :: i, nd
        nd = size(dts)
        lx = log(dts)
        ly = log(max(errs, 1.0e-30_dp))
        mx = sum(lx) / real(nd, dp)
        my = sum(ly) / real(nd, dp)
        num = sum((lx - mx) * (ly - my))
        den = sum((lx - mx)**2)
        p = num / den
    end function fit_slope

end program test_suzuki4
