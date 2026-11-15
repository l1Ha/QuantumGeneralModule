!> \brief Lanczos 短迭代 Krylov 传播器验收测试 (§14.1 路线图项)
!> \details 闭式基准验收:
!>          组 1: HO 谐振子本征基相干态 |alpha> 的解析演化
!>                exp(-iHt)|alpha> = |alpha e^{-i omega t}>, 与 Krylov 传播逐点对照;
!>          组 2: 二能级系统不变子空间精确终止 (m_eff = 2, beta_3 = 0) 与
!>                Rabi 闭式 P_e(t) 到机器精度;
!>          组 3: Lanczos 三对角矩阵结构验证 T = V^T H V 与基正交性 V^T V = I;
!>          组 4: 范数守恒 (长时, 多步) 与后验误差估计随 m 收敛;
!>          组 5: 与库内切比雪夫传播器交叉验证 (DVR 网格 + quartic 非谐).
program test_krylov_propagator
    use general_module
    implicit none

    ! 测试 5 的网格/势参数 (统一常量)
    integer, parameter :: ng = 64
    real(dp), parameter :: dxg = 0.15625_dp      ! 10.0 / 64
    real(dp), parameter :: x0_pkt = 1.0_dp, sig_pkt = 0.8_dp

    integer :: n_tests = 0, n_passed = 0
    integer :: stat, m_used
    real(dp) :: err_est

    call print_banner("GeneralModule: Krylov (Lanczos) Propagator Test", 60)

    ! ================================================================
    ! 组 1: HO 相干态解析演化 (闭式基准)
    ! ================================================================
    block
        integer, parameter :: nh = 48
        real(dp), parameter :: omega = 0.5_dp, alpha_coh = 0.7_dp
        real(dp) :: dt_h, t_fin, t_now, fid, fid_min
        complex(dp) :: psi_h(nh), psi_exact(nh)
        real(dp) :: c_n(0:nh - 1)
        integer :: n, steps

        ! 相干态 |alpha> = e^{-|a|^2/2} sum a^n/sqrt(n!) |n> (H = omega(a^dag a + 1/2))
        c_n = 0.0_dp
        do n = 0, nh - 1
            c_n(n) = exp(-alpha_coh**2 / 2.0_dp) * alpha_coh**n / sqrt(factorial(n))
        end do
        psi_h = cmplx(c_n, 0.0_dp, kind=dp)

        ! 解析演化 |psi(t)> = |alpha e^{-i omega t}>: 分量相位 e^{-i omega(n+1/2)t}
        dt_h = 0.1_dp
        t_fin = 2.0_dp * PI / omega   ! 一个完整周期
        steps = nint(t_fin / dt_h)
        t_now = 0.0_dp
        fid_min = 1.0_dp
        do n = 1, steps
            call krylov_propagate_step(psi_h, dt_h, 24, h_mult_ho, err_est, m_used, stat)
            t_now = t_now + dt_h
            call ho_coherent_exact(alpha_coh, omega, t_now, nh, psi_exact)
            fid = abs(sum(conjg(psi_exact) * psi_h))
            fid_min = min(fid_min, fid)
        end do
        n_tests = n_tests + 1
        if (stat == 0 .and. fid_min > 1.0_dp - 1.0e-10_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4)', " [PASS] HO coherent state: min fidelity vs closed form =", fid_min
        else
            print '(A, ES12.4)', " [FAIL] HO coherent state fidelity", fid_min
        end if

        ! <x>(t) = sqrt(2/omega) * alpha cos(omega t) 闭式对照 (x = sqrt(1/2w)(a + a^dag))
        block
            real(dp) :: x_expect, x_closed, x_err_max
            integer :: k
            complex(dp) :: psi_t(nh)
            psi_t = cmplx(c_n, 0.0_dp, kind=dp)
            t_now = 0.0_dp
            x_err_max = 0.0_dp
            do k = 1, 20
                call krylov_propagate_step(psi_t, 0.25_dp, 22, h_mult_ho, err_est, m_used, stat)
                t_now = t_now + 0.25_dp
                x_expect = sqrt(0.5_dp / omega) * real(sum(conjg(psi_t) * x_op_apply(psi_t)), dp)
                x_closed = sqrt(2.0_dp / omega) * alpha_coh * cos(omega * t_now)
                x_err_max = max(x_err_max, abs(x_expect - x_closed))
            end do
            n_tests = n_tests + 1
            if (x_err_max < 1.0e-9_dp) then
                n_passed = n_passed + 1
                print '(A, ES12.4)', " [PASS] HO <x>(t) vs closed form: max dev =", x_err_max
            else
                print '(A, ES12.4)', " [FAIL] HO <x>(t) max deviation", x_err_max
            end if
        end block
    end block

    ! ================================================================
    ! 组 2: 二能级不变子空间精确终止 + Rabi 闭式
    ! ================================================================
    block
        real(dp), parameter :: om_r = 0.6_dp, det_r = 0.3_dp
        real(dp) :: omega_r, dt_r, p_e, p_closed, err_max
        complex(dp) :: psi_r(2)
        integer :: k
        real(dp) :: alpha_lan(4), beta_lan(5), v_bas(2, 4)
        integer :: m_lan

        omega_r = sqrt(om_r**2 + det_r**2)
        ! H = [[0, om/2],[om/2, det]]: P_e(t) = (om^2/Omega^2) sin^2(Omega t/2)
        psi_r = cmplx([1.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], kind=dp)

        ! Krylov 空间在 d = 2 处精确终止 (不变子空间)
        call lanczos_tridiag(2, 4, [1.0_dp, 0.0_dp], h_mult_rabi, alpha_lan, beta_lan, &
                             v_bas, m_lan, stat)
        n_tests = n_tests + 1
        if (stat == 0 .and. m_lan == 2 .and. beta_lan(3) < 1.0e-14_dp) then
            n_passed = n_passed + 1
            print '(A, I2, A, ES10.2)', " [PASS] two-level invariant subspace: m_eff =", m_lan, &
                  ", beta_3 =", beta_lan(3)
        else
            print '(A, I2, ES12.4)', " [FAIL] two-level termination m_eff, beta_3 =", m_lan, beta_lan(3)
        end if

        ! Rabi 振荡到机器精度 (m = 8 请求, 实际 2 步后即精确)
        dt_r = 0.05_dp
        err_max = 0.0_dp
        do k = 1, 400
            call krylov_propagate_step(psi_r, dt_r, 8, h_mult_rabi, err_est, m_used, stat)
            p_e = abs(psi_r(2))**2
            p_closed = (om_r**2 / omega_r**2) * sin(0.5_dp * omega_r * k * dt_r)**2
            err_max = max(err_max, abs(p_e - p_closed))
        end do
        n_tests = n_tests + 1
        if (err_max < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4)', " [PASS] Rabi P_e(t) vs closed form: max err =", err_max
        else
            print '(A, ES12.4)', " [FAIL] Rabi max err", err_max
        end if
    end block

    ! ================================================================
    ! 组 3: Lanczos 三对角矩阵结构 (T = V^T H V) 与基正交性
    ! ================================================================
    block
        integer, parameter :: nq = 40
        real(dp) :: alpha_q(nq), beta_q(nq + 1), v_q(nq, nq)
        real(dp) :: hv(nq, nq), gram(nq, nq), t_recon(nq, nq)
        real(dp) :: err_orth, err_tri, hsum
        integer :: m_l, i, j

        ! 类随机实对称 H: H_ij = cos(i+j)/(1+|i-j|)
        do i = 1, nq
            do j = 1, nq
                hv(i, j) = cos(real(i + j, dp)) / real(1 + abs(i - j), dp)
            end do
        end do
        hv = 0.5_dp * (hv + transpose(hv))

        call lanczos_tridiag(nq, nq, [(1.0_dp, i=1, nq)], h_mult_dense, alpha_q, &
                             beta_q, v_q, m_l, stat)

        ! 正交性 V^T V = I
        gram = matmul(transpose(v_q(:, 1:m_l)), v_q(:, 1:m_l))
        err_orth = 0.0_dp
        do i = 1, m_l
            do j = 1, m_l
                err_orth = max(err_orth, abs(gram(i, j) - merge(1.0_dp, 0.0_dp, i == j)))
            end do
        end do

        ! 三对角重构: 非三对角元应为零, 对角/次对角元应等于 alpha/beta
        t_recon = matmul(transpose(v_q(:, 1:m_l)), matmul(hv, v_q(:, 1:m_l)))
        err_tri = 0.0_dp
        do i = 1, m_l
            do j = 1, m_l
                if (abs(i - j) > 1) err_tri = max(err_tri, abs(t_recon(i, j)))
            end do
        end do
        hsum = 0.0_dp
        do i = 1, m_l
            hsum = hsum + abs(t_recon(i, i) - alpha_q(i))
            if (i < m_l) hsum = hsum + abs(t_recon(i, i + 1) - beta_q(i + 1))
        end do

        n_tests = n_tests + 1
        if (err_orth < 1.0e-13_dp .and. err_tri < 1.0e-12_dp .and. hsum < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] Lanczos orthonormal basis + T = V^T H V: ", &
                  err_orth, " /", err_tri
        else
            print '(A, 3ES12.4)', " [FAIL] orthonormality/tridiagonal reconstruction", &
                  err_orth, err_tri, hsum
        end if
    end block

    ! ================================================================
    ! 组 4: 长时范数守恒 + 后验误差估计随 m 收敛
    ! ================================================================
    block
        integer, parameter :: nh4 = 48
        real(dp), parameter :: omega4 = 0.5_dp, alpha4 = 0.7_dp
        complex(dp) :: psi4(nh4)
        real(dp) :: nrm_drift, nrm0, errs(5)
        integer :: k, idx, n_dec
        real(dp) :: c4(0:nh4 - 1)
        integer :: n4

        do n4 = 0, nh4 - 1
            c4(n4) = exp(-alpha4**2 / 2.0_dp) * alpha4**n4 / sqrt(factorial(n4))
        end do

        ! 2000 步长时范数守恒
        psi4 = cmplx(c4, 0.0_dp, kind=dp)
        nrm0 = sqrt(real(sum(abs(psi4)**2), dp))
        do k = 1, 2000
            call krylov_propagate_step(psi4, 0.1_dp, 20, h_mult_ho, err_est, m_used, stat)
        end do
        nrm_drift = abs(sqrt(real(sum(abs(psi4)**2), dp)) - nrm0) / nrm0
        n_tests = n_tests + 1
        if (nrm_drift < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4)', " [PASS] norm conservation over 2000 steps: drift =", nrm_drift
        else
            print '(A, ES12.4)', " [FAIL] norm drift", nrm_drift
        end if

        ! 后验误差估计 vs 真实误差 (单步 dt = 5.0 大步长, 解析相干态为基准):
        ! 要求估计与真实误差同量级 (比值 ∈ [0.3, 10]) 且 err(m) 下降 > 3 个数量级
        n_tests = n_tests + 1
        block
            real(dp) :: true_err, ratio
            logical :: ok_ratio
            complex(dp) :: psi_ex4(nh4)
            ok_ratio = .true.
            do idx = 1, 5
                psi4 = cmplx(c4, 0.0_dp, kind=dp)
                call krylov_propagate_step(psi4, 5.0_dp, 4 * idx, h_mult_ho, err_est, m_used, stat)
                errs(idx) = err_est
                if (idx >= 2) then
                    call ho_coherent_exact(alpha4, omega4, 5.0_dp, nh4, psi_ex4)
                    true_err = sqrt(real(sum(abs(psi_ex4 - psi4)**2), dp))
                    ratio = err_est / max(true_err, 1.0e-300_dp)
                    if (ratio < 0.3_dp .or. ratio > 10.0_dp) ok_ratio = .false.
                end if
            end do
            n_dec = count(errs(2:5) < errs(1:4))
            if (ok_ratio .and. n_dec >= 3 .and. errs(5) < errs(1) * 1.0e-3_dp) then
                n_passed = n_passed + 1
                print '(A, I1, A, ES10.2, A, ES10.2)', " [PASS] a-posteriori err estimate tracks true error; drops =", &
                      n_dec, "/4, err(4) =", errs(1), " err(20) =", errs(5)
            else
                print '(A, L3, I1, A, 5ES11.3)', " [FAIL] err estimate ratio ok =", ok_ratio, &
                      " drops =", n_dec, " errs =", errs
            end if
        end block
    end block

    ! ================================================================
    ! 组 5: 与切比雪夫传播器交叉验证 (DVR 网格 + quartic 非谐)
    ! ================================================================
    block
        complex(dp) :: psi_k(ng), psi_c(ng), psi0_g(ng)
        real(dp) :: fid_x, e_min_c, e_max_c, v_edge, x_edge
        integer :: i2

        do i2 = 1, ng
            psi0_g(i2) = cmplx(exp(-(x_of(i2) - x0_pkt)**2 / (2.0_dp * sig_pkt**2)), &
                               0.0_dp, kind=dp)
        end do
        psi0_g = psi0_g / sqrt(real(sum(abs(psi0_g)**2), dp) * dxg)

        ! Gershgorin 谱界: 对角 1/dxg^2 + V, 非对角绝对值之和 1/dxg^2
        x_edge = 0.5_dp * real(ng, dp) * dxg
        v_edge = 0.5_dp * x_edge**2 + 0.01_dp * x_edge**4
        e_min_c = -1.0_dp
        e_max_c = v_edge + 2.0_dp / dxg**2 + 1.0_dp

        psi_k = psi0_g
        psi_c = psi0_g
        call krylov_propagate_step(psi_k, 0.02_dp, 24, h_mult_quartic, err_est, m_used, stat)
        call chebyshev_propagate_step(psi_c, 0.02_dp, e_min_c, e_max_c, h_mult_quartic_c, 40)
        fid_x = abs(sum(conjg(psi_c) * psi_k)) * dxg
        n_tests = n_tests + 1
        if (fid_x > 1.0_dp - 1.0e-11_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4)', " [PASS] Krylov vs Chebyshev cross-validation: 1 - F =", 1.0_dp - fid_x
        else
            print '(A, ES12.4)', " [FAIL] Krylov vs Chebyshev fidelity", fid_x
        end if
    end block

    print '(A, I2, A, I2, A)', "Krylov Propagator Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: Krylov propagator test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> n! (整数阶乘, 双精度)
    pure function factorial(n) result(f)
        integer, intent(in) :: n
        real(dp) :: f
        integer :: i
        f = 1.0_dp
        do i = 2, n
            f = f * real(i, dp)
        end do
    end function factorial

    !> DVR 网格坐标: x_i = (i - (ng+1)/2) * dxg
    pure function x_of(i) result(x)
        integer, intent(in) :: i
        real(dp) :: x
        x = (real(i, dp) - 0.5_dp * real(ng + 1, dp)) * dxg
    end function x_of

    !> HO 粒子数基: (H psi)_n = omega(n + 1/2) psi_n
    subroutine h_mult_ho(v_in, v_out)
        real(dp), intent(in) :: v_in(:)
        real(dp), intent(out) :: v_out(:)
        real(dp), parameter :: om = 0.5_dp
        integer :: n
        do n = 1, size(v_in)
            v_out(n) = om * (real(n - 1, dp) + 0.5_dp) * v_in(n)
        end do
    end subroutine h_mult_ho

    !> 相干态解析演化: c_n(t) = c_n(0) exp(-i omega (n+1/2) t)
    subroutine ho_coherent_exact(alp, om, t, nh, psi_ex)
        real(dp), intent(in) :: alp, om, t
        integer, intent(in) :: nh
        complex(dp), intent(out) :: psi_ex(nh)
        integer :: n
        do n = 0, nh - 1
            psi_ex(n + 1) = exp(-alp**2 / 2.0_dp) * alp**n / sqrt(factorial(n)) &
                            * exp(cmplx(0.0_dp, -om * (real(n, dp) + 0.5_dp) * t, kind=dp))
        end do
    end subroutine ho_coherent_exact

    !> 数基中施加 (a + a^dag) (用于 <x> 期望)
    function x_op_apply(v) result(xv)
        complex(dp), intent(in) :: v(:)
        complex(dp) :: xv(size(v))
        integer :: n, nh
        nh = size(v)
        xv = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do n = 1, nh - 1
            xv(n) = xv(n) + sqrt(real(n, dp)) * v(n + 1)
            xv(n + 1) = xv(n + 1) + sqrt(real(n, dp)) * v(n)
        end do
    end function x_op_apply

    !> Rabi 二能级: H = [[0, om/2],[om/2, det]] (实对称)
    subroutine h_mult_rabi(v_in, v_out)
        real(dp), intent(in) :: v_in(:)
        real(dp), intent(out) :: v_out(:)
        real(dp), parameter :: om = 0.6_dp, det = 0.3_dp
        v_out(1) = 0.5_dp * om * v_in(2)
        v_out(2) = 0.5_dp * om * v_in(1) + det * v_in(2)
    end subroutine h_mult_rabi

    !> 稠密实对称矩阵乘向量 (组 3)
    subroutine h_mult_dense(v_in, v_out)
        real(dp), intent(in) :: v_in(:)
        real(dp), intent(out) :: v_out(:)
        integer, parameter :: nq = 40
        integer :: i, j
        v_out = 0.0_dp
        do i = 1, nq
            do j = 1, nq
                v_out(i) = v_out(i) + cos(real(i + j, dp)) / real(1 + abs(i - j), dp) * v_in(j)
            end do
        end do
    end subroutine h_mult_dense

    !> DVR 非谐哈密顿量 (复数接口): H = -1/2 d^2/dx^2 + x^2/2 + 0.01 x^4
    subroutine h_mult_quartic_c(v_in, v_out)
        complex(dp), intent(in) :: v_in(:)
        complex(dp), intent(out) :: v_out(:)
        real(dp) :: kin_pref
        integer :: i
        kin_pref = 0.5_dp / dxg**2
        do i = 1, ng
            v_out(i) = (0.5_dp * x_of(i)**2 + 0.01_dp * x_of(i)**4) * v_in(i) &
                       - 2.0_dp * kin_pref * v_in(i)
            if (i > 1) v_out(i) = v_out(i) + kin_pref * v_in(i - 1)
            if (i < ng) v_out(i) = v_out(i) + kin_pref * v_in(i + 1)
        end do
    end subroutine h_mult_quartic_c

    !> 同上, 实数接口 (Krylov 用)
    subroutine h_mult_quartic(v_in, v_out)
        real(dp), intent(in) :: v_in(:)
        real(dp), intent(out) :: v_out(:)
        complex(dp) :: c_in(size(v_in)), c_out(size(v_in))
        c_in = cmplx(v_in, 0.0_dp, kind=dp)
        call h_mult_quartic_c(c_in, c_out)
        v_out = real(c_out, dp)
    end subroutine h_mult_quartic

end program test_krylov_propagator
