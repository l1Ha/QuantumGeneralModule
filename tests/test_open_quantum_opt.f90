!> \brief 开放量子系统 (Lindblad) 与量子最优控制 (Krotov) 单元测试
program test_open_quantum_opt
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    complex(dp) :: rho_pure(2, 2), rho_mixed(2, 2), rho_t(2, 2)
    complex(dp) :: h_2lvl(2, 2), l_relax(2, 2, 1), l_dephase(2, 2, 1)
    real(dp) :: gamma_r(1), gamma_d(1), entropy, coherence
    real(dp) :: dt
    integer :: step, nt

    ! 最优控制变量
    real(dp) :: h_diag(2), dip_m(2, 2)
    complex(dp) :: psi_ini(2), phi_tgt(2)
    real(dp), allocatable :: e_field(:)
    real(dp) :: fid_final
    integer :: actual_iters, nt_opt

    print '(A)', "=================================================="
    print '(A)', " GeneralModule Unit Tests: Open Quantum & OCT     "
    print '(A)', "=================================================="

    ! ----------------------------------------------------
    ! 1. 开放量子系统基底与状态诊断测试
    ! ----------------------------------------------------
    ! 纯态 rho = |1><1|
    rho_pure = (0.0_dp, 0.0_dp)
    rho_pure(1, 1) = (1.0_dp, 0.0_dp)
    call assert_close("Pure state purity = 1.0", 1.0_dp, calculate_quantum_purity(rho_pure), 1.0e-12_dp)

    call calculate_von_neumann_entropy(rho_pure, entropy)
    call assert_close("Pure state von Neumann entropy = 0.0", 0.0_dp, entropy, 1.0e-10_dp)

    ! 最大混合态 rho = 0.5*I
    rho_mixed = (0.0_dp, 0.0_dp)
    rho_mixed(1, 1) = (0.5_dp, 0.0_dp)
    rho_mixed(2, 2) = (0.5_dp, 0.0_dp)
    call assert_close("2-level mixed state purity = 0.5", 0.5_dp, calculate_quantum_purity(rho_mixed), 1.0e-12_dp)

    call calculate_von_neumann_entropy(rho_mixed, entropy)
    call assert_close("2-level mixed state entropy = ln(2)", log(2.0_dp), entropy, 1.0e-5_dp)

    ! ----------------------------------------------------
    ! 2. Lindblad 自发跃迁衰变演化动力学测试 (|2> -> |1>)
    ! ----------------------------------------------------
    h_2lvl = (0.0_dp, 0.0_dp)
    h_2lvl(1, 1) = (0.0_dp, 0.0_dp)
    h_2lvl(2, 2) = (0.1_dp, 0.0_dp)  ! 能隙 0.1 a.u.

    ! 自发辐射跃迁算符 L = |1><2|
    call create_relaxation_jump_op(2, 1, 2, l_relax(:, :, 1))
    gamma_r(1) = 0.05_dp   ! 衰减率

    ! 初态处于激发态 |2>
    rho_t = (0.0_dp, 0.0_dp)
    rho_t(2, 2) = (1.0_dp, 0.0_dp)

    dt = 0.2_dp
    nt = 500  ! 演化 100 a.u. 时间 (对应 gamma * t = 5, 衰变 > 99%)
    do step = 1, nt
        call rk4_lindblad_step(rho_t, h_2lvl, l_relax, gamma_r, dt)
    end do

    call assert_close("Decay conservation Tr(rho) = 1.0", 1.0_dp, real(rho_t(1,1) + rho_t(2,2), dp), 1.0e-10_dp)
    call assert_true("Decayed into ground state rho(1,1) > 0.99", real(rho_t(1, 1), dp) > 0.99_dp)

    ! ----------------------------------------------------
    ! 3. Lindblad 纯退相位演化测试 (T_2*)
    ! ----------------------------------------------------
    ! 初态为等权重叠加态 (|1> + |2>) / sqrt(2)
    rho_t = cmplx(0.5_dp, 0.0_dp, kind=dp)
    call create_dephasing_jump_op(2, 1, l_dephase(:, :, 1))
    gamma_d(1) = 0.10_dp

    do step = 1, 500
        call rk4_lindblad_step(rho_t, h_2lvl, l_dephase, gamma_d, dt)
    end do

    ! 退相位保持布居不变 rho(1,1) = 0.5, rho(2,2) = 0.5
    call assert_close("Dephasing preserves populations rho(1,1)=0.5", 0.5_dp, real(rho_t(1, 1), dp), 1.0e-8_dp)
    ! 非对角相干度指数衰减趋近于 0
    coherence = calculate_quantum_coherence(rho_t)
    call assert_true("Coherence decayed by pure dephasing (< 0.05)", coherence < 0.05_dp)

    ! ----------------------------------------------------
    ! 4. 量子最优控制理论 (Krotov 算法逆向设计脉冲) 测试
    ! ----------------------------------------------------
    ! 两能级系统: 目标将初态 |1> 完全激发到目标态 |2>
    h_diag(1) = 0.0_dp
    h_diag(2) = 0.05_dp  ! 0.05 a.u.
    dip_m = 0.0_dp
    dip_m(1, 2) = 1.0_dp
    dip_m(2, 1) = 1.0_dp

    psi_ini = [ (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp) ]
    phi_tgt = [ (0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp) ]

    nt_opt = 500
    dt = 0.5_dp * FS2AU
    allocate(e_field(nt_opt))
    do step = 1, nt_opt
        e_field(step) = 0.005_dp * sin(0.05_dp * real(step-1, dp) * dt)
    end do

    call oct_optimize_pulse(n_states=2, nt=nt_opt, dt=dt, h_diag=h_diag, &
                            dip_mat=dip_m, psi_ini=psi_ini, phi_tgt=phi_tgt, &
                            alpha_0=50.0_dp, t_ramp=10.0_dp*FS2AU, &
                            max_iters=25, target_fid=0.90_dp, &
                            e_field=e_field, final_fid=fid_final, &
                            actual_iters=actual_iters)

    call assert_true("Krotov OCT final fidelity > 90%", fid_final > 0.90_dp)
    call assert_true("Krotov actual iterations <= 25", actual_iters <= 25)
    deallocate(e_field)

    ! ----------------------------------------------------
    ! 5. Lindblad Liouvillian 谱方法 (开放系统谱方法, 热化二能级)
    ! ----------------------------------------------------
    ! 精确谱 {0, -Gamma, -(Gamma/2+2gamma_phi) -/+ i delta}; 阻尼基重构
    ! 与 RK4 推进两条独立方法交叉验证。
    block
        complex(dp), allocatable :: rho_spec(:, :), rho_rk(:, :), h_bl(:, :), &
                                    l_ops_sp(:, :, :)
        complex(dp) :: rho_ss2(2, 2), evals4(4), modes4(2, 2, 4), rho_t(2, 2)
        real(dp) :: gam_d, gam_u, gam_p, dlt, gam_t, t_sp, dt_sp
        real(dp) :: gam_rates_sp(4)
        integer :: stat_sp, istep

        gam_d = 0.6_dp; gam_u = 0.15_dp; gam_p = 0.25_dp; dlt = 0.8_dp
        gam_t = gam_d + gam_u

        call lindblad_steady_state_2level(gam_d, gam_u, gam_p, dlt, rho_ss2, stat_sp)
        call assert_true("Steady state augmented solve = closed form", stat_sp == 0)
        call assert_close("Steady state trace = 1.0", 1.0_dp, &
                          real(rho_ss2(1, 1) + rho_ss2(2, 2), dp), 1.0e-12_dp)

        call lindblad_eigenmodes_2level(gam_d, gam_u, gam_p, dlt, rho_ss2, &
                                        evals4, modes4, stat_sp)
        call assert_true("Damping-basis eigen-relation L R = lambda R", stat_sp == 0)
        call assert_close("lambda_2 = -Gamma", -gam_t, real(evals4(2), dp), 1.0e-12_dp)
        call assert_close("Re lambda_3 = -(Gamma/2+2gp)", &
                          -(0.5_dp * gam_t + 2.0_dp * gam_p), &
                          real(evals4(3), dp), 1.0e-12_dp)
        call assert_close("Im lambda_3 = -delta", -dlt, aimag(evals4(3)), 1.0e-12_dp)

        ! 谱重构 vs RK4: 初态 (|e>+|g>)/sqrt2, H = delta |e><e|
        allocate(rho_spec(2, 2), rho_rk(2, 2), h_bl(2, 2), l_ops_sp(2, 2, 4))
        rho_spec = (0.0_dp, 0.0_dp)
        rho_spec(1, 1) = (0.5_dp, 0.0_dp); rho_spec(2, 2) = (0.5_dp, 0.0_dp)
        rho_spec(1, 2) = (0.5_dp, 0.0_dp); rho_spec(2, 1) = (0.5_dp, 0.0_dp)
        h_bl = (0.0_dp, 0.0_dp)
        h_bl(1, 1) = cmplx(dlt, 0.0_dp, dp)
        l_ops_sp(:, :, 1) = (0.0_dp, 0.0_dp); l_ops_sp(2, 1, 1) = (1.0_dp, 0.0_dp)
        l_ops_sp(:, :, 2) = (0.0_dp, 0.0_dp); l_ops_sp(1, 2, 2) = (1.0_dp, 0.0_dp)
        l_ops_sp(:, :, 3) = (0.0_dp, 0.0_dp); l_ops_sp(1, 1, 3) = (1.0_dp, 0.0_dp)
        l_ops_sp(:, :, 4) = (0.0_dp, 0.0_dp); l_ops_sp(2, 2, 4) = (1.0_dp, 0.0_dp)
        gam_rates_sp(1) = gam_d; gam_rates_sp(2) = gam_u
        gam_rates_sp(3) = 2.0_dp * gam_p; gam_rates_sp(4) = 2.0_dp * gam_p

        t_sp = 0.7_dp
        dt_sp = 1.0e-3_dp
        rho_rk = rho_spec
        do istep = 1, nint(t_sp / dt_sp)
            call rk4_lindblad_step(rho_rk, h_bl, l_ops_sp, gam_rates_sp, dt_sp)
        end do
        call lindblad_spectral_reconstruct(rho_spec, rho_ss2, gam_d, gam_u, &
                                           gam_p, dlt, t_sp, rho_t)
        call assert_close("Spectral reconstruction vs RK4 (t=0.7)", &
                          real(rho_rk(1, 1), dp), real(rho_t(1, 1), dp), 1.0e-8_dp)
        call assert_close("Spectral reconstruction coherence", &
                          real(rho_rk(1, 2), dp), real(rho_t(1, 2), dp), 1.0e-8_dp)
        call assert_close("Reconstruction preserves trace", 1.0_dp, &
                          real(rho_t(1, 1) + rho_t(2, 2), dp), 1.0e-12_dp)
        deallocate(rho_spec, rho_rk, h_bl, l_ops_sp)
    end block

    ! ----------------------------------------------------
    ! 6. 驱动耗散二能级: 饱和稳态 + 量子回归定理 + Mollow 三重峰
    ! ----------------------------------------------------
    block
        complex(dp) :: rho_ss3(2, 2), ev3(4), md3(2, 2, 4), rho_t3(2, 2)
        real(dp), allocatable :: om_g(:), s_g(:), tau_g(:), g1_g(:)
        real(dp) :: om_r, dlt, gam, s_sat, om_eff, s_norm, d_om
        integer :: stat3, i3, n_pk, n_om3

        om_r = 5.0_dp; dlt = 0.0_dp; gam = 1.0_dp
        call optical_bloch_steady_state_driven(om_r, dlt, gam, rho_ss3, stat3)
        call assert_true("Driven steady state = saturation closed form", stat3 == 0)
        s_sat = 2.0_dp * om_r**2 / gam**2
        call assert_close("Driven rho_ee = s/(2(1+s))", &
                          s_sat / (2.0_dp * (1.0_dp + s_sat)), &
                          real(rho_ss3(1, 1), dp), 1.0e-10_dp)

        ! Mollow 谱: 归一化 = rho_ee^ss; 边带 = ±Ω_eff
        n_om3 = 400
        allocate(om_g(n_om3), s_g(n_om3))
        call mollow_spectrum_grid(om_r, dlt, gam, rho_ss3, 12.0_dp, 400, &
                                  n_om3, om_g, s_g, stat3)
        s_norm = sum(s_g) * (om_g(2) - om_g(1))
        call assert_close("Mollow spectrum normalization = rho_ee^ss", &
                          real(rho_ss3(1, 1), dp), s_norm, 2.0e-2_dp)
        om_eff = sqrt(om_r**2 - 0.25_dp * gam**2)
        n_pk = 0
        d_om = om_g(2) - om_g(1)
        do i3 = 2, n_om3 - 1
            if (s_g(i3) > s_g(i3 - 1) .and. s_g(i3) >= s_g(i3 + 1) .and. &
                s_g(i3) > 0.05_dp) then
                if (abs(abs(om_g(i3)) - om_eff) < 3.0_dp * d_om .or. &
                    abs(om_g(i3)) < d_om) n_pk = n_pk + 1
            end if
        end do
        call assert_true("Mollow triplet: central + 2 sidebands at ±Ω_eff", n_pk == 3)

        ! 弱驱动: 单 Lorentz 峰在 delta (无边带), 峰位 = 0
        om_r = 0.2_dp
        call optical_bloch_steady_state_driven(om_r, dlt, gam, rho_ss3, stat3)
        call mollow_spectrum_grid(om_r, dlt, gam, rho_ss3, 12.0_dp, 400, &
                                  n_om3, om_g, s_g, stat3)
        n_pk = 0
        do i3 = 2, n_om3 - 1
            if (s_g(i3) > s_g(i3 - 1) .and. s_g(i3) >= s_g(i3 + 1) .and. &
                s_g(i3) > 0.3_dp * maxval(s_g)) n_pk = n_pk + 1
        end do
        call assert_true("Weak drive: single Lorentzian peak (no sidebands)", n_pk == 1)
        deallocate(om_g, s_g)
    end block

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Open Quantum & OCT Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some tests failed!"
        stop 1
    else
        print '(A)', "SUCCESS: All open quantum & OCT tests passed."
    end if

contains

    subroutine assert_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        n_tests = n_tests + 1
        if (condition) then
            n_passed = n_passed + 1
            print '(A, A42, A)', " [PASS] ", name, ""
        else
            print '(A, A42, A)', " [FAIL] ", name, " (Condition was FALSE)"
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
            print '(A, A42, A)', " [PASS] ", name, ""
        else
            print '(A, A42, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_open_quantum_opt
