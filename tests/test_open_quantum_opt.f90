!> \brief 开放量子系统 (Lindblad) 与量子最优控制 (Krotov) 单元测试
program test_open_quantum_opt
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    complex(dp) :: rho_pure(2, 2), rho_mixed(2, 2), rho_t(2, 2)
    complex(dp) :: h_2lvl(2, 2), l_relax(2, 2, 1), l_dephase(2, 2, 1)
    real(dp) :: gamma_r(1), gamma_d(1), purity, entropy, coherence
    real(dp) :: dt, t_end
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
