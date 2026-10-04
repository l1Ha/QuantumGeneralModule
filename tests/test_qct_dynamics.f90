!> \brief 准经典轨迹动力学单元测试 (Quasi-Classical Trajectory Unit Tests)
program test_qct_dynamics
    use mod_constants, only: dp, PI
    use mod_qct_dynamics
    implicit none

    type(qct_leps_param_t) :: par
    type(qct_config_t)     :: cfg
    type(qct_trajectory_t) :: traj
    type(qct_result_t)     :: res

    real(dp) :: q_test(3, 3), force_ana(3, 3), force_num(3, 3)
    real(dp) :: v0, vp, vm, h_step, max_force_diff
    real(dp) :: e_init, e_final, p_init(3), p_final(3)
    integer  :: i_atom, i_dim, n_tests, n_pass, seed
    real(dp) :: sigma, err

    n_tests = 0
    n_pass  = 0
    call init_qct_leps_param(par)

    print '(A)', "=========================================================="
    print '(A)', "       GeneralModule Unit Tests: QCT Dynamics             "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: LEPS 解析受力与中心有限差分梯度一致性检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    h_step = 1.0e-6_dp
    q_test(:, 1) = (/ 0.0_dp,  0.0_dp, 0.0_dp /)
    q_test(:, 2) = (/ 1.4_dp,  0.2_dp, 0.1_dp /)
    q_test(:, 3) = (/ 0.3_dp,  1.5_dp, -0.2_dp /)

    call eval_leps_energy_gradient(q_test, par, v0, force_ana)

    max_force_diff = 0.0_dp
    do i_atom = 1, 3
        do i_dim = 1, 3
            q_test(i_dim, i_atom) = q_test(i_dim, i_atom) + h_step
            call eval_leps_energy_gradient(q_test, par, vp, force_num)

            q_test(i_dim, i_atom) = q_test(i_dim, i_atom) - 2.0_dp * h_step
            call eval_leps_energy_gradient(q_test, par, vm, force_num)

            q_test(i_dim, i_atom) = q_test(i_dim, i_atom) + h_step
            force_num(i_dim, i_atom) = -(vp - vm) / (2.0_dp * h_step)

            max_force_diff = max(max_force_diff, abs(force_ana(i_dim, i_atom) - force_num(i_dim, i_atom)))
        end do
    end do

    if (max_force_diff < 1.0e-5_dp) then
        print '(A, ES10.3)', " [PASS] LEPS analytical forces match numerical gradient: max_diff = ", max_force_diff
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] LEPS analytical forces gradient mismatch: max_diff = ", max_force_diff
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: 单条反应散射轨迹能量守恒与动量守恒检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    seed = 42
    cfg%dt = 2.0_dp
    cfg%max_steps = 1500
    call qct_init_trajectory(cfg, traj, seed, par)
    e_init = traj%e_total
    p_init = sum(traj%p, dim=2)

    do i_atom = 1, 300
        call qct_propagate_step(traj%q, traj%p, cfg%dt, cfg%mass, par)
    end do
    call eval_leps_energy_gradient(traj%q, par, v0, force_ana)
    e_final = sum(traj%p**2) / (2.0_dp * cfg%mass) + v0
    p_final = sum(traj%p, dim=2)

    if (abs(e_final - e_init) < 1.0e-4_dp .and. maxval(abs(p_final - p_init)) < 1.0e-10_dp) then
        print '(A, ES10.3, A, ES10.3)', " [PASS] Trajectory energy conserved: dE = ", abs(e_final - e_init), &
                                        ", dP = ", maxval(abs(p_final - p_init))
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] Trajectory conservation violated: dE = ", abs(e_final - e_init)
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: 终态解离双原子内部能量计算与配对判定
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    ! 人工构造一个已解离的产物构型: 1 与 2 靠近成键 (1.4 bohr), 3 远离 (12.0 bohr)
    traj%q(:, 1) = (/ 0.0_dp,  0.0_dp, 0.0_dp /)
    traj%q(:, 2) = (/ 1.4_dp,  0.0_dp, 0.0_dp /)
    traj%q(:, 3) = (/ 12.0_dp, 0.0_dp, 0.0_dp /)
    traj%p = 0.0_dp
    call qct_analyze_final_state(traj, cfg, par)

    if (traj%reactive .and. traj%e_internal < 0.0_dp) then
        print '(A, F10.4, A)', " [PASS] Product dissociation correctly identified: E_int = ", traj%e_internal, " a.u."
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Product dissociation state analysis failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Monte Carlo 反应截面公式检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call qct_cross_section(n_rxn=25, n_tot=100, b_max=3.0_dp, sigma=sigma, err=err)

    if (abs(sigma - PI * 9.0_dp * 0.25_dp) < 1.0e-12_dp .and. err > 0.0_dp) then
        print '(A, F8.4, A, F8.4)', " [PASS] QCT Monte Carlo cross section formula: sigma = ", sigma, " +/- ", err
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] QCT Monte Carlo cross section formula failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: 小型 QCT 反应系综模拟检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    cfg%n_traj = 30
    cfg%max_steps = 1500
    cfg%seed = 12345
    call run_qct_ensemble(cfg, res, par)

    if (res%n_traj == 30 .and. res%cross_section >= 0.0_dp) then
        print '(A, I3, A, F8.4, A, F8.4)', " [PASS] QCT ensemble run: n_rxn = ", res%n_reactive, &
                                           ", sigma = ", res%cross_section, " +/- ", res%stat_error
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] QCT ensemble run failed"
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  QCT Dynamics Tests: ", n_pass, " / ", n_tests, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_tests) error stop "Test failures detected in test_qct_dynamics."
end program test_qct_dynamics
