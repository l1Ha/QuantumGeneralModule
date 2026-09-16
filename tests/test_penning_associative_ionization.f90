! ==============================================================================
! GeneralModule 单元测试: mod_penning_associative_ionization
! 测试潘宁电离与缔合电离 (Penning & Associative Ionization / Chemi-ionization)
! 验证复光学势转折点、电离概率、分流比、截面与速率、PIES 电子能谱
! ==============================================================================

program test_penning_associative_ionization
    use iso_fortran_env, only: dp => real64
    use mod_penning_associative_ionization
    implicit none

    integer :: num_tests, num_passed
    integer :: stat
    type(penning_system_t) :: sys
    real(dp) :: v_star, v_plus, gamma_inner, gamma_outer
    type(penning_trajectory_result_t) :: traj_res
    type(penning_cross_section_result_t) :: cs_res
    real(dp) :: e_grid(50), pies(50)
    real(dp) :: max_pies, peak_energy
    real(dp) :: a_real, beta_loss, k_loss
    integer  :: i, max_idx

    num_tests = 0
    num_passed = 0

    print *, "=========================================================="
    print *, "  GeneralModule Unit Tests: Penning & Associative Ioniz   "
    print *, "=========================================================="

    ! --------------------------------------------------------------------------
    ! 测试 1: He*(2^3S) + Ar 基准体系初始化与势能参数
    ! --------------------------------------------------------------------------
    num_tests = num_tests + 1
    call init_penning_preset_he_star_ar(sys, stat=stat)
    if (stat == 0 .and. trim(sys%name) == "He*(2^3S) + Ar" .and. &
        abs(sys%excess_energy_ev - 4.06_dp) < 1.0e-2_dp .and. &
        sys%v_star_de_ev > 0.0_dp .and. sys%v_plus_de_ev > 0.0_dp) then
        num_passed = num_passed + 1
        print *, " [PASS] He*(2^3S) + Ar benchmark initialized: E_exc = 19.82 eV, Excess = ", &
                 real(sys%excess_energy_ev, dp), " eV"
    else
        print *, " [FAIL] Benchmark initialization failed."
    end if

    ! --------------------------------------------------------------------------
    ! 测试 2: 相互作用势形态与自电离跃迁宽度指数衰减
    ! --------------------------------------------------------------------------
    num_tests = num_tests + 1
    call eval_penning_potentials(sys, 2.5_dp, v_star, v_plus, gamma_inner)
    call eval_penning_potentials(sys, 6.0_dp, v_star, v_plus, gamma_outer)
    if (gamma_inner > gamma_outer .and. gamma_inner > 0.01_dp .and. &
        gamma_outer < 0.001_dp) then
        num_passed = num_passed + 1
        print *, " [PASS] Autoionization width verified: Gamma(2.5A) = ", &
                 real(gamma_inner, dp), " eV, Gamma(6.0A) = ", real(gamma_outer, dp), " eV"
    else
        print *, " [FAIL] Potential evaluation or width decay failed."
    end if

    ! --------------------------------------------------------------------------
    ! 测试 3: 半经典向径转折点与单碰撞参数电离概率
    ! --------------------------------------------------------------------------
    num_tests = num_tests + 1
    call calc_penning_trajectory_prob(sys, e_coll_ev=0.05_dp, b_ang=2.0_dp, &
                                      res=traj_res, stat=stat)
    if (stat == 0 .and. traj_res%r_turn_ang > 1.5_dp .and. &
        traj_res%p_ion_tot > 0.0_dp .and. traj_res%p_ion_tot <= 1.0_dp .and. &
        traj_res%p_penning > 0.0_dp) then
        num_passed = num_passed + 1
        print *, " [PASS] Classical trajectory: b = 2.0 A, R_turn = ", &
                 real(traj_res%r_turn_ang, dp), " A, P_ion = ", &
                 real(traj_res%p_ion_tot, dp), ", P_AI = ", real(traj_res%p_associative, dp)
    else
        print *, " [FAIL] Trajectory probability calculation failed."
    end if

    ! --------------------------------------------------------------------------
    ! 测试 4: 积分电离截面、缔合电离分支比与反应速率常数
    ! --------------------------------------------------------------------------
    num_tests = num_tests + 1
    call calc_penning_cross_sections(sys, e_coll_ev=0.05_dp, n_b=128, &
                                     b_max_ang=8.0_dp, res=cs_res, stat=stat)
    if (stat == 0 .and. cs_res%sigma_tot_ang2 > 5.0_dp .and. &
        cs_res%sigma_tot_ang2 < 100.0_dp .and. &
        cs_res%ai_branching_ratio >= 0.0_dp .and. &
        cs_res%rate_coeff_cm3_s > 1.0e-12_dp) then
        num_passed = num_passed + 1
        print *, " [PASS] Cross sections verified: sigma_tot = ", &
                 real(cs_res%sigma_tot_ang2, dp), " A^2, AI Branching = ", &
                 real(cs_res%ai_branching_ratio * 100.0_dp, dp), "%, k_tot = ", &
                 real(cs_res%rate_coeff_cm3_s, dp), " cm^3/s"
    else
        print *, " [FAIL] Cross section calculation failed."
    end if

    ! --------------------------------------------------------------------------
    ! 测试 5: PIES 潘宁电子能谱峰值定位与超冷复散射长度
    ! --------------------------------------------------------------------------
    num_tests = num_tests + 1
    call calc_penning_electron_spectrum(sys, e_coll_ev=0.05_dp, n_e=50, &
                                        e_min_ev=3.5_dp, e_max_ev=4.5_dp, &
                                        e_grid_ev=e_grid, pies_intensity=pies, stat=stat)
    max_pies = 0.0_dp
    max_idx = 1
    do i = 1, 50
        if (pies(i) > max_pies) then
            max_pies = pies(i)
            max_idx = i
        end if
    end do
    peak_energy = e_grid(max_idx)

    call calc_ultracold_penning_complex_length(sys, a_real, beta_loss, k_loss, stat=stat)

    if (abs(peak_energy - 4.06_dp) < 0.35_dp .and. beta_loss > 0.0_dp .and. &
        k_loss > 0.0_dp) then
        num_passed = num_passed + 1
        print *, " [PASS] PIES electron peak localized at E_e = ", real(peak_energy, dp), &
                 " eV (E_excess = 4.06 eV), beta_loss = ", real(beta_loss, dp), " A"
    else
        print *, " [FAIL] PIES spectrum or ultracold complex length failed."
    end if

    print *, "=========================================================="
    print '(A, I2, A, I2, A)', "  Test Results: ", num_passed, " / ", num_tests, " passed."
    print *, "=========================================================="

    if (num_passed /= num_tests) then
        error stop "Assertion failure in test_penning_associative_ionization."
    end if

end program test_penning_associative_ionization
