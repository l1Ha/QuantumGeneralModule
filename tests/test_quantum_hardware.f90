!> \brief 量子计算机上的潘宁电离 (Zhuang et al., JPCL 2026, 17, 9513 算法移植) 验收
!> \details 组 1: 膨胀后选择恒等 (M = exp(-Gamma dt/2) 的态矢范数);
!>          组 2: 时间步收敛 (P_ion 随 dt -> 0 收敛, 二阶分裂算符);
!>          组 3: 与库内经典复势轨迹电离概率同量级对表;
!>          组 4: 资源量级 (qubits / ancilla / steps) 与论文口径一致。
program test_quantum_hardware
    use general_module
    implicit none
    integer :: n_tests = 0, n_passed = 0
    integer :: st

    call print_banner("GeneralModule: Penning Ionization on a Quantum Circuit", 68)

    ! ---------------- 组 1: 膨胀后选择恒等 ----------------
    block
        type(penning_system_t) :: sys
        type(qhw_penning_result_t) :: r1, r2
        real(dp) :: p_ref
        call init_penning_preset_he_star_ar(sys, st)
        ! 单步: P_ion 应由 Gamma(R) 在波包重叠区加权
        call qhw_penning_on_circuit(sys, 0.03_dp, 6, 1, 1.5_dp, 2.0_dp, 9.0_dp, &
                                    6.5_dp, 0.35_dp, .false., r1, st)
        p_ref = 1.0_dp - r1%p_surv
        n_tests = n_tests + 1
        if (st == 0 .and. r1%p_surv > 0.0_dp .and. r1%p_surv <= 1.0_dp .and. &
            abs(p_ref - r1%p_ion) < 1e-14_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2)', " [PASS] dilation postselection identity: p_ion(1 step) =", r1%p_ion
        else
            print '(A, ES10.2, ES10.2, I3)', " [FAIL] dilation:", r1%p_ion, r1%p_surv, st
        end if
    end block

    ! ---------------- 组 2: 时间步收敛 ----------------
    block
        type(penning_system_t) :: sys
        type(qhw_penning_result_t) :: ref, ra, rb
        real(dp) :: ea, eb
        call init_penning_preset_he_star_ar(sys, st)
        ! 穿越时间 ~ 1.5e4 a.u. (v ~ 5.8e-4 a.u., 行程 ~7 bohr)
        call qhw_penning_on_circuit(sys, 0.03_dp, 6, 96000, 2.4e4_dp, 2.0_dp, 9.0_dp, &
                                    6.5_dp, 0.35_dp, .false., ref, st)
        call qhw_penning_on_circuit(sys, 0.03_dp, 6, 6000, 2.4e4_dp, 2.0_dp, 9.0_dp, &
                                    6.5_dp, 0.35_dp, .false., ra, st)
        call qhw_penning_on_circuit(sys, 0.03_dp, 6, 12000, 2.4e4_dp, 2.0_dp, 9.0_dp, &
                                    6.5_dp, 0.35_dp, .false., rb, st)
        ea = abs(ra%p_ion - ref%p_ion)
        eb = abs(rb%p_ion - ref%p_ion)
        n_tests = n_tests + 1
        if (st == 0 .and. eb <= ea .and. eb < 5.0e-2_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] dt convergence (200 -> 400 steps):", &
                  ea, " ->", eb, "  P_ref =", ref%p_ion
        else
            print '(A, 2ES10.2)', " [FAIL] dt convergence:", ea, eb
        end if
    end block

    ! ---------------- 组 3: 与经典轨迹电离概率对表 ----------------
    block
        type(penning_system_t) :: sys
        type(penning_trajectory_result_t) :: tr
        type(qhw_penning_result_t) :: qr
        real(dp) :: p_cls
        call init_penning_preset_he_star_ar(sys, st)
        call calc_penning_trajectory_prob(sys, 0.5_dp, 0.03_dp, tr, st)
        p_cls = tr%p_ion_tot
        call qhw_penning_on_circuit(sys, 0.03_dp, 7, 48000, 2.4e4_dp, 2.0_dp, 9.5_dp, &
                                    7.0_dp, 0.45_dp, .false., qr, st)
        n_tests = n_tests + 1
        if (st == 0 .and. qr%p_ion > 0.05_dp .and. qr%p_ion < 0.999_dp .and. &
            qr%p_ion > 0.25_dp * p_cls .and. qr%p_ion < 4.0_dp * p_cls) then
            n_passed = n_passed + 1
            print '(A, ES9.2, A, ES9.2)', " [PASS] quantum vs classical P_ion:", &
                  qr%p_ion, " vs", p_cls
        else
            print '(A, 2ES10.2, I3)', " [FAIL] quantum vs classical:", qr%p_ion, p_cls, st
        end if
    end block

    ! ---------------- 组 4: 资源量级 ----------------
    block
        type(penning_system_t) :: sys
        type(qhw_penning_result_t) :: qr
        call init_penning_preset_he_star_ar(sys, st)
        call qhw_penning_on_circuit(sys, 0.03_dp, 7, 100, 1.0_dp, 2.0_dp, 9.0_dp, &
                                    6.5_dp, 0.35_dp, .true., qr, st)
        n_tests = n_tests + 1
        if (st == 0 .and. qr%n_grid == 128 .and. qr%n_qubits == 7 .and. &
            qr%n_ancilla == 2 .and. qr%n_steps == 100) then
            n_passed = n_passed + 1
            print '(A, I2, A, I2, A, I2, A, I3)', " [PASS] resources: qubits =", &
                  qr%n_qubits, "  + ancilla", qr%n_ancilla, "  grid", qr%n_grid, &
                  "  steps", qr%n_steps
        else
            print '(A, 3I4)', " [FAIL] resources:", qr%n_grid, qr%n_qubits, qr%n_ancilla
        end if
    end block

    print '(A, I2, A, I2, A)', "Quantum Hardware Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: quantum-hardware Penning test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"
end program test_quantum_hardware
