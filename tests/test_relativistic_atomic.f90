!> \file test_relativistic_atomic.f90
!> \brief Unit tests for relativistic atomic structure and radial Dirac solver
!> \author LiHao
program test_relativistic_atomic
    use mod_constants, only: dp
    use mod_relativistic_atomic
    implicit none

    type(dirac_state_t) :: state_1s, state_2p12, state_2p32, state_cs_bare, state_cs_pol
    real(dp) :: delta_fs_ev, delta_fs_cm1, osc_str, osc_zero
    real(dp) :: v_near, v_far
    integer  :: n_pass, n_total, stat

    n_pass = 0
    n_total = 0

    print '(A)', "=========================================================="
    print '(A)', "  GeneralModule Unit Tests: Relativistic Atomic & Dirac   "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Dirac Hydrogen 1s Ground State
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call solve_radial_dirac_eigenvalue(z_nuclear=1.0_dp, z_ion=1.0_dp, &
                                       alpha_core=0.0_dp, r_cut=0.0_dp, &
                                       n_princ=1, kappa=-1, state=state_1s, stat=stat)

    if (stat == 0 .and. abs(state_1s%energy_ev - (-13.606_dp)) < 0.05_dp .and. &
        state_1s%l_orb == 0 .and. state_1s%two_j == 1) then
        print '(A, F9.4, A)', " [PASS] Dirac Hydrogen 1s1/2 ground state energy: ", &
            state_1s%energy_ev, " eV (Sommerfeld exact)"
        n_pass = n_pass + 1
    else
        print '(A, F9.4)', " [FAIL] Dirac 1s1/2 energy error: ", state_1s%energy_ev
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Dirac Fine-Structure Splitting (2p1/2 vs 2p3/2 in Hydrogen)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call solve_radial_dirac_eigenvalue(z_nuclear=1.0_dp, z_ion=1.0_dp, &
                                       alpha_core=0.0_dp, r_cut=0.0_dp, &
                                       n_princ=2, kappa=1, state=state_2p12, stat=stat)
    call solve_radial_dirac_eigenvalue(z_nuclear=1.0_dp, z_ion=1.0_dp, &
                                       alpha_core=0.0_dp, r_cut=0.0_dp, &
                                       n_princ=2, kappa=-2, state=state_2p32, stat=stat)

    call calc_dirac_fine_structure_splitting(state_2p12, state_2p32, delta_fs_ev, delta_fs_cm1)

    ! Hydrogen 2p fine-structure splitting is ~0.365 cm^-1 (4.53e-5 eV)
    if (state_2p32%energy_au > state_2p12%energy_au .and. &
        delta_fs_cm1 > 0.30_dp .and. delta_fs_cm1 < 0.40_dp) then
        print '(A, ES12.4, A, F7.4, A)', " [PASS] Fine structure splitting 2p3/2 - 2p1/2: ", &
            delta_fs_ev, " eV (", delta_fs_cm1, " cm^-1)"
        n_pass = n_pass + 1
    else
        print '(A, F10.5)', " [FAIL] Fine structure splitting error: ", delta_fs_cm1
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Model Potential Asymptotics (Near-core -Z/r vs Asymptotic -1/r)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    v_near = calc_dirac_model_potential(z_nuclear=55.0_dp, z_ion=1.0_dp, &
                                        alpha_core=19.0_dp, r_cut=2.0_dp, r=0.01_dp)
    v_far  = calc_dirac_model_potential(z_nuclear=55.0_dp, z_ion=1.0_dp, &
                                        alpha_core=19.0_dp, r_cut=2.0_dp, r=50.0_dp)

    ! Near: ~ -55 / 0.01 = -5500; Far: ~ -1 / 50 = -0.02
    if (v_near < -4000.0_dp .and. abs(v_far - (-0.02_dp)) < 0.005_dp) then
        print '(A, F9.1, A, F8.4, A)', " [PASS] Model potential asymptotics: V(0.01) = ", &
            v_near, " a.u., V(50.0) = ", v_far, " a.u."
        n_pass = n_pass + 1
    else
        print '(A, 2F12.4)', " [FAIL] Model potential error: ", v_near, v_far
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Core Polarization Effect and Quantum Defect
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call solve_radial_dirac_eigenvalue(z_nuclear=55.0_dp, z_ion=1.0_dp, &
                                       alpha_core=0.0_dp, r_cut=2.0_dp, &
                                       n_princ=6, kappa=-1, state=state_cs_bare)
    call solve_radial_dirac_eigenvalue(z_nuclear=55.0_dp, z_ion=1.0_dp, &
                                       alpha_core=19.0_dp, r_cut=2.0_dp, &
                                       n_princ=6, kappa=-1, state=state_cs_pol)

    ! Polarization stabilizes state -> lower energy -> positive quantum defect
    if (state_cs_pol%energy_au < state_cs_bare%energy_au .and. &
        state_cs_pol%quantum_defect > state_cs_bare%quantum_defect) then
        print '(A, F7.4, A, F7.4)', " [PASS] Core polarization shifts: mu(bare) = ", &
            state_cs_bare%quantum_defect, " -> mu(pol) = ", state_cs_pol%quantum_defect
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Core polarization shift error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Relativistic E1 Transition and Oscillator Strength
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call calc_dirac_e1_matrix_element(state_1s, state_2p32, r_overlap_au=1.20_dp, &
                                      osc_strength=osc_str)
    ! Reverse transition should yield zero absorption oscillator strength
    call calc_dirac_e1_matrix_element(state_2p32, state_1s, r_overlap_au=1.20_dp, &
                                      osc_strength=osc_zero)

    if (osc_str > 0.0_dp .and. abs(osc_zero) < 1.0e-12_dp) then
        print '(A, F7.4, A)', " [PASS] Relativistic E1 oscillator strength: f(1s1/2 -> 2p3/2) = ", &
            osc_str, " (emission f = 0.0)"
        n_pass = n_pass + 1
    else
        print '(A, 2F10.4)', " [FAIL] E1 oscillator strength error: ", osc_str, osc_zero
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  Test Results: ", n_pass, " / ", n_total, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_total) stop 1

end program test_relativistic_atomic
