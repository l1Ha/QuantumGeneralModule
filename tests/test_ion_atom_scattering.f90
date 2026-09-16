!> \file test_ion_atom_scattering.f90
!> \brief Unit tests for cold ion-atom hybrid scattering and Langevin dynamics
!> \author LiHao
program test_ion_atom_scattering
    use mod_constants, only: dp
    use mod_ion_atom_scattering
    implicit none

    type(ion_atom_system_t) :: sys
    real(dp) :: e_coll_au, b_c1, b_c2, sig1, sig2, k_rate
    real(dp) :: delta0_mere, delta0_num
    real(dp) :: t_limit, heat_rate
    integer  :: n_pass, n_total, stat

    n_pass = 0
    n_total = 0

    print '(A)', "=================================================="
    print '(A)', "  GeneralModule Unit Tests: Ion-Atom Scattering   "
    print '(A)', "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: System Initialization & Characteristic Scales (Yb+ + 6Li)
    ! --------------------------------------------------------------------------
    call init_ion_atom_system(sys, m_ion_amu=171.0_dp, m_atom_amu=6.015_dp, &
                              charge_e=1.0_dp, alpha_au=164.1_dp, stat=stat)
    n_total = n_total + 1
    if (stat == 0 .and. sys%r_star_bohr > 1800.0_dp .and. sys%r_star_bohr < 1900.0_dp .and. &
        sys%e_star_kelvin > 1.0e-7_dp .and. sys%e_star_kelvin < 1.0e-5_dp) then
        print '(A, F8.1, A, ES12.4, A)', " [PASS] Yb+/6Li scale parameters: R* = ", &
            sys%r_star_bohr, " a0, E* = ", sys%e_star_kelvin, " K"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Ion-atom system initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Langevin Critical Impact Parameter & Capture Scaling
    ! --------------------------------------------------------------------------
    e_coll_au = 1.0e-6_dp
    b_c1 = calc_langevin_critical_impact_parameter(sys, e_coll_au)
    b_c2 = calc_langevin_critical_impact_parameter(sys, 16.0_dp * e_coll_au)

    sig1 = calc_langevin_cross_section(sys, e_coll_au)
    sig2 = calc_langevin_cross_section(sys, 16.0_dp * e_coll_au)

    n_total = n_total + 1
    ! E -> 16*E => b_c -> b_c / 2, sigma -> sigma / 4
    if (abs(b_c1 / b_c2 - 2.0_dp) < 1.0e-12_dp .and. abs(sig1 / sig2 - 4.0_dp) < 1.0e-12_dp) then
        print '(A, F10.2, A, F10.2, A)', " [PASS] Langevin scaling verified: b_c(E) = ", &
            b_c1, " a0, b_c(16E) = ", b_c2, " a0"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Langevin scaling violation"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Langevin Energy-Independent Reaction Rate Constant
    ! --------------------------------------------------------------------------
    k_rate = calc_langevin_rate_coefficient(sys)
    n_total = n_total + 1
    ! Typical cold ion-atom rate constant ~ 1e-9 cm^3/s
    if (k_rate > 1.0e-10_dp .and. k_rate < 1.0e-8_dp) then
        print '(A, ES14.4, A)', " [PASS] Langevin reaction rate coefficient: K_L = ", &
            k_rate, " cm^3/s"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Langevin rate coefficient out of physical range"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: MERE & Numerov Ultracold Phase Shift
    ! --------------------------------------------------------------------------
    e_coll_au = 1.0e-8_dp  ! ~3 uK
    call calc_mere_phase_shift_s_wave(sys, e_coll_au, a_s_bohr=150.0_dp, delta0=delta0_mere, stat=stat)
    call calc_ion_atom_phase_shift(sys, e_coll_au, l=0, r_match=5000.0_dp, r_core=20.0_dp, &
                                   phase_shift=delta0_num, stat=stat)

    n_total = n_total + 1
    if (stat == 0 .and. delta0_mere > 0.0_dp .and. delta0_num > 0.0_dp) then
        print '(A, F10.6, A, F10.6, A)', " [PASS] Ultracold s-wave phase shift: delta_MERE = ", &
            delta0_mere, " rad, delta_num = ", delta0_num, " rad"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Quantum phase shift computation failure"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Paul Trap Micromotion-Induced Heating
    ! --------------------------------------------------------------------------
    call calc_rf_micromotion_heating(trap_q=0.25_dp, m_ion=171.0_dp, m_atom=6.015_dp, &
                                    temp_atom_k=1.0e-6_dp, coll_rate_hz=100.0_dp, &
                                    t_limit_k=t_limit, heating_rate_k_s=heat_rate)

    n_total = n_total + 1
    if (t_limit > 1.0e-6_dp .and. heat_rate > 0.0_dp) then
        print '(A, ES12.4, A, ES12.4, A)', " [PASS] RF micromotion heating: T_limit = ", &
            t_limit, " K, dE/dt = ", heat_rate, " K/s"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Micromotion heating evaluation error"
    end if

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', " Ion-Atom Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', " SUCCESS: All cold ion-atom scattering tests passed."
    else
        error stop " Test failures detected in test_ion_atom_scattering."
    end if

end program test_ion_atom_scattering
