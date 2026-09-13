! ==============================================================================
! GeneralModule: test_hyperspherical_reactive.f90
!
! Unit Test Suite: Reactive Triatomic Quantum Dynamics in Hyperspherical Coordinates
!
! Standard: Fortran 2008
! ==============================================================================

program test_hyperspherical_reactive
    use mod_constants, only: dp, PI, AMU2AU
    use mod_hyperspherical_reactive
    implicit none

    integer :: n_pass, n_total
    type(reaction_mass_t)    :: rmass_h3, rmass_fh2
    type(transition_state_t) :: ts_h3
    real(dp) :: rho, alpha, r_out, r_at_out, p_half, p_low, p_high
    real(dp) :: n_e1, n_e2, k_rate_300k, k_rate_500k, k_tst, kappa_w

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Hyperspherical Reactive"
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: H + H2 Symmetric Kinematics and Skew Angle = 60 degrees
    ! --------------------------------------------------------------------------
    call init_reaction_mass(1.00784_dp * AMU2AU, 1.00784_dp * AMU2AU, 1.00784_dp * AMU2AU, rmass_h3)

    n_total = n_total + 1
    if (abs(rmass_h3%skew_angle_deg - 60.0_dp) < 1.0e-10_dp) then
        print *, " [PASS] H + H2 skew angle beta = 60.000 deg (Exact: 60 deg)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] H + H2 skew angle mismatch: ", rmass_h3%skew_angle_deg
    end if

    n_total = n_total + 1
    if (abs(rmass_h3%scale_factor_d - (4.0_dp / 3.0_dp)**0.25_dp) < 1.0e-10_dp) then
        print *, " [PASS] Delves scaling factor d = (4/3)^(1/4) = ", rmass_h3%scale_factor_d
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Delves scaling factor error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Jacobi <-> Delves Hyperspherical Coordinate Round-Trip
    ! --------------------------------------------------------------------------
    call jacobi_to_hyperspherical(r_diatom=1.40_dp, r_atom_diatom=3.0_dp, &
                                  d_scale=rmass_h3%scale_factor_d, rho=rho, alpha=alpha)
    call hyperspherical_to_jacobi(rho=rho, alpha=alpha, d_scale=rmass_h3%scale_factor_d, &
                                  r_diatom=r_out, r_atom_diatom=r_at_out)

    n_total = n_total + 1
    if (abs(r_out - 1.40_dp) < 1.0e-12_dp .and. abs(r_at_out - 3.0_dp) < 1.0e-12_dp) then
        print *, " [PASS] Jacobi <-> Delves hyperspherical roundtrip exact (< 1e-12)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Hyperspherical roundtrip error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Quantum Tunneling Transmission P(E) through Saddle Barrier
    ! --------------------------------------------------------------------------
    ! At E = V_b, parabolic barrier transmission must strictly equal 1/2
    p_half = calc_eckart_transmission(energy_au=0.0150_dp, v_barrier_au=0.0150_dp, omega_im_au=0.0050_dp)
    n_total = n_total + 1
    if (abs(p_half - 0.5_dp) < 1.0e-12_dp) then
        print *, " [PASS] Quantum transmission at barrier peak P(V_b) = 0.50000 (Exact)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Transmission at barrier peak != 0.5: ", p_half
    end if

    p_low  = calc_eckart_transmission(energy_au=0.0050_dp, v_barrier_au=0.0150_dp, omega_im_au=0.0050_dp)
    p_high = calc_eckart_transmission(energy_au=0.0250_dp, v_barrier_au=0.0150_dp, omega_im_au=0.0050_dp)

    n_total = n_total + 1
    if (p_low < 0.01_dp .and. p_high > 0.99_dp) then
        print *, " [PASS] Tunneling decay (P_low = ", p_low, ") & over-barrier (P_high = ", p_high, ")"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Barrier transmission limits unexpected"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Cumulative Reaction Probability N(E)
    ! --------------------------------------------------------------------------
    ts_h3%v_barrier_au   = 0.0150_dp   ! ~ 9.4 kcal/mol
    ts_h3%omega_im_au    = 0.0070_dp   ! Imaginary reaction frequency
    ts_h3%omega_bend_au  = 0.0040_dp   ! Bending frequency ~ 880 cm^-1
    ts_h3%omega_symm_au  = 0.0090_dp   ! Symmetric stretch ~ 2000 cm^-1
    ts_h3%n_trans_states = 3

    n_e1 = calc_cumulative_reaction_probability(ts_h3, energy_au=0.020_dp)
    n_e2 = calc_cumulative_reaction_probability(ts_h3, energy_au=0.040_dp)

    n_total = n_total + 1
    if (n_e2 > n_e1 .and. n_e1 > 0.0_dp) then
        print *, " [PASS] Cumulative reaction probability monotonic: N(E1) = ", n_e1, " < N(E2) = ", n_e2
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Non-monotonic reaction probability"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Canonical Reaction Rate Constant k(T) & Wigner Tunneling Factor
    ! --------------------------------------------------------------------------
    k_rate_300k = calc_canonical_rate_constant(ts_h3, rmass_h3, temp_kelvin=300.0_dp, n_e_steps=500)
    k_rate_500k = calc_canonical_rate_constant(ts_h3, rmass_h3, temp_kelvin=500.0_dp, n_e_steps=500)

    n_total = n_total + 1
    if (k_rate_500k > k_rate_300k .and. k_rate_300k > 0.0_dp) then
        print *, " [PASS] Thermal reaction rate scales with T: k(300K) = ", k_rate_300k, &
                 " < k(500K) = ", k_rate_500k, " cm^3/s"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Rate constant temperature scaling failure"
    end if

    k_tst = calc_tst_wigner_rate(ts_h3%v_barrier_au, ts_h3%omega_im_au, 300.0_dp, 1.0e-10_dp)
    n_total = n_total + 1
    if (k_tst > 0.0_dp) then
        print *, " [PASS] TST rate with Wigner tunneling evaluated successfully: k_TST = ", k_tst
        n_pass = n_pass + 1
    else
        print *, " [FAIL] TST rate failed"
    end if

    print *, "--------------------------------------------------"
    print *, "Hyperspherical Reactive Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All hyperspherical reactive scattering tests passed."
    else
        stop 1
    end if

end program test_hyperspherical_reactive
