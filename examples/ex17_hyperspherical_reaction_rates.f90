! ==============================================================================
! GeneralModule Example 17: Hyperspherical Reactive Scattering & Rate Constants
!
! Features:
!  1. Mass-scaled Delves hyperspherical coordinates and reaction skew angle beta_skew
!  2. Kinematics of isotopic reactions H + H2 -> H2 + H and D + H2 -> HD + H
!  3. Quantum tunneling through Eckart potential barrier and cumulative probability N(E)
!  4. Canonical thermal reaction rate constants k(T) from 200 K to 1200 K
!  5. Comparison between exact quantum rate and Wigner-corrected transition state theory
!
! Standard: Fortran 2008
! ==============================================================================
program ex17_hyperspherical_reaction_rates
    use mod_constants, only: dp, PI, EV2AU, AU2EV, CM2AU
    use mod_hyperspherical_reactive
    implicit none

    type(reaction_mass_t)     :: mass_h_h2, mass_d_h2
    type(transition_state_t)  :: ts_h3
    real(dp) :: e_min, e_max, e_val, prob_t, n_e
    real(dp) :: temp_k, k_exact, k_tst
    integer  :: u_out, i, n_pts

    print *, "================================================================"
    print *, " Example 17: Hyperspherical Reactive Scattering & Rate Constants"
    print *, "================================================================"

    ! 1. Initialize mass kinematics for H + H2 and D + H2
    call init_reaction_mass(1.007825_dp, 1.007825_dp, 1.007825_dp, mass_h_h2)
    call init_reaction_mass(2.014102_dp, 1.007825_dp, 1.007825_dp, mass_d_h2)

    write(*, '(A, F8.3, A)') " H + H2 reaction skew angle     : ", mass_h_h2%skew_angle_deg, " deg (Exact: 60 deg)"
    write(*, '(A, F8.3, A)') " D + H2 reaction skew angle     : ", mass_d_h2%skew_angle_deg, " deg"
    write(*, '(A, F8.4)')    " Delves scaling factor d(H+H2)  : ", mass_h_h2%scale_factor_d
    write(*, '(A, F8.4)')    " Delves scaling factor d(D+H2)  : ", mass_d_h2%scale_factor_d
    print *, "----------------------------------------------------------------"

    ! 2. Setup H3 transition state barrier (V_b = 0.425 eV, omega_imag = 1500 cm^-1)
    ts_h3%v_barrier_au   = 0.425_dp * EV2AU
    ts_h3%omega_im_au    = 1500.0_dp * CM2AU
    ts_h3%omega_bend_au  = 900.0_dp * CM2AU
    ts_h3%omega_symm_au  = 2050.0_dp * CM2AU
    ts_h3%n_trans_states = 5

    ! 3. Scan cumulative reaction probability N(E)
    open(newunit=u_out, file="ex17_hyperspherical_reaction_rates.dat", status="replace", action="write")
    write(u_out, '(A)') "# GeneralModule Example 17: Hyperspherical Reactive Dynamics"
    write(u_out, '(A)') "# Reaction: H + H2 -> H2 + H (Eckart barrier V_b = 0.425 eV)"
    write(u_out, '(A)') "# Part 1: Cumulative Reaction Probability N(E)"
    write(u_out, '(A)') "# Col 1: Total Energy E (eV)"
    write(u_out, '(A)') "# Col 2: Barrier Tunneling Transmission P(E)"
    write(u_out, '(A)') "# Col 3: Cumulative Reaction Probability N(E)"

    n_pts = 50
    e_min = 0.10_dp
    e_max = 1.00_dp
    do i = 1, n_pts
        e_val  = (e_min + real(i - 1, dp) * (e_max - e_min) / real(n_pts - 1, dp)) * EV2AU
        prob_t = calc_eckart_transmission(e_val, ts_h3%v_barrier_au, ts_h3%omega_im_au)
        n_e    = calc_cumulative_reaction_probability(ts_h3, e_val)
        write(u_out, '(3ES16.6)') e_val * AU2EV, prob_t, n_e
    end do

    write(u_out, '(A)') ""
    write(u_out, '(A)') "# Part 2: Thermal Rate Constants k(T)"
    write(u_out, '(A)') "# Col 1: Temperature T (K)"
    write(u_out, '(A)') "# Col 2: Exact Quantum Rate k_exact (cm^3 / s)"
    write(u_out, '(A)') "# Col 3: TST Rate with Wigner Tunneling k_TST (cm^3 / s)"

    write(*, '(A)') " Temperature (K) | k_exact (cm^3/s)     | k_TST (cm^3/s)"
    write(*, '(A)') "-----------------+----------------------+---------------------"

    do i = 1, 10
        temp_k  = 200.0_dp + real(i - 1, dp) * 100.0_dp
        k_exact = calc_canonical_rate_constant(ts_h3, mass_h_h2, temp_k, 500)
        k_tst   = calc_tst_wigner_rate(ts_h3%v_barrier_au, ts_h3%omega_im_au, temp_k, 1.0e-10_dp)
        write(*, '(F16.1, " | ", ES20.6, " | ", ES20.6)') temp_k, k_exact, k_tst
        write(u_out, '(3ES16.6)') temp_k, k_exact, k_tst
    end do

    close(u_out)
    print *, "----------------------------------------------------------------"
    print *, " Output data saved to ex17_hyperspherical_reaction_rates.dat"
    print *, " Example 17 completed successfully."

end program ex17_hyperspherical_reaction_rates
