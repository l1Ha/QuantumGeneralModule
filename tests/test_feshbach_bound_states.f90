! ==============================================================================
! GeneralModule: test_feshbach_bound_states.f90
!
! Unit Test Suite: Magnetic & Optical Feshbach Resonances
!
! Standard: Fortran 2008
! ==============================================================================

program test_feshbach_bound_states
    use mod_constants, only: dp
    use mod_feshbach_bound_states
    implicit none

    integer :: n_pass, n_total, stat
    type(mfr_param_t) :: li6_mfr, rb87_mfr
    type(ofr_param_t) :: ofr
    real(dp) :: a_scat, b_zero, eb_univ, eb_coupled, z_frac, z_far
    real(dp) :: a_re, a_im, k2_peak, k2_detuned

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Feshbach Resonances    "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: 6Li & 87Rb Resonance Classification (s_res >> 1 vs s_res << 1)
    ! --------------------------------------------------------------------------
    call init_mfr_preset("6Li", li6_mfr, stat)
    n_total = n_total + 1
    if (stat == 0 .and. li6_mfr%s_res > 10.0_dp) then
        print *, " [PASS] 6Li broad resonance confirmed: s_res = ", li6_mfr%s_res, " >> 1 (open-channel dominated)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] 6Li s_res mismatch: ", li6_mfr%s_res
    end if

    call init_mfr_preset("87Rb", rb87_mfr, stat)
    n_total = n_total + 1
    if (stat == 0 .and. rb87_mfr%s_res < 0.5_dp) then
        print *, " [PASS] 87Rb narrow resonance: s_res = ", rb87_mfr%s_res, " < 0.5 (closed-channel)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] 87Rb s_res mismatch: ", rb87_mfr%s_res
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Scattering Length Pole & Zero Crossing
    ! --------------------------------------------------------------------------
    ! For 6Li: B0 = 832.18 G, Delta B = -262.3 G -> zero crossing at B = B0 + Delta B = 569.88 G
    b_zero = li6_mfr%b0_gauss + li6_mfr%delta_b_gauss
    a_scat = calc_mfr_scattering_length(li6_mfr, b_zero)

    n_total = n_total + 1
    if (abs(a_scat) < 1.0e-5_dp) then
        print *, " [PASS] Scattering length zero crossing at B = ", b_zero, " G: a(B) = ", a_scat, " a0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Zero crossing mismatch: a = ", a_scat
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Halo Dimer Binding Energy (Universal vs Coupled-Channel)
    ! --------------------------------------------------------------------------
    ! Near 6Li resonance at B = 800 G (a(B) > 0, strongly interacting BEC side)
    a_scat = calc_mfr_scattering_length(li6_mfr, 800.0_dp)
    eb_univ    = calc_mfr_bound_energy_universal(li6_mfr, 800.0_dp)
    eb_coupled = calc_mfr_bound_energy_coupled(li6_mfr, 800.0_dp)

    n_total = n_total + 1
    ! For broad 6Li resonance, R* << a, so eb_coupled should closely match eb_univ
    if (a_scat > 0.0_dp .and. eb_univ > 0.0_dp .and. &
        abs(eb_coupled - eb_univ) / eb_univ < 0.05_dp) then
        print *, " [PASS] 6Li dimer binding energy: E_coupled = ", eb_coupled, " a.u. ~ E_univ = ", eb_univ
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Binding energy mismatch: E_coup = ", eb_coupled, ", E_univ = ", eb_univ
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Closed-Channel Fraction Z(B)
    ! --------------------------------------------------------------------------
    ! Near resonance (B = 830 G), Z(B) << 1 (halo dimer is open-channel dominated)
    z_frac = calc_mfr_closed_channel_fraction(li6_mfr, 830.0_dp)
    ! Far from resonance for narrow 87Rb resonance, Z should approach 1
    z_far  = calc_mfr_closed_channel_fraction(rb87_mfr, 1000.0_dp)

    n_total = n_total + 1
    if (z_frac < 0.05_dp .and. z_far > 0.7_dp) then
        print *, " [PASS] Closed channel fraction: Z(6Li, 830G) = ", z_frac, " << 1, Z(87Rb, 1000G) = ", z_far, " > 0.7"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Closed channel fraction error: Z_6Li = ", z_frac, ", Z_87Rb = ", z_far
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Optical Feshbach Resonance (OFR) Dispersive Tuning & Loss Peak
    ! --------------------------------------------------------------------------
    call init_ofr_param("87Rb", mass_amu=86.91_dp, a_bg_bohr=100.0_dp, &
                        gamma_hz=1.0e7_dp, l_opt_bohr=50.0_dp, ofr=ofr, stat=stat)

    ! On resonance (Delta = 0): delta_a = 0, loss is maximal
    call calc_ofr_complex_scattering_length(ofr, 0.0_dp, a_re, a_im)
    k2_peak = calc_ofr_inelastic_loss_rate(ofr, 0.0_dp)

    n_total = n_total + 1
    if (stat == 0 .and. abs(a_re - ofr%a_bg_au) < 1.0e-10_dp .and. k2_peak > 0.0_dp) then
        print *, " [PASS] OFR on-resonance: Re(a) = a_bg = ", a_re, " a0, peak K2 = ", k2_peak, " cm^3/s"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] OFR resonance check failure"
    end if

    ! Off-resonance (Delta = +gamma / 2): positive scattering length modification delta_a > 0
    call calc_ofr_complex_scattering_length(ofr, 0.5_dp * ofr%gamma_mol_hz, a_re, a_im)
    k2_detuned = calc_ofr_inelastic_loss_rate(ofr, 0.5_dp * ofr%gamma_mol_hz)

    n_total = n_total + 1
    if (a_re > ofr%a_bg_au .and. k2_detuned < k2_peak) then
        print *, " [PASS] OFR detuned: Re(a) = ", a_re, " a0 > a_bg, reduced loss K2 = ", k2_detuned
        n_pass = n_pass + 1
    else
        print *, " [FAIL] OFR detuned behavior unexpected: Re(a) = ", a_re
    end if

    print *, "--------------------------------------------------"
    print *, "Feshbach Resonances Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All Feshbach resonance tests passed."
    else
        stop 1
    end if

end program test_feshbach_bound_states
