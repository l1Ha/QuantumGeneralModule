! ==============================================================================
! GeneralModule: test_bicircular_pecd.f90
!
! Unit Test Suite: Bicircular Laser Fields and Photoelectron Circular Dichroism
!
! Standard: Fortran 2008
! ==============================================================================

program test_bicircular_pecd
    use mod_constants, only: dp, PI
    use mod_bicircular_pecd
    implicit none

    integer :: n_pass, n_total
    type(bicircular_field_t) :: field_counter, field_co
    type(chiral_tetrahedral_molecule_t) :: mol_r, mol_s, mol_achiral
    real(dp) :: chi_r, chi_s, chi_achiral
    real(dp) :: beta1_r, beta1_s, beta1_achiral, g_r, g_s
    integer :: sym_counter, sym_co
    real(dp) :: ex, ey, ax, ay
    real(dp) :: th_grid(15), ph_grid(12), pad(15, 12)

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Bicircular PECD        "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Dynamical C_N Rotational Symmetry Fold
    ! --------------------------------------------------------------------------
    ! Counter-rotating omega + 2*omega (h1 = +1, h2 = -1) -> C3 symmetry
    sym_counter = calc_dynamical_symmetry_fold(h1=1, h2=-1, freq_ratio=2)
    ! Co-rotating omega + 2*omega (h1 = +1, h2 = +1) -> C1 symmetry
    sym_co = calc_dynamical_symmetry_fold(h1=1, h2=1, freq_ratio=2)

    n_total = n_total + 1
    if (sym_counter == 3 .and. sym_co == 1) then
        print *, " [PASS] Dynamical symmetry folds verified: Counter-rotating C", sym_counter, &
                 ", Co-rotating C", sym_co
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Symmetry fold error: counter = ", sym_counter, " co = ", sym_co
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Bicircular Field and Vector Potential Calculation
    ! --------------------------------------------------------------------------
    ! omega1 = 0.057 a.u. (~800 nm), ratio 2, I1 = 1e14 W/cm^2, I2 = 5e13 W/cm^2
    call init_bicircular_field(field_counter, omega1_au=0.057_dp, r_freq=2.0_dp, &
                               i1_wcm2=1.0e14_dp, i2_wcm2=5.0e13_dp, &
                               h1=1, h2=-1, phi1=0.0_dp, phi2=0.0_dp, &
                               fwhm_fs=30.0_dp, envelope_type=1)
    call calc_bicircular_field_at_t(field_counter, t_au=0.0_dp, ex=ex, ey=ey, ax=ax, ay=ay)

    n_total = n_total + 1
    if (abs(ex) > 1.0e-3_dp .and. abs(ay) > 1.0e-3_dp) then
        print *, " [PASS] Bicircular field evaluated at t=0: Ex = ", ex, " a.u., Ay = ", ay, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected bicircular field amplitude at t=0: ", ex, ay
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Chirality Measure Invariance and Mirror Inversion
    ! --------------------------------------------------------------------------
    call init_chiral_tetrahedral_molecule(mol_r, "R")
    call init_chiral_tetrahedral_molecule(mol_s, "S")
    call init_chiral_tetrahedral_molecule(mol_achiral, "ACHIRAL")

    chi_r = calc_chirality_measure(mol_r)
    chi_s = calc_chirality_measure(mol_s)
    chi_achiral = calc_chirality_measure(mol_achiral)

    n_total = n_total + 1
    if (abs(chi_r) > 1.0e-5_dp .and. abs(chi_r + chi_s) < 1.0e-12_dp .and. &
        abs(chi_achiral) < 1.0e-12_dp) then
        print *, " [PASS] Chirality pseudoscalar invariant verified:"
        print *, "        chi(R) = ", chi_r, ", chi(S) = ", chi_s, ", chi(Achiral) = ", chi_achiral
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Chirality measure error: chi(R) = ", chi_r, " chi(S) = ", chi_s
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Enantiomeric PECD Beta1 & Forward-Backward Asymmetry Inversion
    ! --------------------------------------------------------------------------
    beta1_r = calc_chiral_beta1_model(mol_r, energy_ev=5.0_dp, photon_energy_ev=10.0_dp)
    beta1_s = calc_chiral_beta1_model(mol_s, energy_ev=5.0_dp, photon_energy_ev=10.0_dp)
    beta1_achiral = calc_chiral_beta1_model(mol_achiral, energy_ev=5.0_dp, photon_energy_ev=10.0_dp)

    g_r = calc_forward_backward_asymmetry(beta1_r)
    g_s = calc_forward_backward_asymmetry(beta1_s)

    n_total = n_total + 1
    if (abs(beta1_r) > 1.0e-4_dp .and. abs(beta1_r + beta1_s) < 1.0e-12_dp .and. &
        abs(beta1_achiral) < 1.0e-12_dp .and. abs(g_r - 0.5_dp * beta1_r) < 1.0e-12_dp) then
        print *, " [PASS] Enantiomeric PECD asymmetry inversion verified:"
        print *, "        beta1(R) = ", beta1_r, " -> G_PECD(R) = ", g_r
        print *, "        beta1(S) = ", beta1_s, " -> G_PECD(S) = ", g_s
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Enantiomeric beta1 failure: beta1(R)=", beta1_r, " beta1(S)=", beta1_s
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: 2D PAD Spectrum & Forward/Backward Difference
    ! --------------------------------------------------------------------------
    ! Generate PAD for R-enantiomer with beta1 > 0
    call calc_pecd_pad_spectrum(beta1=0.10_dp, beta2=0.20_dp, gamma33=0.01_dp, &
                                n_theta=15, n_phi=12, theta_grid=th_grid, &
                                phi_grid=ph_grid, pad_2d=pad)

    n_total = n_total + 1
    ! Forward is index 1 (theta = 0, P1 = +1), Backward is index 15 (theta = pi, P1 = -1)
    ! Since beta1 > 0, forward emission should be strictly larger than backward emission
    if (pad(1, 1) > pad(15, 1) .and. abs(pad(1, 1) - pad(15, 1)) > 0.01_dp) then
        print *, " [PASS] 2D PAD forward-backward difference verified:"
        print *, "        I(theta=0) = ", pad(1, 1), " > I(theta=pi) = ", pad(15, 1)
        n_pass = n_pass + 1
    else
        print *, " [FAIL] PAD forward/backward unexpected: I(0) = ", pad(1, 1), " I(pi) = ", pad(15, 1)
    end if

    print *, "--------------------------------------------------"
    print *, "Bicircular PECD Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All bicircular PECD tests passed."
    else
        stop 1
    end if

end program test_bicircular_pecd
