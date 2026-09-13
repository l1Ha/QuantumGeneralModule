! ==============================================================================
! GeneralModule: test_three_body_recombination.f90
!
! Unit Test Suite: Ultracold Three-Body Recombination & Efimov Physics
!
! Standard: Fortran 2008
! ==============================================================================

program test_three_body_recombination
    use mod_constants, only: dp, PI, AMU2AU
    use mod_three_body_recombination
    implicit none

    integer :: n_pass, n_total
    integer :: stat
    real(dp) :: s0_boson, lambda_scale, s0_hetero1, s0_hetero2
    type(efimov_param_t) :: param
    type(three_body_loss_t) :: res_pos, res_neg, res_zero, res_min, res_peak
    real(dp) :: mass_rb87, a_pos, a_neg, a_min, a_peak
    real(dp) :: k3_t1_au, k3_t1_si, k3_t2_au, k3_t2_si
    real(dp) :: a_grid(5), k3_au_arr(5), k3_si_arr(5)

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, "  GeneralModule Unit Tests: Three-Body & Efimov   "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Transcendental Efimov Root s_0 for Identical Bosons
    ! --------------------------------------------------------------------------
    s0_boson = solve_efimov_s0_identical_bosons()
    n_total = n_total + 1
    if (abs(s0_boson - 1.0062378_dp) < 1.0e-5_dp) then
        print *, " [PASS] Identical boson s_0 = ", s0_boson, " (Exact: 1.00624)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] s_0 unexpected: ", s0_boson
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Discrete Scaling Factor lambda = exp(pi / s_0) ~ 22.694
    ! --------------------------------------------------------------------------
    lambda_scale = calc_efimov_scale_factor(s0_boson)
    n_total = n_total + 1
    if (abs(lambda_scale - 22.694_dp) < 0.05_dp) then
        print *, " [PASS] Discrete scale factor lambda = ", lambda_scale, " ~ 22.694"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] lambda scale factor error: ", lambda_scale
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Heteronuclear Root with Mass Ratios (AAB systems)
    ! --------------------------------------------------------------------------
    ! K + Rb + Rb (beta = 40 / 87 ~ 0.46)
    s0_hetero1 = solve_efimov_s0_heteronuclear(40.0_dp / 87.0_dp)
    ! Li + Rb + Rb (beta = 6 / 87 ~ 0.069)
    s0_hetero2 = solve_efimov_s0_heteronuclear(6.0_dp / 87.0_dp)

    n_total = n_total + 1
    if (s0_hetero1 > 0.0_dp .and. s0_hetero2 > 0.0_dp) then
        print *, " [PASS] Heteronuclear root solved: s0(40/87)=", s0_hetero1, &
                 ", s0(6/87)=", s0_hetero2
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Heteronuclear root failed"
    end if

    n_total = n_total + 1
    if (s0_hetero2 > s0_hetero1) then
        print *, " [PASS] Light-heavy mass ratio enhances Efimov attraction: ", &
                 s0_hetero2, " > ", s0_hetero1
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Heteronuclear mass ratio scaling violation"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Positive Scattering Length K_3(a > 0)
    ! --------------------------------------------------------------------------
    mass_rb87 = 87.0_dp * AMU2AU
    param%s0 = s0_boson
    param%scale_factor = lambda_scale
    param%a_star = 200.0_dp   ! a_* in a_0
    param%a_minus = 100.0_dp  ! a_- in a_0
    param%eta_star = 0.06_dp  ! Inelasticity parameter

    a_pos = 500.0_dp  ! a = 500 a_0
    call calc_three_body_recombination_a_positive(a_pos, mass_rb87, param, res_pos, stat)

    n_total = n_total + 1
    if (stat == 0) then
        print *, " [PASS] K_3(a > 0) status = 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] K_3(a > 0) calculation error"
    end if

    n_total = n_total + 1
    if (res_pos%k3_au > 0.0_dp .and. res_pos%k3_si > 0.0_dp) then
        print *, " [PASS] K_3(a=500 a0) = ", res_pos%k3_si, " cm^6/s > 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] K_3(a > 0) non-positive"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Efimov Interference Minimum (a > 0)
    ! --------------------------------------------------------------------------
    ! When s_0 * ln(a / a_*) = pi => a = a_* * exp(pi / s_0) = a_* * lambda
    a_min = param%a_star * exp(PI / param%s0)
    call calc_three_body_recombination_a_positive(a_min, mass_rb87, param, res_min, stat)

    n_total = n_total + 1
    if (res_min%is_minimum) then
        print *, " [PASS] Efimov interference minimum detected at a = ", a_min
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Interference minimum detection failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: Negative Scattering Length K_3(a < 0) & Resonance Peak
    ! --------------------------------------------------------------------------
    a_neg = -140.0_dp
    call calc_three_body_recombination_a_negative(a_neg, mass_rb87, param, res_neg, stat)

    n_total = n_total + 1
    if (stat == 0 .and. res_neg%k3_si > 0.0_dp) then
        print *, " [PASS] K_3(a < 0) = ", res_neg%k3_si, " cm^6/s > 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] K_3(a < 0) calculation error"
    end if

    ! Resonance pole at a = -a_minus
    a_peak = - param%a_minus
    call calc_three_body_recombination_a_negative(a_peak, mass_rb87, param, res_peak, stat)

    n_total = n_total + 1
    if (res_peak%is_resonance) then
        print *, " [PASS] Efimov resonance pole detected at a = ", a_peak
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Resonance pole detection failed"
    end if

    n_total = n_total + 1
    if (res_peak%k3_si > 3.0_dp * res_neg%k3_si) then
        print *, " [PASS] Efimov resonance enhancement: K_3,peak / K_3,bg = ", &
                 res_peak%k3_si / res_neg%k3_si, " > 3"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Resonance enhancement insufficient"
    end if

    ! --------------------------------------------------------------------------
    ! Test 7: Universal Dispatcher & Zero Scattering Length Protection
    ! --------------------------------------------------------------------------
    call calc_three_body_recombination_universal(0.0_dp, mass_rb87, param, res_zero, stat)
    n_total = n_total + 1
    if (stat == 0 .and. res_zero%k3_au == 0.0_dp) then
        print *, " [PASS] Zero scattering length safely returns zero loss"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Zero scattering length protection failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 8: Unitary Limit Finite-Temperature Scaling K_3(T) ~ T^(-2)
    ! --------------------------------------------------------------------------
    call calc_unitary_three_body_loss_temperature(1.0e-6_dp, mass_rb87, param%eta_star, &
                                                 k3_t1_au, k3_t1_si, stat)
    call calc_unitary_three_body_loss_temperature(2.0e-6_dp, mass_rb87, param%eta_star, &
                                                 k3_t2_au, k3_t2_si, stat)

    n_total = n_total + 1
    if (stat == 0 .and. k3_t1_si > 0.0_dp) then
        print *, " [PASS] Unitary finite-T loss K_3(1 uK) = ", k3_t1_si, " cm^6/s"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unitary finite-T loss calculation error"
    end if

    n_total = n_total + 1
    ! Ratio should be (2.0)^2 = 4.0
    if (abs(k3_t1_si / k3_t2_si - 4.0_dp) < 1.0e-4_dp) then
        print *, " [PASS] Exact T^(-2) power-law scaling confirmed: K_3(1uK)/K_3(2uK) = ", &
                 k3_t1_si / k3_t2_si, " ~ 4.0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] T^(-2) power-law scaling violation: ", k3_t1_si / k3_t2_si
    end if

    ! --------------------------------------------------------------------------
    ! Test 9: Spectrum Scanning Routine
    ! --------------------------------------------------------------------------
    a_grid = [-500.0_dp, -200.0_dp, 0.0_dp, 200.0_dp, 500.0_dp]
    call scan_efimov_recombination_spectrum(a_grid, 5, mass_rb87, param, &
                                           k3_au_arr, k3_si_arr, stat)

    n_total = n_total + 1
    if (stat == 0 .and. k3_si_arr(1) > 0.0_dp .and. k3_si_arr(5) > 0.0_dp) then
        print *, " [PASS] Efimov spectrum scan executed successfully across 5 points"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Spectrum scan failed"
    end if

    print *, "--------------------------------------------------"
    print *, "Three-Body & Efimov Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All three-body recombination tests passed."
    else
        stop 1
    end if

end program test_three_body_recombination
