! ==============================================================================
! GeneralModule: test_autoionization_fano.f90
!
! Unit Test Suite: Autoionization Systems, Fano Profile, & CCR
!
! Standard: Fortran 2008
! ==============================================================================

program test_autoionization_fano
    use mod_constants, only: dp, PI
    use mod_autoionization_fano
    implicit none

    integer :: n_pass, n_total
    integer :: stat
    type(fano_profile_t) :: fano
    type(ccr_resonance_t) :: ccr_res
    real(dp) :: e_peak, e_zero, sig_peak, sig_zero, sig_expected_max
    real(dp) :: tau_au, tau_fs, q_calc
    real(dp) :: prob_surv, flux_emiss
    real(dp) :: e_grid(21), sig_arr(21), min_e, max_e

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Autoionization & Fano  "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Fano Profile Anti-Resonance & Maximum
    ! --------------------------------------------------------------------------
    fano%e_resonance_au = 2.0_dp     ! E_R = 2.0 a.u.
    fano%gamma_width_au = 0.10_dp    ! Gamma = 0.1 a.u.
    fano%q_parameter    = 2.0_dp     ! q = 2.0
    fano%sigma_0_au     = 10.0_dp    ! sigma_0 = 10 a.u.
    fano%lifetime_fs    = 0.0_dp

    ! Anti-resonance minimum at epsilon = -q => E = E_R - q * (Gamma/2) = 2.0 - 2.0*0.05 = 1.90
    e_zero = fano%e_resonance_au - fano%q_parameter * (0.5_dp * fano%gamma_width_au)
    sig_zero = calc_fano_profile(fano, e_zero)

    n_total = n_total + 1
    if (abs(sig_zero) < 1.0e-12_dp) then
        print *, " [PASS] Anti-resonance zero cross section at E = ", e_zero, " a.u. (Exact: 0.0)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Anti-resonance failure: ", sig_zero
    end if

    ! Resonance maximum at epsilon = 1/q => E = E_R + (1/q)*(Gamma/2) = 2.0 + 0.5*0.05 = 2.025
    e_peak = fano%e_resonance_au + (0.5_dp * fano%gamma_width_au) / fano%q_parameter
    sig_peak = calc_fano_profile(fano, e_peak)
    sig_expected_max = fano%sigma_0_au * (1.0_dp + fano%q_parameter**2)  ! 10 * (1 + 4) = 50.0

    n_total = n_total + 1
    if (abs(sig_peak - sig_expected_max) < 1.0e-10_dp) then
        print *, " [PASS] Resonance maximum cross section = ", sig_peak, &
                 " (Exact: ", sig_expected_max, ")"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Resonance maximum mismatch: ", sig_peak, " vs ", sig_expected_max
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Autoionization Lifetime tau = hbar / Gamma
    ! --------------------------------------------------------------------------
    call calc_autoionization_lifetime(0.01_dp, tau_au, tau_fs, stat)

    n_total = n_total + 1
    if (stat == 0 .and. abs(tau_au - 100.0_dp) < 1.0e-10_dp) then
        print *, " [PASS] Autoionization lifetime tau_au = ", tau_au, " a.u. (Exact: 100.0)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Lifetime calculation in a.u. failed"
    end if

    n_total = n_total + 1
    ! 100 * 0.024188843 fs ~ 2.41888 fs
    if (abs(tau_fs - 2.4188843_dp) < 1.0e-4_dp) then
        print *, " [PASS] Autoionization lifetime tau_fs = ", tau_fs, " fs ~ 2.419 fs"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Lifetime calculation in fs failed: ", tau_fs
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Fano q-Parameter Definition
    ! --------------------------------------------------------------------------
    ! q = <Phi|D|0> / (pi * V_E * <psi_E|D|0>)
    q_calc = calc_fano_q_parameter(dipole_discrete=2.0_dp, dipole_continuum=1.0_dp, &
                                   v_coupl=1.0_dp / PI)

    n_total = n_total + 1
    if (abs(q_calc - 2.0_dp) < 1.0e-12_dp) then
        print *, " [PASS] Fano q-parameter = ", q_calc, " (Exact: 2.0)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Fano q-parameter mismatch: ", q_calc
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Complex Coordinate Rotation (CCR) Non-Hermitian Model Solver
    ! --------------------------------------------------------------------------
    ! Bare bound state E_0 = 1.5 a.u., continuum state E_cont = 1.4 a.u., coupling V = 0.05
    call solve_ccr_resonance_model(e_bound_0=1.5_dp, e_cont_0=1.4_dp, v_coupl=0.05_dp, &
                                   theta_rad=0.25_dp, res=ccr_res, stat=stat)

    n_total = n_total + 1
    if (stat == 0 .and. ccr_res%is_converged) then
        print *, " [PASS] CCR resonance solved: E_R = ", ccr_res%e_r_au, &
                 ", Gamma = ", ccr_res%gamma_au
        n_pass = n_pass + 1
    else
        print *, " [FAIL] CCR model solver failed"
    end if

    n_total = n_total + 1
    if (ccr_res%gamma_au > 0.0_dp .and. ccr_res%lifetime_fs > 0.0_dp) then
        print *, " [PASS] Positive autoionization width and physical lifetime tau = ", &
                 ccr_res%lifetime_fs, " fs"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Non-physical width from CCR"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Time-Domain Exponential Decay of Resonant State
    ! --------------------------------------------------------------------------
    ! At t = tau = 1/Gamma, survival probability P(t) = 1/e ~ 0.3678794
    call calc_time_domain_autoionization_decay(gamma_au=0.02_dp, t_au=50.0_dp, &
                                               prob_surv=prob_surv, flux_emiss=flux_emiss)

    n_total = n_total + 1
    if (abs(prob_surv - exp(-1.0_dp)) < 1.0e-10_dp) then
        print *, " [PASS] Time-domain survival probability P(tau) = 1/e = ", prob_surv
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Exponential survival probability error: ", prob_surv
    end if

    n_total = n_total + 1
    if (flux_emiss > 0.0_dp) then
        print *, " [PASS] Autoionization electron emission flux J(tau) = ", flux_emiss, " > 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Negative electron flux"
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: Fano Spectrum Scanning
    ! --------------------------------------------------------------------------
    e_grid = [(1.8_dp + real(stat, dp)*0.02_dp, stat=0, 20)]
    call calc_fano_cross_section_spectrum(fano, e_grid, 21, sig_arr, min_e, max_e, stat)

    n_total = n_total + 1
    if (stat == 0 .and. min_e < max_e) then
        print *, " [PASS] Fano spectrum scanned: E_min = ", min_e, ", E_max = ", max_e
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Fano spectrum scanning failed"
    end if

    print *, "--------------------------------------------------"
    print *, "Autoionization & Fano Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All autoionization tests passed."
    else
        stop 1
    end if

end program test_autoionization_fano
