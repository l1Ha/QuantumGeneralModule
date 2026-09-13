! ==============================================================================
! GeneralModule: test_strong_field_nsdi.f90
!
! Unit Test Suite: Strong-Field Non-Sequential Double Ionization & Recollision
!
! Standard: Fortran 2008
! ==============================================================================

program test_strong_field_nsdi
    use mod_constants, only: dp
    use mod_strong_field_nsdi
    implicit none

    integer :: n_pass, n_total, stat
    type(nsdi_laser_t)  :: laser
    type(nsdi_target_t) :: he_target, ar_target
    real(dp) :: gamma_k, w1, w2, w_zero
    real(dp) :: phi_0, phi_r, e_rec, ratio_up
    real(dp) :: sigma_sub, sigma_above
    real(dp) :: pz1, pz2, corr_c
    real(dp), allocatable :: grid_p(:), dist_2d(:,:)
    real(dp) :: intensities(5), y_nsdi(5), y_sdi(5)
    logical  :: ok

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Strong-Field NSDI      "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Laser Initialization & Ponderomotive Energy Up
    ! --------------------------------------------------------------------------
    ! 800 nm, 2.0e14 W/cm^2: Up ~ 0.44 a.u. (~ 12 eV), F0 ~ 0.075 a.u.
    call init_nsdi_laser(800.0_dp, 2.0e14_dp, laser, stat)
    n_total = n_total + 1
    if (stat == 0 .and. abs(laser%up_au - 0.439_dp) < 0.05_dp .and. laser%field_peak_au > 0.07_dp) then
        print *, " [PASS] Laser params: F0 = ", laser%field_peak_au, " a.u., Up = ", laser%up_au, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected laser parameters: Up = ", laser%up_au
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Target Atom & Keldysh Parameter
    ! --------------------------------------------------------------------------
    call init_nsdi_target("He", he_target, stat)
    gamma_k = calc_keldysh_gamma(laser, he_target%ip1_au)
    n_total = n_total + 1
    if (stat == 0 .and. gamma_k > 0.8_dp .and. gamma_k < 1.3_dp) then
        print *, " [PASS] Helium target initialized: Ip1 = ", he_target%ip1_au, " a.u., gamma = ", gamma_k
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Keldysh parameter error: gamma = ", gamma_k
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: ADK Tunneling Rate Monotonicity & Vanishing at Zero Field
    ! --------------------------------------------------------------------------
    w1 = calc_adk_rate(laser%field_peak_au, he_target%ip1_au, 1.0_dp)
    w2 = calc_adk_rate(laser%field_peak_au, he_target%ip2_au, 2.0_dp)
    w_zero = calc_adk_rate(0.0_dp, he_target%ip1_au, 1.0_dp)

    n_total = n_total + 1
    if (w1 > 0.0_dp .and. w1 > 1.0e4_dp * w2 .and. w_zero < 1.0e-30_dp) then
        print *, " [PASS] ADK tunneling: w(Ip1) = ", w1, " >> w(Ip2) = ", w2, ", w(0) = 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] ADK tunneling rate failure: w1 = ", w1, ", w2 = ", w2
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Corkum 3.17 Up Recollision Cutoff
    ! --------------------------------------------------------------------------
    ! Maximum return energy occurs at birth phase phi_0 ~ 0.297 rad (~ 17 deg)
    phi_0 = 0.2972_dp
    call calc_recollision_trajectory(phi_0, laser%up_au, phi_r, e_rec, ok)
    ratio_up = e_rec / laser%up_au

    n_total = n_total + 1
    if (ok .and. abs(ratio_up - 3.173_dp) < 0.05_dp .and. phi_r > 4.3_dp .and. phi_r < 4.6_dp) then
        print *, " [PASS] Corkum recollision cutoff: E_rec = ", ratio_up, " Up (Theory: 3.173 Up), phi_r = ", phi_r
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Corkum cutoff mismatch: ratio = ", ratio_up, ", phi_r = ", phi_r
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Lotz (e,2e) Impact Ionization Threshold
    ! --------------------------------------------------------------------------
    sigma_sub   = calc_lotz_cross_section(1.5_dp, 2.0_dp)
    sigma_above = calc_lotz_cross_section(3.0_dp, 2.0_dp)

    n_total = n_total + 1
    if (sigma_sub < 1.0e-30_dp .and. sigma_above > 0.0_dp) then
        print *, " [PASS] Lotz threshold: sigma(E < Ip2) = 0, sigma(3.0 a.u. > 2.0 a.u.) = ", sigma_above
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Lotz cross section behavior error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: Correlated Longitudinal Drift Momenta
    ! --------------------------------------------------------------------------
    ! For phi_r ~ 4.45 rad (sin(phi_r) < 0), both electrons acquire same-sign drift
    call calc_nsdi_drift_momenta(phi_r, 3.173_dp * laser%up_au, 1.0_dp, &
                                 laser%field_peak_au, laser%omega_au, &
                                 0.5_dp, pz1, pz2)

    n_total = n_total + 1
    if (pz1 * pz2 > 0.0_dp) then
        print *, " [PASS] Drift momenta correlated in same hemisphere: pz1 = ", pz1, ", pz2 = ", pz2
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Momenta not correlated in same hemisphere: pz1 = ", pz1, ", pz2 = ", pz2
    end if

    ! --------------------------------------------------------------------------
    ! Test 7: 2D Correlated Momentum Spectrum P(pz1, pz2) & Positive Correlation
    ! --------------------------------------------------------------------------
    call init_nsdi_target("Ar", ar_target, stat)
    allocate(grid_p(21), dist_2d(21, 21))
    call calc_nsdi_2d_momentum_dist(laser, ar_target, 21, 2.5_dp, grid_p, dist_2d, corr_c)

    n_total = n_total + 1
    ! In direct (e,2e) NSDI, corr_c is strongly positive (Weber et al. Nature 2000)
    if (corr_c > 0.3_dp) then
        print *, " [PASS] 2D momentum correlation coefficient C_corr = ", corr_c, " > 0.3 (COLTRIMS signature)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Insufficient momentum correlation: C_corr = ", corr_c
    end if
    deallocate(grid_p, dist_2d)

    ! --------------------------------------------------------------------------
    ! Test 8: Non-Sequential "Knee Structure" (NSDI >> SDI at Moderate Intensity)
    ! --------------------------------------------------------------------------
    call calc_double_ion_yield_curve(800.0_dp, ar_target, 5, 1.0e14_dp, 5.0e14_dp, &
                                     intensities, y_nsdi, y_sdi)

    n_total = n_total + 1
    ! At 1.0e14 W/cm^2, NSDI should exceed sequential SDI by multiple orders of magnitude
    if (y_nsdi(1) > 1.0e2_dp * y_sdi(1)) then
        print *, " [PASS] Knee structure verified: At 1e14 W/cm2, Y_NSDI / Y_SDI = ", y_nsdi(1) / y_sdi(1), " >> 1"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Knee structure missing: Y_NSDI = ", y_nsdi(1), ", Y_SDI = ", y_sdi(1)
    end if

    print *, "--------------------------------------------------"
    print *, "Strong-Field NSDI Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All strong-field NSDI tests passed."
    else
        stop 1
    end if

end program test_strong_field_nsdi
