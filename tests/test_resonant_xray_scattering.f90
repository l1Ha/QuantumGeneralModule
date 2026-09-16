!> \file test_resonant_xray_scattering.f90
!> \brief Unit tests for Resonant Inelastic X-ray Scattering (RIXS) and core spectroscopy
!> \author LiHao
program test_resonant_xray_scattering
    use mod_constants, only: dp
    use mod_resonant_xray_scattering
    implicit none

    type(rixs_system_t) :: sys
    real(dp), dimension(1) :: e_core, gamma_core, d_in
    real(dp), dimension(2) :: e_fin, gamma_fin
    real(dp), dimension(2, 1) :: d_out
    real(dp) :: xas_res, xas_off1, xas_off2
    real(dp) :: rixs_res, rixs_off
    real(dp) :: rixs_peak, rixs_left, rixs_right
    real(dp), dimension(5) :: w_in_grid, w_loss_grid
    real(dp), dimension(5, 5) :: rixs_map
    real(dp), dimension(0:4) :: vib_loss
    integer  :: n_pass, n_total, stat

    n_pass = 0
    n_total = 0

    print '(A)', "=========================================================="
    print '(A)', "  GeneralModule Unit Tests: RIXS & Core-Level Spectroscopy "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Setup: Model Cu L3-edge RIXS System (Ground, 2p Core-Hole, dd-Excitation)
    ! --------------------------------------------------------------------------
    e_core = [931.5_dp]        ! Cu 2p3/2 core hole resonance at 931.5 eV
    gamma_core = [0.40_dp]     ! Core hole lifetime HWHM = 0.4 eV
    d_in = [1.0_dp]            ! Incident dipole transition amplitude

    e_fin = [0.0_dp, 1.80_dp]  ! Final states: elastic (0 eV) and dd-loss (1.8 eV)
    gamma_fin = [0.05_dp, 0.05_dp] ! Final state lifetime/resolution HWHM = 0.05 eV
    d_out(1, 1) = 0.80_dp      ! Elastic emission dipole amplitude
    d_out(2, 1) = 0.60_dp      ! Inelastic dd emission dipole amplitude

    call init_rixs_system(sys, e_init=0.0_dp, e_inter=e_core, gamma_core=gamma_core, &
                          d_in=d_in, e_fin=e_fin, gamma_fin=gamma_fin, d_out=d_out, stat=stat)

    ! --------------------------------------------------------------------------
    ! Test 1: System Initialization & XAS Absorption Peak
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    xas_res  = calc_xas_cross_section(sys, 931.5_dp)
    xas_off1 = calc_xas_cross_section(sys, 930.0_dp)
    xas_off2 = calc_xas_cross_section(sys, 933.0_dp)

    if (stat == 0 .and. xas_res > xas_off1 * 5.0_dp .and. xas_res > xas_off2 * 5.0_dp) then
        print '(A, F7.3, A, F6.3, A)', " [PASS] XAS resonance at 931.5 eV verified: peak = ", &
            xas_res, " vs wings ~ ", (xas_off1 + xas_off2) * 0.5_dp
        n_pass = n_pass + 1
    else
        print '(A, 3F10.3)', " [FAIL] XAS peak error: ", xas_res, xas_off1, xas_off2
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Kramers-Heisenberg Resonant Enhancement Factor
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    ! RIXS cross section at dd-excitation loss (1.8 eV) on-resonance vs 15 eV detuned
    rixs_res = calc_kramers_heisenberg_cross_section(sys, 931.5_dp, 1.80_dp)
    rixs_off = calc_kramers_heisenberg_cross_section(sys, 915.0_dp, 1.80_dp)

    if (rixs_res > rixs_off * 100.0_dp .and. rixs_res > 0.0_dp) then
        print '(A, ES11.3, A, ES11.3, A, F8.1, A)', " [PASS] Resonant Kramers-Heisenberg enhancement: ", &
            rixs_res, " / ", rixs_off, " = ", (rixs_res / rixs_off), "x"
        n_pass = n_pass + 1
    else
        print '(A, 2ES12.4)', " [FAIL] Resonant enhancement error: ", rixs_res, rixs_off
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Inelastic Loss Peak Localization
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    rixs_peak  = calc_kramers_heisenberg_cross_section(sys, 931.5_dp, 1.80_dp)
    rixs_left  = calc_kramers_heisenberg_cross_section(sys, 931.5_dp, 1.60_dp)
    rixs_right = calc_kramers_heisenberg_cross_section(sys, 931.5_dp, 2.00_dp)

    if (rixs_peak > rixs_left * 5.0_dp .and. rixs_peak > rixs_right * 5.0_dp) then
        print '(A, F6.2, A, F7.3, A)', " [PASS] Inelastic dd-loss peak localized at Omega = ", &
            1.80_dp, " eV (peak = ", rixs_peak, ")"
        n_pass = n_pass + 1
    else
        print '(A, 3F10.3)', " [FAIL] Inelastic peak localization error: ", &
            rixs_peak, rixs_left, rixs_right
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: 2D RIXS Map Matrix Generation
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    w_in_grid = [930.5_dp, 931.0_dp, 931.5_dp, 932.0_dp, 932.5_dp]
    w_loss_grid = [0.0_dp, 0.9_dp, 1.8_dp, 2.7_dp, 3.6_dp]

    call calc_rixs_2d_map(sys, 5, w_in_grid, 5, w_loss_grid, rixs_map)

    ! Matrix element at (loss=1.8, w_in=931.5) is index (3, 3)
    if (minval(rixs_map) >= 0.0_dp .and. abs(rixs_map(3, 3) - maxval(rixs_map(3, :))) < 1.0e-12_dp) then
        print '(A, F7.3, A)', " [PASS] 2D RIXS intensity map verified: max peak at (w1=931.5, loss=1.8) = ", &
            rixs_map(3, 3), " a.u."
        n_pass = n_pass + 1
    else
        print '(A, 2F10.3)', " [FAIL] 2D RIXS map error: min = ", minval(rixs_map), &
            ", (3,3) = ", rixs_map(3, 3)
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Electron-Phonon Huang-Rhys Vibrational Progression
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call calc_huang_rhys_vibrational_rixs(omega_0_ev=0.060_dp, s_factor=0.50_dp, &
                                         gamma_core_ev=0.35_dp, detuning_ev=0.0_dp, &
                                         n_max_loss=4, loss_intensity=vib_loss)

    ! For S = 0.5, elastic peak (n=0) is dominant, followed by n=1 > n=2 > n=3
    if (vib_loss(0) > vib_loss(1) .and. vib_loss(1) > vib_loss(2) .and. &
        vib_loss(2) > vib_loss(3) .and. minval(vib_loss) > 0.0_dp) then
        print '(A, ES10.3, A, ES10.3, A, ES10.3)', " [PASS] Huang-Rhys phonon progression: I_0 = ", &
            vib_loss(0), ", I_1 = ", vib_loss(1), ", I_2 = ", vib_loss(2)
        n_pass = n_pass + 1
    else
        print '(A, 4ES12.4)', " [FAIL] Huang-Rhys progression error: ", vib_loss(0:3)
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  Test Results: ", n_pass, " / ", n_total, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_total) stop 1

end program test_resonant_xray_scattering
