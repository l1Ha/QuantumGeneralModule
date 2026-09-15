! ==============================================================================
! GeneralModule: test_grazing_fast_atom_diffraction.f90
!
! Unit Test Suite: Grazing Incidence Fast Atom Diffraction (GIFAD) &
!                  Surface Rainbow Scattering
!
! Standard: Fortran 2008
! ==============================================================================

program test_grazing_fast_atom_diffraction
    use mod_constants, only: dp
    use mod_grazing_fast_atom_diffraction
    implicit none

    integer :: n_pass, n_total
    type(gifad_experiment_t) :: exp_cfg, exp_high_e
    type(gifad_spectrum_t)   :: spec
    real(dp) :: e_perp, lambda_perp, th_rainbow, zeta_reconstructed
    real(dp) :: sum_int, diff_sym
    integer  :: m, idx_pos, idx_neg

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: GIFAD & Rainbow Scatt  "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: GIFAD Initialization & Transverse Fast-Decoupling Kinematics
    ! --------------------------------------------------------------------------
    ! 1.0 keV He beam grazing LiF(001) <110> channel at theta_in = 1.0 degree:
    ! Channel spacing a_x = 2.84 A, corrugation zeta = 0.05 A
    call init_gifad_experiment(exp_cfg, projectile="He", mass_amu=4.0026_dp, &
                              e_kev=1.0_dp, theta_deg=1.0_dp, &
                              ax_ang=2.84_dp, corrugation_ang=0.05_dp)
    call calc_gifad_transverse_kinematics(exp_cfg, e_perp, lambda_perp)

    n_total = n_total + 1
    ! Transverse energy E_perp = 1000 * sin^2(1 deg) ~ 0.3046 eV
    ! Wavelength lambda_perp ~ 0.26 A for He (de Broglie wave regime)
    if (abs(e_perp - 0.3046_dp) < 0.01_dp .and. &
        lambda_perp > 0.20_dp .and. lambda_perp < 0.35_dp) then
        print *, " [PASS] Fast-channeling decoupling verified:"
        print *, "        E_tot = 1.0 keV, E_perp = ", e_perp, " eV"
        print *, "        Transverse de Broglie wavelength lambda_perp = ", lambda_perp, " A"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] GIFAD kinematics error: E_perp = ", e_perp, " lambda = ", lambda_perp
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Classical Rainbow Deflection Angle
    ! --------------------------------------------------------------------------
    th_rainbow = calc_gifad_rainbow_angle(exp_cfg)
    n_total = n_total + 1
    ! slope_max = 2*pi * 0.05 / 2.84 = 0.1106 rad -> theta_R ~ 6.31 deg
    if (th_rainbow > 6.0_dp .and. th_rainbow < 6.6_dp) then
        print *, " [PASS] Classical surface rainbow angle verified: theta_R = ", th_rainbow, " deg"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected rainbow deflection angle: ", th_rainbow
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Inverse Reconstruction of Sub-Picometer Surface Corrugation
    ! --------------------------------------------------------------------------
    zeta_reconstructed = calc_surface_corrugation_from_rainbow(exp_cfg%ax_channel_ang, th_rainbow)
    n_total = n_total + 1
    if (abs(zeta_reconstructed - exp_cfg%corrugation_ang) < 1.0e-5_dp) then
        print *, " [PASS] Sub-picometer corrugation reconstruction accurate:"
        print *, "        zeta_input = ", exp_cfg%corrugation_ang, " A, zeta_recon = ", &
                 zeta_reconstructed, " A"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Corrugation reconstruction mismatch: ", zeta_reconstructed
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Diffraction Spectrum Unitarity & Spatial Inversion Symmetry
    ! --------------------------------------------------------------------------
    call calc_gifad_diffraction_spectrum(exp_cfg, max_order=5, spec=spec)

    sum_int = 0.0_dp
    diff_sym = 0.0_dp
    do m = 1, spec%n_open_orders
        sum_int = sum_int + spec%intensities(m)
    end do

    ! Check inversion symmetry: I(+m) == I(-m) for m in 1..5
    ! For max_order=5, indices run from -5 (1) to +5 (11), center m=0 is at 6.
    do m = 1, 5
        idx_neg = 6 - m
        idx_pos = 6 + m
        diff_sym = max(diff_sym, abs(spec%intensities(idx_pos) - spec%intensities(idx_neg)))
    end do

    n_total = n_total + 1
    if (abs(sum_int - 1.0_dp) < 1.0e-5_dp .and. diff_sym < 1.0e-10_dp) then
        print *, " [PASS] GIFAD diffraction spectrum unitarity & symmetry verified:"
        print *, "        Sum I(m) = ", sum_int, ", Max |I(+m) - I(-m)| = ", diff_sym
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Spectrum unitarity or symmetry failure: sum = ", sum_int, &
                 " diff_sym = ", diff_sym
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Incident Energy & Grazing Angle Scaling
    ! --------------------------------------------------------------------------
    call init_gifad_experiment(exp_high_e, projectile="He", mass_amu=4.0026_dp, &
                              e_kev=3.0_dp, theta_deg=1.5_dp, &
                              ax_ang=2.84_dp, corrugation_ang=0.05_dp)

    n_total = n_total + 1
    ! E_perp must increase with both e_kev and theta_deg:
    ! E_perp(3 keV, 1.5 deg) = 3000 * sin^2(1.5 deg) ~ 2.05 eV > 0.30 eV
    ! lambda_perp must decrease correspondingly
    if (exp_high_e%e_perp_ev > exp_cfg%e_perp_ev .and. &
        exp_high_e%lambda_perp_ang < exp_cfg%lambda_perp_ang) then
        print *, " [PASS] Fast-atom grazing energy scaling verified:"
        print *, "        E_perp(3keV, 1.5deg) = ", exp_high_e%e_perp_ev, &
                 " eV > E_perp(1keV, 1.0deg) = ", exp_cfg%e_perp_ev, " eV"
        print *, "        lambda_perp(high) = ", exp_high_e%lambda_perp_ang, &
                 " A < lambda_perp(low) = ", exp_cfg%lambda_perp_ang, " A"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Energy scaling error"
    end if

    print *, "--------------------------------------------------"
    print *, "GIFAD Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All GIFAD tests passed."
    else
        stop 1
    end if

end program test_grazing_fast_atom_diffraction
