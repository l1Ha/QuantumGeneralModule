! ==============================================================================
! GeneralModule Example 28: Grazing Incidence Fast Atom Diffraction (GIFAD)
!
! Scientific Engineering Application:
!  1. keV He beam grazing LiF(001) along <110> channel
!  2. Fast axial channeling decoupling: E_tot = 1.0 keV, E_perp = 0.305 eV
!  3. Transverse quantum diffraction spectrum with classical surface rainbow envelope
!  4. Inverse sub-picometer surface corrugation reconstruction zeta from theta_R
!  5. Systematic scan of rainbow angles vs grazing angles and incident beam energies
!
! Reference:
!   Rousseau et al., Phys. Rev. Lett. 98, 016104 (2007);
!   Schuller et al., Phys. Rev. Lett. 99, 136106 (2007).
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================
program ex28_grazing_fast_atom_diffraction
    use mod_constants, only: dp
    use mod_grazing_fast_atom_diffraction
    implicit none

    type(gifad_experiment_t) :: exp_cfg, scan_cfg
    type(gifad_spectrum_t)   :: spec
    real(dp) :: e_perp, lambda_perp, th_rainbow, zeta_recon
    integer  :: u_spec, u_scan, m, step
    real(dp) :: th_in, e_kev

    print *, "================================================================"
    print *, " Example 28: Grazing Incidence Fast Atom Diffraction (GIFAD)   "
    print *, "================================================================"

    ! 1. Initialize GIFAD experiment: 1.0 keV He grazing LiF(001) <110>
    ! Channel spacing a_x = 2.84 A, corrugation amplitude zeta = 0.05 A
    call init_gifad_experiment(exp_cfg, projectile="He", mass_amu=4.0026_dp, &
                              e_kev=1.0_dp, theta_deg=1.0_dp, &
                              ax_ang=2.84_dp, corrugation_ang=0.05_dp)

    call calc_gifad_transverse_kinematics(exp_cfg, e_perp, lambda_perp)
    th_rainbow = calc_gifad_rainbow_angle(exp_cfg)
    zeta_recon = calc_surface_corrugation_from_rainbow(exp_cfg%ax_channel_ang, th_rainbow)

    write(*, '(A, A)')        " Projectile Species          : ", trim(exp_cfg%projectile_name)
    write(*, '(A, F8.2, A)')  " Total Beam Energy E_tot     : ", exp_cfg%e_beam_kev, " keV"
    write(*, '(A, F8.2, A)')  " Grazing Angle theta_in      : ", exp_cfg%theta_in_deg, " deg"
    write(*, '(A, F8.4, A)')  " Transverse Energy E_perp    : ", e_perp, " eV"
    write(*, '(A, F8.4, A)')  " Transverse de Broglie Wave  : ", lambda_perp, " A"
    write(*, '(A, F8.3, A)')  " Classical Rainbow Angle     : ", th_rainbow, " deg"
    write(*, '(A, F8.4, A)')  " Input Corrugation zeta      : ", exp_cfg%corrugation_ang, " A"
    write(*, '(A, F8.4, A)')  " Reconstructed Corrugation   : ", zeta_recon, " A"
    print *, "----------------------------------------------------------------"

    ! 2. Compute 1D Transverse Diffraction Spectrum
    call calc_gifad_diffraction_spectrum(exp_cfg, max_order=15, spec=spec)

    open(newunit=u_spec, file="ex28_gifad_diffraction_spectrum.dat", &
         status="replace", action="write")
    write(u_spec, '(A)') "# GeneralModule Example 28: GIFAD Transverse Diffraction Spectrum"
    write(u_spec, '(A)') "# System: He (1 keV) / LiF(001) <110>, theta_in = 1.0 deg"
    write(u_spec, '(A)') "# Col 1: Diffraction order m"
    write(u_spec, '(A)') "# Col 2: Deflection angle theta_m (deg)"
    write(u_spec, '(A)') "# Col 3: Normalized diffraction intensity I(m)"

    do m = 1, spec%n_open_orders
        write(u_spec, '(I4, 2ES16.6)') spec%orders(m), spec%angles_deg(m), spec%intensities(m)
    end do
    close(u_spec)
    print *, " Saved diffraction spectrum to ex28_gifad_diffraction_spectrum.dat"

    ! 3. Parametric Scan: Rainbow Angle vs Grazing Angle & Beam Energy
    open(newunit=u_scan, file="ex28_gifad_rainbow_scan.dat", &
         status="replace", action="write")
    write(u_scan, '(A)') "# GeneralModule Example 28: GIFAD Kinematic Scaling Scan"
    write(u_scan, '(A)') "# Col 1: Grazing angle theta_in (deg)"
    write(u_scan, '(A)') "# Col 2: Transverse energy E_perp (eV)"
    write(u_scan, '(A)') "# Col 3: Transverse wavelength lambda_perp (A)"
    write(u_scan, '(A)') "# Col 4: Beam energy E_tot (keV)"
    write(u_scan, '(A)') "# Col 5: Classical rainbow angle theta_R (deg)"

    do step = 1, 50
        th_in = 0.2_dp + real(step, dp) * 0.05_dp     ! 0.25 to 2.70 deg
        e_kev = 0.5_dp + real(step, dp) * 0.10_dp     ! 0.60 to 5.50 keV

        call init_gifad_experiment(scan_cfg, projectile="He", mass_amu=4.0026_dp, &
                                  e_kev=e_kev, theta_deg=th_in, &
                                  ax_ang=2.84_dp, corrugation_ang=0.05_dp)

        call calc_gifad_transverse_kinematics(scan_cfg, e_perp, lambda_perp)
        th_rainbow = calc_gifad_rainbow_angle(scan_cfg)

        write(u_scan, '(5ES16.6)') th_in, e_perp, lambda_perp, e_kev, th_rainbow
    end do
    close(u_scan)
    print *, " Saved kinematic scaling scan to ex28_gifad_rainbow_scan.dat"
    print *, " Example 28 completed successfully."

end program ex28_grazing_fast_atom_diffraction
