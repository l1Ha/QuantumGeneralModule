! ==============================================================================
! GeneralModule Example 13: Autoionization and Fano Resonances in Atoms
!
! Features:
!  1. Configuration interaction between discrete quasi-bound state and continuum
!  2. Asymmetric Fano absorption lineshape sigma(E) with q-parameter
!  3. Complex Coordinate Rotation (CCR) to extract resonance pole E_R - i*Gamma/2
!  4. Time-domain autoionization decay and electron emission flux J(t)
!
! Standard: Fortran 2008
! ==============================================================================
program ex13_autoionization_fano_resonance
    use mod_constants, only: dp, PI, FS2AU, AU2FS
    use mod_autoionization_fano
    implicit none

    type(fano_profile_t)  :: fano
    type(ccr_resonance_t) :: ccr_res
    real(dp) :: e_min, e_max, energy, sigma, e_anti, e_peak, sig_peak
    real(dp) :: t_fs, t_au, p_surv, flux_em, tau_au, tau_fs
    integer  :: u_out, i, n_pts, stat

    print *, "================================================================"
    print *, " Example 13: Autoionization and Fano Resonance Spectroscopy     "
    print *, "================================================================"

    ! 1. Initialize Helium-like 2s2p 1P^o autoionizing state parameters:
    !    Resonance energy E_0 = 2.220 a.u. (~60.4 eV above ground state)
    !    Decay width Gamma = 0.00137 a.u. (~37 meV)
    !    Fano asymmetry parameter q = -2.80
    fano%e_resonance_au = 2.220_dp
    fano%gamma_width_au = 0.00137_dp
    fano%q_parameter    = -2.80_dp
    fano%sigma_0_au     = 1.0_dp

    call calc_autoionization_lifetime(fano%gamma_width_au, tau_au, tau_fs, stat)
    fano%lifetime_fs = tau_fs

    write(*, '(A, F10.4, A)') " Autoionizing resonance energy E_0: ", fano%e_resonance_au, " a.u."
    write(*, '(A, ES12.4, A)') " Autoionization decay width Gamma : ", fano%gamma_width_au, " a.u."
    write(*, '(A, F10.2)')   " Fano asymmetry shape parameter q: ", fano%q_parameter
    write(*, '(A, F10.2, A)') " Autoionization lifetime tau     : ", fano%lifetime_fs, " fs"
    print *, "----------------------------------------------------------------"

    ! 2. Compute Fano profile key positions
    e_anti   = fano%e_resonance_au - fano%q_parameter * (0.5_dp * fano%gamma_width_au)
    e_peak   = fano%e_resonance_au + (0.5_dp * fano%gamma_width_au) / fano%q_parameter
    sig_peak = fano%sigma_0_au * (1.0_dp + fano%q_parameter**2)

    write(*, '(A, F10.5, A)') " Anti-resonance zero point E_zero : ", e_anti, " a.u."
    write(*, '(A, F10.5, A)') " Resonant absorption maximum E_max: ", e_peak, " a.u."
    write(*, '(A, F10.2, A)') " Peak cross section enhancement   : ", sig_peak, " sigma_0"

    open(newunit=u_out, file="autoionization_fano_resonance.dat", status="replace", action="write")
    write(u_out, '(A)') "# Autoionization and Fano Absorption Profile"
    write(u_out, '(A)') "# Col 1: Energy E (a.u.)"
    write(u_out, '(A)') "# Col 2: Reduced energy epsilon = 2*(E - E_0) / Gamma"
    write(u_out, '(A)') "# Col 3: Photoabsorption cross section sigma(E) / sigma_0"
    write(u_out, '(A)') "# Col 4: Time t (fs)"
    write(u_out, '(A)') "# Col 5: Survival probability P(t) = exp(-t/tau)"
    write(u_out, '(A)') "# Col 6: Electron emission flux J(t) (fs^-1)"

    n_pts = 100
    e_min = fano%e_resonance_au - 6.0_dp * fano%gamma_width_au
    e_max = fano%e_resonance_au + 6.0_dp * fano%gamma_width_au

    do i = 1, n_pts
        energy = e_min + real(i - 1, dp) * ((e_max - e_min) / real(n_pts - 1, dp))
        sigma  = calc_fano_profile(fano, energy)

        ! Time domain from 0 to 50 fs
        t_fs = real(i - 1, dp) * (50.0_dp / real(n_pts - 1, dp))
        t_au = t_fs * FS2AU
        call calc_time_domain_autoionization_decay(fano%gamma_width_au, t_au, p_surv, flux_em)

        write(u_out, '(6ES18.8)') energy, 2.0_dp * (energy - fano%e_resonance_au) / fano%gamma_width_au, &
                                  sigma, t_fs, p_surv, flux_em / AU2FS
    end do
    close(u_out)
    print *, " Fano spectroscopy and time-domain decay written to: autoionization_fano_resonance.dat"

    ! 3. Complex Coordinate Rotation (CCR) Verification
    print *, "----------------------------------------------------------------"
    print *, " Complex Coordinate Rotation (CCR) Resonance Extraction:"
    call solve_ccr_resonance_model(e_bound_0=2.22_dp, e_cont_0=2.20_dp, v_coupl=0.015_dp, &
                                   theta_rad=0.30_dp, res=ccr_res, stat=stat)
    write(*, '(A, F10.4, A, ES12.4, A)') " CCR Pole Energy: E = ", ccr_res%e_r_au, &
                                         " - i*", ccr_res%gamma_au / 2.0_dp, " a.u."
    write(*, '(A, F10.2, A)') " Extracted resonance lifetime tau: ", ccr_res%lifetime_fs, " fs"

    print *, "================================================================"
    print *, " Example 13 completed successfully.                             "
    print *, "================================================================"
end program ex13_autoionization_fano_resonance
