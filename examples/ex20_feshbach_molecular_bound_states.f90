! ==============================================================================
! GeneralModule Example 20: Magnetic & Optical Feshbach Molecular States
!
! Features:
!  1. 6Li broad magnetic Feshbach resonance (B0 = 832.2 G) scattering length a(B)
!  2. Universal halo dimer vs coupled-channel molecular binding energy Eb(B)
!  3. Hellmann-Feynman closed-channel fraction Z(B) across the resonance
!  4. 87Rb Optical Feshbach Resonance (OFR) dispersive scattering length tuning
!  5. Photostimulated two-body inelastic loss rate coefficient K_2(Delta_L)
!
! Standard: Fortran 2008
! ==============================================================================
program ex20_feshbach_molecular_bound_states
    use mod_constants, only: dp, AU2EV, PI
    use mod_feshbach_bound_states
    implicit none

    type(mfr_param_t) :: li6_mfr
    type(ofr_param_t) :: rb87_ofr
    real(dp) :: b_field, a_scat, eb_univ, eb_coup, z_frac
    real(dp) :: delta_hz, a_re, a_im, k2_loss
    integer  :: u_out, i, stat

    print *, "================================================================"
    print *, " Example 20: Magnetic & Optical Feshbach Molecular States       "
    print *, "================================================================"

    ! 1. Initialize 6Li magnetic resonance (832 G)
    call init_mfr_preset("6Li", li6_mfr, stat)

    write(*, '(A, F8.2, A)') " 6Li resonance pole B0         : ", li6_mfr%b0_gauss, " G"
    write(*, '(A, F8.1, A)') " Resonance width Delta B       : ", li6_mfr%delta_b_gauss, " G"
    write(*, '(A, F8.1, A)') " Background scattering length  : ", li6_mfr%a_bg_au, " a0"
    write(*, '(A, F8.2)')    " Resonance strength s_res      : ", li6_mfr%s_res
    write(*, '(A, F8.4, A)') " Characteristic length R*      : ", li6_mfr%r_star_au, " a0"
    print *, "----------------------------------------------------------------"

    ! 2. Scan MFR across BEC side (700 G to 830 G) where bound state exists
    open(newunit=u_out, file="ex20_mfr_bound_states.dat", status="replace", action="write")
    write(u_out, '(A)') "# GeneralModule Example 20: 6Li Magnetic Feshbach Molecular Bound States"
    write(u_out, '(A)') "# Col 1: Magnetic Field B (Gauss)"
    write(u_out, '(A)') "# Col 2: Scattering Length a(B) (Bohr)"
    write(u_out, '(A)') "# Col 3: Universal Binding Energy Eb_univ (eV)"
    write(u_out, '(A)') "# Col 4: Coupled-Channel Binding Energy Eb_coupled (eV)"
    write(u_out, '(A)') "# Col 5: Closed-Channel Fraction Z(B)"

    write(*, '(A)') " B (Gauss) |   a(B) (a0)   |  Eb_coup (eV)  |   Eb_univ (eV) |     Z(B)"
    write(*, '(A)') "-----------+---------------+----------------+----------------+--------------"

    do i = 1, 27
        b_field = 700.0_dp + real(i - 1, dp) * 5.0_dp
        a_scat  = calc_mfr_scattering_length(li6_mfr, b_field)
        eb_univ = calc_mfr_bound_energy_universal(li6_mfr, b_field)
        eb_coup = calc_mfr_bound_energy_coupled(li6_mfr, b_field)
        z_frac  = calc_mfr_closed_channel_fraction(li6_mfr, b_field)

        write(u_out, '(5ES16.6)') b_field, a_scat, eb_univ * AU2EV, eb_coup * AU2EV, z_frac
        if (mod(i, 5) == 1) then
            write(*, '(F10.1, " | ", ES13.4, " | ", ES14.4, " | ", ES14.4, " | ", ES12.4)') &
                b_field, a_scat, eb_coup * AU2EV, eb_univ * AU2EV, z_frac
        end if
    end do
    close(u_out)

    ! 3. Optical Feshbach Resonance (OFR) in 87Rb: laser detuning scan
    call init_ofr_param("87Rb", 86.91_dp, 100.0_dp, 1.0e7_dp, 50.0_dp, rb87_ofr, stat)

    open(newunit=u_out, file="ex20_ofr_loss.dat", status="replace", action="write")
    write(u_out, '(A)') "# Optical Feshbach Resonance (OFR) Dispersive Scattering & Loss"
    write(u_out, '(A)') "# Col 1: Laser Detuning Delta_L (MHz)"
    write(u_out, '(A)') "# Col 2: Real Scattering Length Re(a) (Bohr)"
    write(u_out, '(A)') "# Col 3: Scattering Length Shift delta_a (Bohr)"
    write(u_out, '(A)') "# Col 4: Two-Body Inelastic Loss Rate K_2 (cm^3 / s)"

    do i = 1, 41
        delta_hz = (-20.0_dp + real(i - 1, dp) * 1.0_dp) * 1.0e6_dp
        call calc_ofr_complex_scattering_length(rb87_ofr, delta_hz, a_re, a_im)
        k2_loss  = calc_ofr_inelastic_loss_rate(rb87_ofr, delta_hz)

        write(u_out, '(4ES16.6)') delta_hz / 1.0e6_dp, a_re, a_re - rb87_ofr%a_bg_au, k2_loss
    end do
    close(u_out)

    print *, "----------------------------------------------------------------"
    print *, " Output data saved to ex20_mfr_bound_states.dat and ex20_ofr_loss.dat"
    print *, " Example 20 completed successfully."

end program ex20_feshbach_molecular_bound_states
