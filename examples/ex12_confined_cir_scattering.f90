! ==============================================================================
! GeneralModule Example 12: Confinement-Induced Resonances (CIR) in 1D Waveguides
!
! Features:
!  1. Quasi-1D atomic waveguide with transverse trap frequency omega_perp
!  2. Olshanii CIR resonance pole extraction: a_CIR = 0.9684 * a_perp
!  3. Effective 1D coupling strength g_1D(a_s) and molecular binding energy E_b
!  4. Tonks-Girardeau strongly-correlated gas regime parameter gamma_LL
!
! Standard: Fortran 2008
! ==============================================================================
program ex12_confined_cir_scattering
    use mod_constants, only: dp, PI, TWOPI, AMU2AU
    use mod_confined_scattering
    implicit none

    real(dp) :: mu_rb87, mass_rb87, omega_perp_si, omega_perp_au, a_perp, a_cir, a_s, e_b, gamma_ll
    real(dp) :: n_1d_au
    type(waveguide_1d_t) :: wg
    type(cir_result_t)   :: res
    integer :: u_out, i, n_pts, stat

    print *, "================================================================"
    print *, " Example 12: Confinement-Induced Resonances (CIR) in 1D Guides  "
    print *, "================================================================"

    ! 1. Initialize 87Rb + 87Rb in transverse trap omega_perp = 2*pi * 20 kHz
    mass_rb87 = 86.90918_dp * AMU2AU
    mu_rb87   = 0.5_dp * mass_rb87
    omega_perp_si = TWOPI * 20000.0_dp  ! 20 kHz in rad/s
    ! Convert rad/s to a.u.: 1 s = 1.0 / (2.4188843265857e-17) a.u.
    omega_perp_au = omega_perp_si * 2.4188843265857e-17_dp

    call init_waveguide_1d(omega_perp_au, mu_rb87, wg, stat)

    a_perp = wg%a_perp_au
    a_cir  = wg%a_cir_au

    write(*, '(A, F10.2, A)') " Transverse oscillator length a_perp: ", a_perp, " a_0"
    write(*, '(A, F10.2, A)') " Olshanii CIR resonance pole a_CIR  : ", a_cir, " a_0"
    write(*, '(A, F10.5)')   " Ratio a_CIR / a_perp               : ", a_cir / a_perp
    print *, "----------------------------------------------------------------"

    ! 2. Open output file for CIR dispersion curve
    open(newunit=u_out, file="confined_cir_scattering.dat", status="replace", action="write")
    write(u_out, '(A)') "# Confinement-Induced Resonance in 1D Waveguide"
    write(u_out, '(A)') "# Col 1: 3D Scattering length a_s (a_0)"
    write(u_out, '(A)') "# Col 2: Ratio a_s / a_perp"
    write(u_out, '(A)') "# Col 3: 1D Effective coupling g_1D (a.u.)"
    write(u_out, '(A)') "# Col 4: Dimer Binding Energy E_b (a.u.)"
    write(u_out, '(A)') "# Col 5: Lieb-Liniger parameter gamma_LL"

    ! Density: 1 atom per micrometer = 1.0 / (1.0e-6 m / 5.291772109e-11 m/a0) = 5.29177e-5 a.u.^-1
    n_1d_au = 5.291772109e-5_dp
    n_pts = 80

    do i = 1, n_pts
        ! Scan a_s from 100 to 2500 a_0 around a_cir
        a_s = 100.0_dp + real(i - 1, dp) * (2400.0_dp / real(n_pts - 1, dp))
        call calc_olshanii_cir_parameters(wg, a_s, n_1d_au, res, stat)
        e_b  = calc_confined_dimer_binding_energy(wg, a_s)
        gamma_ll = calc_lieb_liniger_parameter(res%g_1d_au, n_1d_au, mass_rb87)

        write(u_out, '(5ES18.8)') a_s, a_s / a_perp, res%g_1d_au, e_b, gamma_ll
    end do
    close(u_out)
    print *, " Confinement-induced resonance data written to: confined_cir_scattering.dat"

    ! 3. Summary of Tonks-Girardeau crossover
    print *, "----------------------------------------------------------------"
    print *, " Lieb-Liniger Correlation Parameter at Key Scattering Lengths:"
    a_s = 100.0_dp
    call calc_olshanii_cir_parameters(wg, a_s, n_1d_au, res, stat)
    gamma_ll = calc_lieb_liniger_parameter(res%g_1d_au, n_1d_au, mass_rb87)
    write(*, '(A, F8.1, A, ES14.4, A)') " Weak interaction (a_s = ", a_s, " a0): gamma_LL = ", gamma_ll, " (Mean Field)"

    a_s = a_cir * 0.99_dp
    call calc_olshanii_cir_parameters(wg, a_s, n_1d_au, res, stat)
    gamma_ll = calc_lieb_liniger_parameter(res%g_1d_au, n_1d_au, mass_rb87)
    write(*, '(A, F8.1, A, ES14.4, A)') " Near CIR pole     (a_s = ", a_s, " a0): gamma_LL = ", gamma_ll, " (Tonks-Girardeau Gas)"

    print *, "================================================================"
    print *, " Example 12 completed successfully.                             "
    print *, "================================================================"
end program ex12_confined_cir_scattering
