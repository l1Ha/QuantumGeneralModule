! ==============================================================================
! GeneralModule Example 21: Attosecond Transient Absorption Spectroscopy (ATAS)
!
! Features:
!  1. Helium 2s2p (^1P) classic autoionizing state benchmark (60.15 eV)
!  2. Phase-Perturbation Model (PPM) & laser-dressed dynamic Fano q(tau)
!  3. Light-Induced States (LIS) and quantum beat periods
!  4. 2D transient absorption spectrogram Delta OD(omega, tau) generation
!
! Standard: Fortran 2008
! ==============================================================================
program ex21_attosecond_transient_absorption
    use mod_constants, only: dp, PI
    use mod_attosecond_transient_absorption
    implicit none

    type(atas_state_t) :: he_state
    integer  :: u_out, ie, it, stat
    real(dp) :: e_lis, t_beat
    integer,  parameter :: N_E = 31, N_TAU = 25
    real(dp) :: e_grid(N_E), tau_grid(N_TAU), spec_2d(N_E, N_TAU)

    print *, "================================================================"
    print *, " Example 21: Attosecond Transient Absorption Spectroscopy (ATAS)"
    print *, "================================================================"

    ! 1. Initialize Helium 2s2p benchmark
    call init_atas_helium_benchmark(he_state, stat)

    write(*, '(A, F8.2, A)') " Autoionizing Resonance Energy E0: ", he_state%energy_ev, " eV"
    write(*, '(A, F8.4, A)') " Autoionization Width Gamma      : ", he_state%gamma_ev, " eV"
    write(*, '(A, F8.2)')    " Field-free Fano q-parameter     : ", he_state%q_fano
    write(*, '(A, F8.2, A)') " Autoionization Lifetime (hbar/G): ", (0.6582119_dp / he_state%gamma_ev), " fs"
    print *, "----------------------------------------------------------------"

    ! 2. Light-Induced State (LIS) and Quantum Beat
    e_lis = calc_light_induced_state_energy(he_state%energy_ev, 58.60_dp, 1.55_dp, 0.15_dp)
    t_beat = calc_quantum_beat_period_fs(abs(he_state%energy_ev - 58.60_dp))

    write(*, '(A, F8.3, A)') " Light-Induced State (LIS) Energy: ", e_lis, " eV"
    write(*, '(A, F8.3, A)') " Quantum Beat Period (Delta E)   : ", t_beat, " fs"
    print *, "----------------------------------------------------------------"

    ! 3. Compute 2D ATAS Spectrogram Delta OD(omega, tau)
    call calc_atas_spectrum(he_state, nir_intensity_w_cm2=2.0e12_dp, nir_wavelength_nm=800.0_dp, &
                            n_energy=N_E, e_min_ev=59.5_dp, e_max_ev=60.8_dp, &
                            n_delay=N_TAU, tau_min_fs=-30.0_dp, tau_max_fs=30.0_dp, &
                            e_grid_ev=e_grid, tau_grid_fs=tau_grid, spec_2d=spec_2d)

    ! 4. Save 2D Matrix Data
    open(newunit=u_out, file="ex21_atas_spectrogram.dat", status="replace", action="write")
    write(u_out, '(A)') "# GeneralModule Example 21: 2D ATAS Spectrogram Delta OD(omega, tau)"
    write(u_out, '(A)') "# Col 1: Photon Energy E (eV)"
    write(u_out, '(A)') "# Col 2: Time Delay tau (fs)"
    write(u_out, '(A)') "# Col 3: Differential Optical Density Delta OD"

    write(*, '(A)') " tau (fs) | E = 59.8 eV | E = 60.15 eV (Peak) | E = 60.5 eV"
    write(*, '(A)') "----------+-------------+---------------------+------------"

    do it = 1, N_TAU
        if (mod(it, 4) == 1) then
            write(*, '(F9.2, " | ", ES11.3, " | ", ES19.3, " | ", ES10.3)') &
                tau_grid(it), spec_2d(8, it), spec_2d(16, it), spec_2d(24, it)
        end if
        do ie = 1, N_E
            write(u_out, '(3ES16.6)') e_grid(ie), tau_grid(it), spec_2d(ie, it)
        end do
        write(u_out, *) ""
    end do
    close(u_out)

    print *, "----------------------------------------------------------------"
    print *, " SUCCESS: ATAS 2D spectrogram saved to ex21_atas_spectrogram.dat"
    print *, "================================================================"

end program ex21_attosecond_transient_absorption
