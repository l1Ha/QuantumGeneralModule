! ==============================================================================
! GeneralModule: test_attosecond_transient_absorption.f90
!
! Unit Test Suite: Attosecond Transient Absorption Spectroscopy (ATAS)
!
! Standard: Fortran 2008
! ==============================================================================

program test_attosecond_transient_absorption
    use mod_constants, only: dp, PI
    use mod_attosecond_transient_absorption
    implicit none

    integer :: n_pass, n_total, stat
    type(atas_state_t) :: he_state
    real(dp) :: q0, q_eff_zero, q_eff_shifted
    real(dp) :: e_lis, t_beat
    real(dp) :: e_grid(15), tau_grid(11), spec_2d(15, 11)

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: ATAS & Autoionization  "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Helium 2s2p State Initialization
    ! --------------------------------------------------------------------------
    call init_atas_helium_benchmark(he_state, stat)
    n_total = n_total + 1
    if (stat == 0 .and. abs(he_state%energy_ev - 60.15_dp) < 1.0e-3_dp .and. &
        abs(he_state%gamma_ev - 0.037_dp) < 1.0e-4_dp) then
        print *, " [PASS] Helium 2s2p benchmark initialized: E0 = ", &
            he_state%energy_ev, " eV, Gamma = ", he_state%gamma_ev
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Helium benchmark initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Dynamic Fano q-Parameter Under Phase Perturbation
    ! --------------------------------------------------------------------------
    q0 = he_state%q_fano
    q_eff_zero = calc_laser_dressed_fano_q(q0, 0.0_dp)
    q_eff_shifted = calc_laser_dressed_fano_q(q0, 0.5_dp * PI)

    n_total = n_total + 1
    if (abs(q_eff_zero - q0) < 1.0e-12_dp .and. abs(q_eff_shifted - q0) > 0.1_dp) then
        print *, " [PASS] Dynamic Fano q: q(0) = ", q_eff_zero, ", dressed q(pi/2) = ", q_eff_shifted
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Laser-dressed Fano q failure: ", q_eff_zero, q_eff_shifted
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Light-Induced State (LIS) Energy Calculation
    ! --------------------------------------------------------------------------
    ! Bright at 60.15 eV, Dark at 58.60 eV, NIR photon = 1.55 eV (800 nm), Rabi = 0.1 eV
    e_lis = calc_light_induced_state_energy(60.15_dp, 58.60_dp, 1.55_dp, 0.10_dp)
    n_total = n_total + 1
    if (e_lis > 59.5_dp .and. e_lis < 60.5_dp) then
        print *, " [PASS] Light-Induced State (LIS) energy: E_LIS = ", e_lis, " eV"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected LIS energy: ", e_lis
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Quantum Beat Period for 1.0 eV State Splitting
    ! --------------------------------------------------------------------------
    t_beat = calc_quantum_beat_period_fs(1.0_dp)
    n_total = n_total + 1
    ! T = h / (1 eV) = 4.135667 fs
    if (abs(t_beat - 4.135667_dp) < 0.01_dp) then
        print *, " [PASS] Quantum beat period for Delta E = 1 eV: T_beat = ", t_beat, " fs (Theory: ~4.14 fs)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Quantum beat period error: ", t_beat
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: 2D ATAS Matrix Generation & Delay Asymmetry
    ! --------------------------------------------------------------------------
    call calc_atas_spectrum(he_state, nir_intensity_w_cm2=1.0e12_dp, nir_wavelength_nm=800.0_dp, &
                            n_energy=15, e_min_ev=59.8_dp, e_max_ev=60.5_dp, &
                            n_delay=11, tau_min_fs=-30.0_dp, tau_max_fs=30.0_dp, &
                            e_grid_ev=e_grid, tau_grid_fs=tau_grid, spec_2d=spec_2d)

    n_total = n_total + 1
    ! For negative delay (tau = -30 fs, index 1), perturbation is much smaller than near overlap (tau ~ 0 fs, index 6)
    if (abs(spec_2d(8, 1)) < abs(spec_2d(8, 6)) .and. abs(spec_2d(8, 6)) > 1.0e-5_dp) then
        print *, " [PASS] ATAS delay asymmetry verified: Delta OD(-30fs) = ", spec_2d(8, 1), &
                 " < Delta OD(0fs) = ", spec_2d(8, 6)
        n_pass = n_pass + 1
    else
        print *, " [FAIL] ATAS delay spectrum unexpected: OD(-30) = ", spec_2d(8, 1), ", OD(0) = ", spec_2d(8, 6)
    end if

    print *, "--------------------------------------------------"
    print *, "ATAS Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All ATAS and autoionization tests passed."
    else
        stop 1
    end if

end program test_attosecond_transient_absorption
