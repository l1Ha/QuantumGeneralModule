! ==============================================================================
! GeneralModule: test_crossed_field_scattering.f90
!
! Unit Test Suite: Crossed Electric & Magnetic Fields Quantum Dynamics
!
! Standard: Fortran 2008
! ==============================================================================

program test_crossed_field_scattering
    use mod_constants, only: dp, PI
    use mod_crossed_field_scattering
    implicit none

    integer :: n_pass, n_total
    integer :: stat
    type(crossed_field_config_t) :: cfg
    type(crossed_field_state_t)  :: state
    real(dp), allocatable        :: h_par(:,:), h_tilt(:,:)
    integer                      :: dim_tot
    real(dp)                     :: theta_scan(5), energies_scan(5, 4)

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Crossed Fields Dynamics"
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Configuration Initialization
    ! --------------------------------------------------------------------------
    ! KRb molecule: B_rot = 1.114 GHz, d_0 = 0.574 Debye, E = 10 kV/cm, B = 100 G, theta = 45 deg
    call init_crossed_field_config(e_field_kv_cm=10.0_dp, b_field_gauss=100.0_dp, &
                                  theta_eb_deg=45.0_dp, rot_ghz=1.114_dp, &
                                  dipole_d=0.574_dp, j_max=2, cfg=cfg, stat=stat)

    n_total = n_total + 1
    if (stat == 0 .and. cfg%e_field_au > 0.0_dp .and. cfg%b_field_au > 0.0_dp) then
        print *, " [PASS] Crossed field configuration initialized successfully"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Configuration initialization failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Parallel Fields (theta = 0): Transverse Zeeman Decoupling
    ! --------------------------------------------------------------------------
    cfg%theta_eb_rad = 0.0_dp
    call build_crossed_field_hamiltonian(cfg, h_par, dim_tot, stat)

    n_total = n_total + 1
    ! For theta = 0, off-diagonal elements coupling M_S=+1/2 and M_S=-1/2 with same J, M must be strictly zero
    ! Element (1, 2) is (J=0, M=0, Ms=+1/2) and (J=0, M=0, Ms=-1/2)
    if (stat == 0 .and. abs(h_par(1, 2)) < 1.0e-14_dp) then
        print *, " [PASS] Parallel fields: Transverse coupling strictly zero (< 1e-14)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Parallel field decoupling violation: ", h_par(1, 2)
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Tilted Fields (theta = 45 deg): Active Transverse Coupling
    ! --------------------------------------------------------------------------
    cfg%theta_eb_rad = PI * 0.25_dp
    call build_crossed_field_hamiltonian(cfg, h_tilt, dim_tot, stat)

    n_total = n_total + 1
    if (stat == 0 .and. abs(h_tilt(1, 2)) > 1.0e-8_dp) then
        print *, " [PASS] Tilted fields (45 deg): Transverse coupling active = ", h_tilt(1, 2)
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Tilted field coupling missing"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Crossed Field Eigenstate Solver & Spectrum
    ! --------------------------------------------------------------------------
    call solve_crossed_field_eigenstates(cfg, state, stat)

    n_total = n_total + 1
    if (stat == 0 .and. state%dim_tot == 18) then
        print *, " [PASS] Crossed field states diagonalized: Hilbert dim = ", state%dim_tot
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Crossed field diagonalization failed"
    end if

    n_total = n_total + 1
    ! Spectrum must be monotonically ascending
    if (state%energies_au(1) <= state%energies_au(2) .and. &
        state%energies_au(2) <= state%energies_au(3)) then
        print *, " [PASS] Eigenenergies sorted: E_0 = ", state%energies_au(1), &
                 ", E_1 = ", state%energies_au(2)
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Energy sorting violation"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Electric Field Induced Orientation
    ! --------------------------------------------------------------------------
    ! Ground state should have non-zero orientation <cos(theta_R)> due to Stark mixing
    n_total = n_total + 1
    if (abs(state%orientation(1)) > 1.0e-3_dp) then
        print *, " [PASS] Ground state Stark orientation <cos(theta)> = ", state%orientation(1)
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Zero orientation in DC electric field"
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: Tilt Angle Scanning (0 to 90 degrees)
    ! --------------------------------------------------------------------------
    theta_scan = [0.0_dp, 22.5_dp, 45.0_dp, 67.5_dp, 90.0_dp]
    call scan_tilt_angle_spectrum(cfg, theta_scan, 5, energies_scan, stat)

    n_total = n_total + 1
    if (stat == 0) then
        print *, " [PASS] Tilt angle scan from 0 to 90 deg executed smoothly"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Tilt angle scan failed"
    end if

    print *, "--------------------------------------------------"
    print *, "Crossed Fields Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All crossed field dynamics tests passed."
    else
        stop 1
    end if

end program test_crossed_field_scattering
