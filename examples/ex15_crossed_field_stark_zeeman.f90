! ==============================================================================
! GeneralModule Example 15: Molecular Orientation in Crossed Electric & Magnetic Fields
!
! Features:
!  1. Polar paramagnetic molecule in non-collinear static E and B fields
!  2. Arbitrary tilt angle theta_EB between DC electric and magnetic fields
!  3. Off-diagonal Delta_M = +/- 1 transitions driven by transverse B_x and E_x
!  4. Field-induced molecular orientation <cos theta> and Stark-Zeeman avoided crossings
!
! Standard: Fortran 2008
! ==============================================================================
program ex15_crossed_field_stark_zeeman
    use mod_constants, only: dp, PI
    use mod_crossed_field_scattering
    implicit none

    type(crossed_field_config_t) :: cfg
    type(crossed_field_state_t)  :: state
    real(dp) :: e_kv_cm, b_gauss, rot_ghz, dipole_d
    integer  :: j_max, stat, u_out, i, n_angles
    real(dp), allocatable :: angles_deg(:), energies_au(:), orientations(:)

    print *, "================================================================"
    print *, " Example 15: Molecular Orientation in Crossed E x B Fields     "
    print *, "================================================================"

    ! 1. Initialize KRb polar paramagnetic molecule parameters
    !    B_rot = 1.11395 GHz, d = 0.574 Debye
    e_kv_cm  = 12.0_dp         ! 12 kV/cm
    b_gauss  = 1000.0_dp       ! 1000 Gauss (0.1 T)
    rot_ghz  = 1.11395_dp
    dipole_d = 0.574_dp
    j_max    = 2               ! J = 0, 1, 2 -> 9 spatial states x 2 spin = 18 states

    ! 2. Solve parallel fields (theta_EB = 0 deg)
    call init_crossed_field_config(e_kv_cm, b_gauss, 0.0_dp, rot_ghz, dipole_d, j_max, cfg, stat)
    call solve_crossed_field_eigenstates(cfg, state, stat)
    call calc_crossed_field_observables(cfg, state)

    write(*, '(A, I4)')       " Total Hilbert space dimension   : ", state%dim_tot
    write(*, '(A, F10.4, A)') " Parallel fields Ground Energy    : ", state%energies_au(1) * 1.0e6_dp, " micro-a.u."
    write(*, '(A, F10.4)')   " Parallel fields Ground <cos theta>: ", state%orientation(1)
    print *, "----------------------------------------------------------------"

    ! 3. Scan tilt angle theta_EB from 0 to 90 degrees
    n_angles = 46
    allocate(angles_deg(n_angles), energies_au(n_angles), orientations(n_angles))

    open(newunit=u_out, file="crossed_field_orientation.dat", status="replace", action="write")
    write(u_out, '(A)') "# Crossed Field Stark-Zeeman Angle Scan (KRb molecule)"
    write(u_out, '(A)') "# Col 1: Tilt Angle theta_EB (degrees)"
    write(u_out, '(A)') "# Col 2: Ground State Energy E_0 (micro-a.u.)"
    write(u_out, '(A)') "# Col 3: Orientation <cos theta>"

    do i = 1, n_angles
        angles_deg(i) = real(i - 1, dp) * (90.0_dp / real(n_angles - 1, dp))
        cfg%theta_eb_rad = angles_deg(i) * (PI / 180.0_dp)
        call solve_crossed_field_eigenstates(cfg, state, stat)
        call calc_crossed_field_observables(cfg, state)

        energies_au(i)  = state%energies_au(1)
        orientations(i) = state%orientation(1)
        write(u_out, '(3ES18.8)') angles_deg(i), energies_au(i) * 1.0e6_dp, orientations(i)
    end do
    close(u_out)
    print *, " Crossed-field orientation scan written to: crossed_field_orientation.dat"

    print *, "----------------------------------------------------------------"
    write(*, '(A, F10.4)') " Orientation at theta =  0 deg (Parallel)     : ", orientations(1)
    write(*, '(A, F10.4)') " Orientation at theta = 45 deg (Tilted)       : ", orientations(23)
    write(*, '(A, F10.4)') " Orientation at theta = 90 deg (Perpendicular): ", orientations(46)

    print *, "================================================================"
    print *, " Example 15 completed successfully.                             "
    print *, "================================================================"
end program ex15_crossed_field_stark_zeeman
