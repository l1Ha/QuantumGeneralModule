! ==============================================================================
! GeneralModule Example 27: Surface Electronic Friction & GLE Trajectory
!
! Scientific Engineering Application:
!  1. NO / Au(111) non-adiabatic gas-surface scattering dynamics
!  2. Local Density Friction Approximation (LDFA) electronic friction profile
!  3. Generalized Langevin Equation (GLE) trajectory and turning point
!  4. Non-adiabatic electron-hole pair energy loss Delta E_loss(E_incident)
!  5. Surface vibrational relaxation rate Gamma_vib and lifetime tau_vib
!
! Reference:
!   Head-Gordon & Tully, J. Chem. Phys. 103, 10137 (1995);
!   Wodtke et al., Science 290, 1585 (2000).
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================
program ex27_surface_electronic_friction_gle
    use mod_constants, only: dp, AU2ANG
    use mod_surface_electronic_friction
    implicit none

    type(metal_surface_t)          :: au111
    type(scattering_loss_result_t) :: res
    integer, parameter :: N_PTS = 200
    real(dp) :: t_fs(N_PTS), z_ang(N_PTS), v_ms(N_PTS)
    real(dp) :: m_no, e_in, gamma_vib, tau_ps, z_probe_bohr, z_probe_ang, eta_val
    integer  :: u_traj, u_loss, i, step

    print *, "================================================================"
    print *, " Example 27: Surface Electronic Friction & GLE Dissipation      "
    print *, "================================================================"

    ! 1. Initialize Au(111) metal surface: E_Fermi = 5.53 eV, peak friction = 1.2e-3 a.u.
    call init_metal_surface(au111, "Au(111)", 300.0_dp)

    m_no = 30.006_dp ! NO molecule mass [amu]
    write(*, '(A, A)')        " Metal Substrate             : ", trim(au111%name)
    write(*, '(A, F8.2, A)')  " Substrate Temperature       : ", au111%temp_k, " K"
    write(*, '(A, ES12.4, A)')" Peak Friction Coeff eta_0   : ", au111%eta0_au, " a.u."
    write(*, '(A, F8.2, A)')  " Projectile Mass (NO)        : ", m_no, " amu"
    print *, "----------------------------------------------------------------"

    ! 2. Calculate single trajectory at E_incident = 0.50 eV
    e_in = 0.50_dp
    call integrate_gle_scattering_trajectory(au111, mass_amu=m_no, e_incident_ev=e_in, &
                                            dt_fs=0.5_dp, n_steps=N_PTS, t_arr_fs=t_fs, &
                                            z_arr_ang=z_ang, v_arr_ms=v_ms, loss_res=res)

    write(*, '(A, F8.2, A)')  " Incident Energy E_in        : ", res%e_incident_ev, " eV"
    write(*, '(A, ES14.6, A)')" Non-adiabatic Energy Loss   : ", res%e_lost_ev, " eV"
    write(*, '(A, F10.4, A)') " Final Scattered Energy      : ", res%e_final_ev, " eV"
    write(*, '(A, F8.3, A)')  " Closest Approach Distance   : ", res%z_turnaround_ang, " A"

    open(newunit=u_traj, file="ex27_gle_trajectory.dat", status="replace", action="write")
    write(u_traj, '(A)') "# GeneralModule Example 27: NO/Au(111) GLE Collision Trajectory"
    write(u_traj, '(A)') "# E_in = 0.50 eV, dt = 0.5 fs"
    write(u_traj, '(A)') "# Col 1: Time (fs)"
    write(u_traj, '(A)') "# Col 2: Surface-Molecule Distance z (Angstrom)"
    write(u_traj, '(A)') "# Col 3: Velocity v_z (m/s)"

    do i = 1, N_PTS
        write(u_traj, '(3ES16.6)') t_fs(i), z_ang(i), v_ms(i)
    end do
    close(u_traj)
    print *, " Saved collision trajectory to ex27_gle_trajectory.dat"

    ! 3. Energy Loss vs Incident Energy & Vibrational Relaxation Rates
    open(newunit=u_loss, file="ex27_electronic_energy_loss.dat", &
         status="replace", action="write")
    write(u_loss, '(A)') "# GeneralModule Example 27: Electronic Friction Energy Loss Scan"
    write(u_loss, '(A)') "# Col 1: Incident energy E_in (eV)"
    write(u_loss, '(A)') "# Col 2: Electron-hole pair energy loss Delta E (eV)"
    write(u_loss, '(A)') "# Col 3: Distance z from surface (Bohr)"
    write(u_loss, '(A)') "# Col 4: Electronic friction coefficient eta(z) (a.u.)"
    write(u_loss, '(A)') "# Col 5: Vibrational relaxation rate Gamma_vib (ps^-1)"

    do step = 1, 50
        e_in = 0.05_dp + real(step, dp) * 0.03_dp     ! 0.08 to 1.55 eV
        z_probe_bohr = 1.0_dp + real(step, dp) * 0.15_dp ! 1.15 to 8.5 Bohr
        z_probe_ang  = z_probe_bohr * AU2ANG

        call integrate_gle_scattering_trajectory(au111, mass_amu=m_no, e_incident_ev=e_in, &
                                                dt_fs=0.5_dp, n_steps=N_PTS, t_arr_fs=t_fs, &
                                                z_arr_ang=z_ang, v_arr_ms=v_ms, loss_res=res)

        eta_val   = calc_electronic_friction_coeff(au111, z_probe_bohr)
        gamma_vib = calc_vibrational_relaxation_rate(au111, z_ang=z_probe_ang, mass_amu=m_no)

        write(u_loss, '(5ES16.6)') e_in, res%e_lost_ev, z_probe_bohr, eta_val, gamma_vib
    end do
    close(u_loss)

    ! Molecular vibrational lifetime near surface
    gamma_vib = calc_vibrational_relaxation_rate(au111, z_ang=2.0_dp * AU2ANG, mass_amu=m_no)
    tau_ps = 1.0_dp / max(1.0e-12_dp, gamma_vib)
    write(*, '(A, ES12.4, A)')" Vibrational Rate Gamma_vib  : ", gamma_vib, " ps^-1"
    write(*, '(A, F8.2, A)')  " Vibrational Lifetime tau    : ", tau_ps, " ps"
    print *, " Saved energy loss scan to ex27_electronic_energy_loss.dat"
    print *, " Example 27 completed successfully."

end program ex27_surface_electronic_friction_gle
