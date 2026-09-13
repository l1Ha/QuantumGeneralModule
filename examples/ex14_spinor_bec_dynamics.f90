! ==============================================================================
! GeneralModule Example 14: Ultracold Spinor BEC Macroscopic Spin Dynamics
!
! Features:
!  1. Spin-1 (F=1) Bose-Einstein condensates: 87Rb (ferromagnetic) vs 23Na (polar)
!  2. Single-Mode Approximation (SMA) nonlinear spin-exchange dynamics
!  3. Quadratic Zeeman effect (QZE) tuning of dynamical instabilities
!  4. Coherent macroscopic spin oscillations and exact conservation of m_z & norm
!
! Standard: Fortran 2008
! ==============================================================================
program ex14_spinor_bec_dynamics
    use mod_constants, only: dp, PI, FS2AU, AU2FS
    use mod_spinor_bec
    implicit none

    type(spinor_param_t) :: rb87_param, na23_param
    type(spinor_state_t) :: init_state, state_rb
    real(dp) :: b_gauss, density_cm3, t_total_s, t_total_au, dt_s, dt_au
    integer  :: u_out, stat, i, n_steps
    real(dp) :: t_ms, p_p1, p_0, p_m1, mz, norm_val

    print *, "================================================================"
    print *, " Example 14: Ultracold Spin-1 BEC Macroscopic Spin Mixing       "
    print *, "================================================================"

    ! 1. Initialize 87Rb (Ferromagnetic) and 23Na (Antiferromagnetic) parameters
    b_gauss = 0.25_dp           ! Magnetic field 250 mG
    density_cm3 = 1.0e14_dp     ! Typical BEC condensate density 10^14 cm^-3

    call init_spinor_preset("87Rb", b_gauss, density_cm3, rb87_param, stat)
    call init_spinor_preset("23Na", b_gauss, density_cm3, na23_param, stat)

    write(*, '(A, ES14.4, A)') " 87Rb Spin-exchange coupling c_2: ", rb87_param%c2_au, " a.u. (Ferromagnetic)"
    write(*, '(A, ES14.4, A)') " 23Na Spin-exchange coupling c_2: ", na23_param%c2_au, " a.u. (Antiferromagnetic)"
    write(*, '(A, ES14.4, A)') " Quadratic Zeeman shift q_Z(B)  : ", rb87_param%q_zeeman_au, " a.u."
    print *, "----------------------------------------------------------------"

    ! 2. Prepare initial state with high m=0 population and small seeding:
    !    zeta_0 = sqrt(0.96), zeta_+1 = zeta_-1 = sqrt(0.02) -> m_z = 0
    init_state%zeta_p1 = cmplx(sqrt(0.02_dp), 0.0_dp, kind=dp)
    init_state%zeta_0  = cmplx(sqrt(0.96_dp), 0.0_dp, kind=dp)
    init_state%zeta_m1 = cmplx(sqrt(0.02_dp), 0.0_dp, kind=dp)
    init_state%pop_p1 = 0.02_dp
    init_state%pop_0  = 0.96_dp
    init_state%pop_m1 = 0.02_dp
    init_state%total_norm = 1.0_dp
    init_state%magnetization_mz = 0.0_dp

    ! 3. Evolve 87Rb spinor dynamics over 200 ms
    t_total_s = 0.200_dp               ! 200 ms
    ! 1 s in a.u. = 1.0 / (2.4188843265857e-17)
    t_total_au = t_total_s / 2.4188843265857e-17_dp
    n_steps = 200
    dt_s  = t_total_s / real(n_steps, dp)
    dt_au = t_total_au / real(n_steps, dp)

    open(newunit=u_out, file="spinor_bec_dynamics.dat", status="replace", action="write")
    write(u_out, '(A)') "# Spin-1 BEC Macroscopic Spin Mixing Dynamics (87Rb)"
    write(u_out, '(A)') "# Col 1: Time t (ms)"
    write(u_out, '(A)') "# Col 2: Population rho_0 (m = 0)"
    write(u_out, '(A)') "# Col 3: Population rho_+1 (m = +1)"
    write(u_out, '(A)') "# Col 4: Population rho_-1 (m = -1)"
    write(u_out, '(A)') "# Col 5: Magnetization m_z"
    write(u_out, '(A)') "# Col 6: Total Condensate Norm"

    state_rb = init_state
    do i = 0, n_steps
        t_ms = real(i, dp) * dt_s * 1000.0_dp
        p_0  = state_rb%pop_0
        p_p1 = state_rb%pop_p1
        p_m1 = state_rb%pop_m1
        mz   = state_rb%magnetization_mz
        norm_val = state_rb%total_norm

        write(u_out, '(6ES18.8)') t_ms, p_0, p_p1, p_m1, mz, norm_val

        if (i < n_steps) then
            call propagate_spinor_sma_rk4(rb87_param, dt_au, state_rb)
        end if
    end do
    close(u_out)
    print *, " Spinor BEC spin mixing time evolution written to: spinor_bec_dynamics.dat"

    print *, "----------------------------------------------------------------"
    print *, " Final Conservation Checks after 200 ms dynamics:"
    write(*, '(A, F14.10, A)') " Final total norm   : ", state_rb%total_norm, " (Strictly = 1.0)"
    write(*, '(A, ES14.6, A)') " Final magnetization: ", state_rb%magnetization_mz, " (Strictly = 0.0)"
    write(*, '(A, F10.4)')     " Final rho_0        : ", state_rb%pop_0
    write(*, '(A, F10.4)')     " Final rho_+1       : ", state_rb%pop_p1

    print *, "================================================================"
    print *, " Example 14 completed successfully.                             "
    print *, "================================================================"
end program ex14_spinor_bec_dynamics
