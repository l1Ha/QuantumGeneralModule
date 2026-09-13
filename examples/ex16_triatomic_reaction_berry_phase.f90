! ==============================================================================
! GeneralModule Example 16: Triatomic Reactive PES and Conical Intersection Berry Phase
!
! Features:
!  1. Standard benchmark H + H_2 -> H_3^\ddagger -> H_2 + H LEPS potential energy surface
!  2. Jacobi to internuclear coordinate transformations
!  3. Collinear minimum energy reaction path (MEP) and saddle-point barrier
!  4. Conical intersection adiabatic potential cones and topological Berry phase Phi = pi
!
! Standard: Fortran 2008
! ==============================================================================
program ex16_triatomic_reaction_berry_phase
    use mod_constants, only: dp, PI, AU2EV, AU2CM
    use mod_triatomic_geometry
    implicit none

    type(leps_param_t)           :: leps_h3
    type(internuclear_dist_t)    :: dist
    type(conical_intersection_t) :: ci
    real(dp) :: v_asympt, v_saddle, barrier_kcal_mol, berry_phase
    real(dp) :: r1, r2, v_pot, e_low, e_up, gap
    integer  :: u_out, i, j, n_grid

    print *, "================================================================"
    print *, " Example 16: Triatomic Reaction PES & Conical Intersection Phase"
    print *, "================================================================"

    ! 1. Initialize LEPS surface for H + H2 -> H2 + H
    call init_default_h3_leps(leps_h3)

    ! Asymptotic valley (isolated H2 + H)
    dist%r12 = 20.0_dp
    dist%r23 = 1.401_dp
    dist%r31 = 20.0_dp
    v_asympt = calc_leps_potential(dist, leps_h3)

    ! Collinear symmetric saddle point (Transition State H---H---H)
    dist%r12 = 1.760_dp
    dist%r23 = 1.760_dp
    dist%r31 = 3.520_dp
    v_saddle = calc_leps_potential(dist, leps_h3)

    ! 1 a.u. = 627.5095 kcal/mol
    barrier_kcal_mol = (v_saddle - v_asympt) * 627.5095_dp

    write(*, '(A, F10.4, A)') " Asymptotic H2 + H energy (well depth): ", v_asympt, " a.u."
    write(*, '(A, F10.4, A)') " Transition state H3 saddle energy    : ", v_saddle, " a.u."
    write(*, '(A, F10.2, A)') " Classical reaction barrier height    : ", barrier_kcal_mol, " kcal/mol"
    print *, "----------------------------------------------------------------"

    ! 2. Generate collinear PES 2D grid for contour plotting
    open(newunit=u_out, file="triatomic_reaction_leps.dat", status="replace", action="write")
    write(u_out, '(A)') "# Collinear H + H2 -> H2 + H LEPS Potential Energy Surface"
    write(u_out, '(A)') "# Col 1: r12 (a0)"
    write(u_out, '(A)') "# Col 2: r23 (a0)"
    write(u_out, '(A)') "# Col 3: V_LEPS (a.u., relative to asymptotic valley)"

    n_grid = 40
    do i = 1, n_grid
        r1 = 1.0_dp + real(i - 1, dp) * (3.5_dp / real(n_grid - 1, dp))
        do j = 1, n_grid
            r2 = 1.0_dp + real(j - 1, dp) * (3.5_dp / real(n_grid - 1, dp))
            dist%r12 = r1
            dist%r23 = r2
            dist%r31 = r1 + r2  ! Collinear
            v_pot = calc_leps_potential(dist, leps_h3) - v_asympt
            write(u_out, '(3ES18.8)') r1, r2, v_pot
        end do
        write(u_out, *) ""
    end do
    close(u_out)
    print *, " Collinear LEPS PES grid written to: triatomic_reaction_leps.dat"

    ! 3. Conical Intersection (CI) and Berry Geometric Phase
    print *, "----------------------------------------------------------------"
    print *, " Conical Intersection (CI) and Topological Geometric Phase:"
    ci%x_ci = 0.0_dp
    ci%y_ci = 0.0_dp
    ci%kappa_tuning = 0.5_dp
    ci%lambda_coupl = 0.5_dp
    ci%e_ci = 0.10_dp

    call calc_conical_intersection_adiabats(ci, 0.0_dp, 0.0_dp, e_low, e_up, gap)
    write(*, '(A, F10.4, A)') " Energy degeneracy at CI center: Gap = ", gap, " a.u. (Exact zero)"

    berry_phase = calc_berry_phase_around_ci(ci, radius=0.25_dp, n_steps=1000)
    write(*, '(A, F10.6, A, F10.6, A)') " Topological Berry Phase Phi_B  : ", berry_phase, &
                                       " rad (Exact pi = ", PI, " rad)"

    print *, "================================================================"
    print *, " Example 16 completed successfully.                             "
    print *, "================================================================"
end program ex16_triatomic_reaction_berry_phase
