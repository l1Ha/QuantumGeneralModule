! ==============================================================================
! GeneralModule: test_triatomic_geometry.f90
!
! Unit Test Suite: Triatomic Jacobi Coordinates, LEPS PES, & Geometric Phase
!
! Standard: Fortran 2008
! ==============================================================================

program test_triatomic_geometry
    use mod_constants, only: dp, PI, AMU2AU
    use mod_triatomic_geometry
    implicit none

    integer :: n_pass, n_total
    type(jacobi_coord_t)      :: jac_in, jac_out
    type(internuclear_dist_t) :: dist
    type(leps_param_t)        :: leps_h3
    type(conical_intersection_t) :: ci
    real(dp) :: v_asympt, v_saddle, e_low, e_up, gap, berry_phi

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Triatomic & PES       "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Jacobi <-> Internuclear Round-Trip Transformation
    ! --------------------------------------------------------------------------
    jac_in%r_diatom = 1.40_dp         ! r = 1.4 a0
    jac_in%r_atom_diatom = 3.50_dp    ! R = 3.5 a0
    jac_in%gamma_rad = PI / 3.0_dp    ! gamma = 60 deg
    jac_in%mass_a = 1.0_dp * AMU2AU
    jac_in%mass_b = 1.0_dp * AMU2AU
    jac_in%mass_c = 1.0_dp * AMU2AU

    call jacobi_to_internuclear(jac_in, dist)
    call internuclear_to_jacobi(dist, jac_in%mass_a, jac_in%mass_b, jac_in%mass_c, jac_out)

    n_total = n_total + 1
    if (abs(jac_in%r_diatom - jac_out%r_diatom) < 1.0e-11_dp .and. &
        abs(jac_in%r_atom_diatom - jac_out%r_atom_diatom) < 1.0e-11_dp .and. &
        abs(jac_in%gamma_rad - jac_out%gamma_rad) < 1.0e-11_dp) then
        print *, " [PASS] Jacobi <-> Internuclear round-trip transformation exact (< 1e-11)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Coordinate transformation round-trip error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Collinear Geometry Consistency
    ! --------------------------------------------------------------------------
    jac_in%gamma_rad = 0.0_dp  ! Collinear
    call jacobi_to_internuclear(jac_in, dist)

    n_total = n_total + 1
    ! For collinear A-B-C with equal masses, r12 + r23/2 = R => r12 = R - r/2
    if (abs(dist%r12 - (3.50_dp - 0.70_dp)) < 1.0e-12_dp) then
        print *, " [PASS] Collinear internuclear distance matches exact geometry: r12 = ", dist%r12
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Collinear internuclear distance error: ", dist%r12
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Benchmark H + H_2 -> H_3 LEPS Potential Energy Surface
    ! --------------------------------------------------------------------------
    call init_default_h3_leps(leps_h3)

    ! Asymptotic region: r23 = 1.401 a0 (equilibrium), r12 = r31 = 20.0 a0 (separated H + H2)
    dist%r23 = 1.401_dp
    dist%r12 = 20.0_dp
    dist%r31 = 20.0_dp
    v_asympt = calc_leps_potential(dist, leps_h3)

    n_total = n_total + 1
    ! Asymptotic potential should match diatom well depth -D_e = -0.1744 a.u.
    if (abs(v_asympt - (-0.1744_dp)) < 1.0e-4_dp) then
        print *, " [PASS] LEPS asymptotic potential matches isolated H2 well depth = ", v_asympt
        n_pass = n_pass + 1
    else
        print *, " [FAIL] LEPS asymptotic limit mismatch: ", v_asympt
    end if

    ! Symmetric collinear transition state saddle point: r12 = r23 = 1.76 a0, r31 = 3.52 a0
    dist%r12 = 1.76_dp
    dist%r23 = 1.76_dp
    dist%r31 = 3.52_dp
    v_saddle = calc_leps_potential(dist, leps_h3)

    n_total = n_total + 1
    ! Transition state has a higher energy than asymptotic valley (reaction barrier)
    if (v_saddle > v_asympt) then
        print *, " [PASS] LEPS reactive barrier exists: V_saddle = ", v_saddle, &
                 " > V_asympt = ", v_asympt
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Transition state barrier missing"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Conical Intersection (CI) Adiabatic Degeneracy
    ! --------------------------------------------------------------------------
    ci%x_ci = 1.5_dp
    ci%y_ci = 2.0_dp
    ci%kappa_tuning = 0.8_dp
    ci%lambda_coupl = 0.5_dp
    ci%e_ci = 0.25_dp

    ! At singularity (x_ci, y_ci), energy gap must strictly vanish
    call calc_conical_intersection_adiabats(ci, ci%x_ci, ci%y_ci, e_low, e_up, gap)

    n_total = n_total + 1
    if (abs(gap) < 1.0e-14_dp .and. abs(e_low - ci%e_ci) < 1.0e-14_dp) then
        print *, " [PASS] Conical intersection degeneracy exact: gap = ", gap
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Conical intersection degeneracy failure: ", gap
    end if

    ! Away from singularity, gap opens linearly
    call calc_conical_intersection_adiabats(ci, ci%x_ci + 0.1_dp, ci%y_ci, e_low, e_up, gap)
    n_total = n_total + 1
    if (gap > 0.1_dp) then
        print *, " [PASS] Linear cone splitting around CI: gap = ", gap
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Cone splitting too small"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Berry Geometric Phase around Conical Intersection (Phi = pi)
    ! --------------------------------------------------------------------------
    berry_phi = calc_berry_phase_around_ci(ci, radius=0.2_dp, n_steps=2000)

    n_total = n_total + 1
    if (abs(berry_phi - PI) < 2.0e-3_dp) then
        print *, " [PASS] Geometric Berry phase around CI = ", berry_phi, " (Exact: pi = ", PI, ")"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Berry phase mismatch: ", berry_phi, " vs ", PI
    end if

    print *, "--------------------------------------------------"
    print *, "Triatomic & PES Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All triatomic geometry and PES tests passed."
    else
        stop 1
    end if

end program test_triatomic_geometry
