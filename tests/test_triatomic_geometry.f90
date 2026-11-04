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

    ! --------------------------------------------------------------------------
    ! Test 6: Analytic LEPS gradient vs central finite differences
    ! --------------------------------------------------------------------------
    block
        real(dp) :: g(3), g_col(2), fd(2), hh, v0, vp, vm
        real(dp) :: pts(3, 2), pt(2)
        integer :: k, i
        pts(1, :) = (/ 2.5_dp, 1.5_dp /)
        pts(2, :) = (/ 1.9_dp, 1.65_dp /)
        pts(3, :) = (/ 3.2_dp, 1.401_dp /)
        do k = 1, 3
            pt = pts(k, :)
            dist%r12 = pt(1); dist%r23 = pt(2); dist%r31 = pt(1) + pt(2)
            g = calc_leps_gradient(dist, leps_h3)
            g_col(1) = g(1) + g(3)     ! collinear chain rule (r31 = r12 + r23)
            g_col(2) = g(2) + g(3)
            hh = 1.0e-6_dp
            do i = 1, 2
                v0 = 0.0_dp; vp = 0.0_dp; vm = 0.0_dp
                dist%r12 = pt(1); dist%r23 = pt(2); dist%r31 = pt(1) + pt(2)
                if (i == 1) then
                    dist%r12 = pt(1) + hh; dist%r31 = dist%r12 + dist%r23
                    vp = calc_leps_potential(dist, leps_h3)
                    dist%r12 = pt(1) - hh; dist%r31 = dist%r12 + dist%r23
                    vm = calc_leps_potential(dist, leps_h3)
                else
                    dist%r23 = pt(2) + hh; dist%r31 = dist%r12 + dist%r23
                    vp = calc_leps_potential(dist, leps_h3)
                    dist%r23 = pt(2) - hh; dist%r31 = dist%r12 + dist%r23
                    vm = calc_leps_potential(dist, leps_h3)
                end if
                fd(i) = (vp - vm) / (2.0_dp * hh)
            end do
            n_total = n_total + 1
            if (abs(g_col(1) - fd(1)) < 1.0e-9_dp .and. abs(g_col(2) - fd(2)) < 1.0e-9_dp) then
                print *, " [PASS] LEPS analytic gradient vs FD at point ", k
                n_pass = n_pass + 1
            else
                print *, " [FAIL] LEPS gradient mismatch at point ", k, &
                         ": analytic", g_col, " FD", fd
            end if
        end do
    end block

    ! --------------------------------------------------------------------------
    ! Test 7: Intrinsic reaction coordinate (mass-weighted steepest descent)
    !   (a) both branches descend monotonically and reach the H2 asymptote;
    !   (b) the two branches are mirror images (symmetric H3);
    !   (c) the Wilson G-metric path deviates substantially from the
    !       unweighted-coordinate descent (the metric is not cosmetic).
    ! --------------------------------------------------------------------------
    block
        real(dp), allocatable :: s_a(:), r12_a(:), r23_a(:), v_a(:)
        real(dp), allocatable :: s_b(:), r12_b(:), r23_b(:), v_b(:)
        real(dp) :: masses(3), r_start(2), u(2), r_uw(2), g3(3), gcol(2), gn, du
        real(dp), allocatable :: r_uw_path(:, :)
        integer :: n_out, reason, i, n_uw, n_cmn
        real(dp) :: max_asc, dev_max

        masses = 1.0_dp
        u = (/ 1.0_dp, -1.0_dp /) / sqrt(2.0_dp)   ! antisymmetric stretch (H3 symmetry)
        r_start = (/ 1.774963_dp, 1.774963_dp /) + 1.0e-3_dp * u
        call intrinsic_reaction_coordinate(leps_h3, masses, r_start, 20000, 1.0e-2_dp, &
                                  1.0e-7_dp, 0.05_dp, 50.0_dp, &
                                  s_a, r12_a, r23_a, v_a, n_out, reason)

        n_total = n_total + 1
        if (reason > 0 .and. n_out > 100 .and. abs(v_a(n_out) + 0.1744_dp) < 5.0e-4_dp) then
            print *, " [PASS] IRC branch reaches the H2 asymptote: V_end = ", v_a(n_out), &
                     " after s = ", s_a(n_out)
            n_pass = n_pass + 1
        else
            print *, " [FAIL] IRC branch endpoint: reason =", reason, " V_end =", v_a(n_out)
        end if

        max_asc = 0.0_dp
        do i = 2, n_out
            max_asc = max(max_asc, v_a(i) - v_a(i - 1))
        end do
        n_total = n_total + 1
        if (max_asc < 1.0e-9_dp) then
            print *, " [PASS] IRC energy decreases monotonically (max step rise = ", max_asc, ")"
            n_pass = n_pass + 1
        else
            print *, " [FAIL] IRC ascends along the path: max rise = ", max_asc
        end if

        r_start = (/ 1.774963_dp, 1.774963_dp /) - 1.0e-3_dp * u
        call intrinsic_reaction_coordinate(leps_h3, masses, r_start, 20000, 1.0e-2_dp, &
                                  1.0e-7_dp, 0.05_dp, 50.0_dp, &
                                  s_b, r12_b, r23_b, v_b, n_out, reason)
        n_cmn = min(size(v_a), size(v_b))
        n_total = n_total + 1
        if (maxval(abs(v_a(1:n_cmn) - v_b(1:n_cmn))) < 1.0e-6_dp) then
            print *, " [PASS] IRC branches are mirror images (symmetric H3)"
            n_pass = n_pass + 1
        else
            print *, " [FAIL] Branch asymmetry: ", maxval(abs(v_a(1:n_cmn) - v_b(1:n_cmn)))
        end if

        ! (c) unweighted-coordinate descent deviates from the G-metric path
        r_uw = (/ 1.774963_dp, 1.774963_dp /) + 1.0e-3_dp * u
        allocate(r_uw_path(2, 20000))
        r_uw_path(:, 1) = r_uw
        n_uw = 1
        do i = 2, 20000
            dist%r12 = r_uw(1); dist%r23 = r_uw(2); dist%r31 = r_uw(1) + r_uw(2)
            g3 = calc_leps_gradient(dist, leps_h3)
            gcol(1) = g3(1) + g3(3)
            gcol(2) = g3(2) + g3(3)
            gn = sqrt(gcol(1)**2 + gcol(2)**2)
            if (gn < 1.0e-7_dp) exit
            r_uw = r_uw - 1.0e-2_dp * gcol / gn
            if (r_uw(1) <= 0.05_dp .or. r_uw(1) >= 50.0_dp .or. &
                r_uw(2) <= 0.05_dp .or. r_uw(2) >= 50.0_dp) exit
            n_uw = n_uw + 1
            r_uw_path(:, n_uw) = r_uw
        end do
        du = 0.0_dp
        n_cmn = min(size(s_a), n_uw)
        do i = 1, n_cmn
            du = max(du, sqrt((r12_a(i) - r_uw_path(1, i))**2 + (r23_a(i) - r_uw_path(2, i))**2))
        end do
        n_total = n_total + 1
        if (du > 0.1_dp) then
            print *, " [PASS] G-metric vs unweighted paths deviate: max = ", du, " a0"
            n_pass = n_pass + 1
        else
            print *, " [FAIL] Metric deviation too small: ", du
        end if
        deallocate(r_uw_path)
        if (allocated(s_a)) deallocate(s_a, r12_a, r23_a, v_a)
        if (allocated(s_b)) deallocate(s_b, r12_b, r23_b, v_b)
    end block

    print *, "--------------------------------------------------"
    print *, "Triatomic & PES Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All triatomic geometry and PES tests passed."
    else
        stop 1
    end if

end program test_triatomic_geometry
