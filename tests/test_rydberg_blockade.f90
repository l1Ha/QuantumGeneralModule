! ==============================================================================
! GeneralModule: test_rydberg_blockade.f90
!
! Unit Test Suite: Rydberg Blockade and Many-Body Quantum Dynamics
!
! Standard: Fortran 2008
! ==============================================================================

program test_rydberg_blockade
    use mod_constants, only: dp
    use mod_rydberg_blockade
    implicit none

    integer :: n_pass, n_total
    type(rydberg_atom_t) :: rb70
    type(rydberg_array_config_t) :: array_cfg
    real(dp) :: r_b, o_z2
    real(dp) :: occ(4)
    real(dp) :: t_arr(100), pg(100), ps(100), pd(100)
    real(dp) :: z2_arr(100)
    real(dp) :: max_pd, max_ps, norm_err

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Rydberg Blockade       "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Rydberg Atom Initialization & n^11 Scaling
    ! --------------------------------------------------------------------------
    call init_rydberg_atom(rb70, "87Rb", n_principal=70, l_orbital=0)

    n_total = n_total + 1
    if (rb70%c6_mhz_um6 > 8.0e5_dp .and. rb70%c6_mhz_um6 < 9.5e5_dp .and. &
        rb70%lifetime_us > 100.0_dp) then
        print *, " [PASS] 87Rb 70S initialized: C6/h = ", rb70%c6_mhz_um6, &
                 " MHz*um^6, lifetime = ", rb70%lifetime_us, " us"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Rydberg atom initialization error: C6 = ", rb70%c6_mhz_um6
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Rydberg Blockade Radius R_b
    ! --------------------------------------------------------------------------
    ! For Rabi frequency Omega / 2pi = 2.0 MHz
    r_b = calc_rydberg_blockade_radius(rb70, rabi_mhz=2.0_dp)

    n_total = n_total + 1
    ! Analytical: (8.62e5 / 2)^(1/6) = 8.67 um
    if (abs(r_b - 8.67_dp) < 0.20_dp) then
        print *, " [PASS] Rydberg blockade radius verified: R_b = ", r_b, " um (Expected ~8.67 um)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected blockade radius: R_b = ", r_b
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Two-Atom Blockade Dynamics (Suppression of |rr>)
    ! --------------------------------------------------------------------------
    ! Spacing a = 3.0 um << R_b (8.67 um), deeply in the blockade regime
    call calc_two_atom_dynamics(rb70, spacing_um=3.0_dp, rabi_mhz=2.0_dp, detuning_mhz=0.0_dp, &
                                t_max_us=1.0_dp, n_steps=100, t_arr=t_arr, &
                                p_g=pg, p_single=ps, p_double=pd)

    max_pd = maxval(pd)
    max_ps = maxval(ps)
    norm_err = maxval(abs(pg + ps + pd - 1.0_dp))

    n_total = n_total + 1
    ! In strong blockade: max(p_double) << 1, max(p_single) > 0.8, norm conserved
    if (max_pd < 0.05_dp .and. max_ps > 0.80_dp .and. norm_err < 1.0e-5_dp) then
        print *, " [PASS] Two-atom blockade dynamics verified:"
        print *, "        Max |rr> double excitation = ", max_pd, " (suppressed!)"
        print *, "        Max single excitation W-state = ", max_ps
        print *, "        Probability norm error = ", norm_err
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Blockade dynamics failed: max_pd = ", max_pd, " max_ps = ", max_ps
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Staggered Z2 Order Parameter
    ! --------------------------------------------------------------------------
    ! Test with Neel state |r g r g>: occ = [1, 0, 1, 0]
    ! O_Z2 = (-1*1 + 1*0 - 1*1 + 1*0) / 4 = -0.5
    occ = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    o_z2 = calc_z2_order_parameter(4, occ)

    n_total = n_total + 1
    if (abs(o_z2 - (-0.5_dp)) < 1.0e-12_dp) then
        print *, " [PASS] Staggered Z2 order parameter verified: O_Z2(|rgrg>) = ", o_z2
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Z2 order parameter error: ", o_z2
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Quantum Many-Body Scar Revivals
    ! --------------------------------------------------------------------------
    call init_rydberg_array(array_cfg, n_atoms=8, spacing_um=5.0_dp, &
                            rabi_mhz=2.0_dp, detuning_mhz=0.0_dp, boundary_cond=2)
    call calc_rydberg_scar_dynamics(array_cfg, rb70, t_max_us=1.5_dp, n_steps=100, &
                                    t_arr=t_arr, z2_order_arr=z2_arr)

    n_total = n_total + 1
    ! Verify that Z2 order oscillates (minval < -0.3 and maxval > 0.3)
    if (minval(z2_arr) < -0.3_dp .and. maxval(z2_arr) > 0.3_dp) then
        print *, " [PASS] Quantum many-body scar oscillations verified: range = [", &
                 minval(z2_arr), ", ", maxval(z2_arr), "]"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Scar oscillation range unexpected: [", minval(z2_arr), ", ", maxval(z2_arr), "]"
    end if

    print *, "--------------------------------------------------"
    print *, "Rydberg Blockade Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All Rydberg blockade tests passed."
    else
        stop 1
    end if

end program test_rydberg_blockade
