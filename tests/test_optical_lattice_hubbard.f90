!> \file test_optical_lattice_hubbard.f90
!> \brief Unit tests for optical lattice Bloch bands, Bose-Hubbard parameters, and Bloch oscillations
!> \author LiHao
program test_optical_lattice_hubbard
    use mod_constants, only: dp
    use mod_optical_lattice_hubbard
    implicit none

    type(optical_lattice_t) :: latt_rb, latt_deep, latt_free
    type(bose_hubbard_param_t) :: bh_shallow, bh_deep
    real(dp), dimension(2) :: bands_free_q0, bands_free_q1
    real(dp) :: t_bloch, omega_bloch, p_lz, g_accel, force_grav
    integer  :: n_pass, n_total, stat

    n_pass = 0
    n_total = 0

    print '(A)', "=================================================="
    print '(A)', "  GeneralModule Unit Tests: Optical Lattice & BH  "
    print '(A)', "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: 87Rb Optical Lattice Recoil Scale (1064 nm)
    ! --------------------------------------------------------------------------
    call init_optical_lattice(latt_rb, mass_amu=86.909_dp, lambda_nm=1064.0_dp, s_depth=10.0_dp, stat=stat)
    n_total = n_total + 1
    ! Typical 87Rb 1064 nm recoil frequency E_R/h ~ 2.03 kHz, E_R/kB ~ 97.4 nK
    if (stat == 0 .and. latt_rb%e_recoil_hz > 1900.0_dp .and. latt_rb%e_recoil_hz < 2200.0_dp .and. &
        latt_rb%e_recoil_nkelvin > 90.0_dp .and. latt_rb%e_recoil_nkelvin < 110.0_dp) then
        print '(A, F7.2, A, F6.1, A)', " [PASS] 87Rb lattice initialized: E_R/h = ", &
            latt_rb%e_recoil_hz, " Hz, E_R/kB = ", latt_rb%e_recoil_nkelvin, " nK"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Optical lattice recoil scale initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Free Particle Analytical Band Limit (s = 0)
    ! --------------------------------------------------------------------------
    call init_optical_lattice(latt_free, mass_amu=86.909_dp, lambda_nm=1064.0_dp, s_depth=0.0_dp, stat=stat)
    call calc_bloch_band_energies(latt_free, q_quasi=0.0_dp, n_bands=2, band_energies_er=bands_free_q0)
    call calc_bloch_band_energies(latt_free, q_quasi=1.0_dp, n_bands=2, band_energies_er=bands_free_q1)

    n_total = n_total + 1
    ! For s=0: E_0(q=0) = 0, E_0(q=1) = 1.0 E_R, E_1(q=1) = 1.0 E_R (degenerate at boundary)
    if (abs(bands_free_q0(1)) < 1.0e-8_dp .and. abs(bands_free_q1(1) - 1.0_dp) < 1.0e-8_dp) then
        print '(A, F6.3, A, F6.3, A)', " [PASS] Free particle band limits: E0(q=0) = ", &
            bands_free_q0(1), " E_R, E0(q=1) = ", bands_free_q1(1), " E_R"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Free particle band limit error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Bose-Hubbard Parameters (Superfluid vs Mott Insulator)
    ! --------------------------------------------------------------------------
    call init_optical_lattice(latt_rb, mass_amu=86.909_dp, lambda_nm=1064.0_dp, s_depth=4.0_dp, stat=stat)
    call calc_bose_hubbard_parameters(latt_rb, a_s_bohr=100.0_dp, omega_perp_hz=1500.0_dp, &
                                     bh=bh_shallow, stat=stat)

    call init_optical_lattice(latt_deep, mass_amu=86.909_dp, lambda_nm=1064.0_dp, s_depth=16.0_dp, stat=stat)
    call calc_bose_hubbard_parameters(latt_deep, a_s_bohr=100.0_dp, omega_perp_hz=1500.0_dp, &
                                     bh=bh_deep, stat=stat)

    n_total = n_total + 1
    ! At s=4: J is large, U/J ~ 0.5 - 2.0 (Superfluid)
    ! At s=16: J is tiny, U/J > 3.3 (Mott Insulator candidate)
    if (bh_shallow%j_hopping_hz > bh_deep%j_hopping_hz .and. &
        bh_deep%is_mott_candidate .and. .not. bh_shallow%is_mott_candidate) then
        print '(A, F6.2, A, F6.2)', " [PASS] Superfluid vs Mott regimes: (U/J)_s=4 = ", &
            bh_shallow%u_over_j_ratio, ", (U/J)_s=16 = ", bh_deep%u_over_j_ratio
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Bose-Hubbard parameter scaling error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Band Gap Opening with Increasing Lattice Depth
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    if (bh_deep%band_gap_er > bh_shallow%band_gap_er .and. bh_deep%band_gap_er > 5.0_dp) then
        print '(A, F6.2, A, F6.2, A)', " [PASS] Band gap opening: Gap(s=4) = ", &
            bh_shallow%band_gap_er, " E_R, Gap(s=16) = ", bh_deep%band_gap_er, " E_R"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Band gap scaling error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Gravity-Driven Bloch Oscillations
    ! --------------------------------------------------------------------------
    g_accel = 9.80665_dp
    force_grav = (86.909_dp * 1.66053906660e-27_dp) * g_accel

    call calc_bloch_oscillation_dynamics(latt_deep, force_grav, t_bloch_ms=t_bloch, &
                                        omega_bloch_hz=omega_bloch, p_lz_tunnel=p_lz, stat=stat)
    n_total = n_total + 1
    ! For 87Rb in 1064 nm vertical lattice (d = 532 nm), T_B ~ 0.88 ms, P_LZ ~ 0
    if (stat == 0 .and. t_bloch > 0.70_dp .and. t_bloch < 1.10_dp .and. p_lz < 1.0e-5_dp) then
        print '(A, F6.3, A, ES10.3)', " [PASS] Gravitational Bloch oscillations: T_B = ", &
            t_bloch, " ms, P_LZ = ", p_lz
        n_pass = n_pass + 1
    else
        print '(A, F10.4, A, ES14.4)', " [FAIL] Bloch oscillation dynamics error: T_B = ", &
            t_bloch, " ms, P_LZ = ", p_lz
    end if

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', " Optical Lattice Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', " SUCCESS: All optical lattice & Bose-Hubbard tests passed."
    else
        error stop " Test failures detected in test_optical_lattice_hubbard."
    end if

end program test_optical_lattice_hubbard
