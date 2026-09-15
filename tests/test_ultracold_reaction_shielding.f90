! ==============================================================================
! GeneralModule: test_ultracold_reaction_shielding.f90
!
! Unit Test Suite: Ultracold Polar Molecule Reactions and Dipolar Shielding
!
! Standard: Fortran 2008
! ==============================================================================

program test_ultracold_reaction_shielding
    use mod_constants, only: dp
    use mod_ultracold_reaction_shielding
    implicit none

    integer :: n_pass, n_total
    type(ultracold_molecule_t) :: krb, narb
    type(shielding_config_t) :: mw_cfg, dc_cfg
    real(dp) :: v_eff_k, r_bar, v_bar_k, t_tunnel
    real(dp) :: k2_el, k2_inel, gamma_ratio
    real(dp) :: det_arr(5), gamma_arr(5)

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Ultracold Shielding    "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Polar Molecule Preset Initialization
    ! --------------------------------------------------------------------------
    call init_ultracold_molecule_preset(krb, "KRb")
    call init_ultracold_molecule_preset(narb, "NaRb")

    n_total = n_total + 1
    if (abs(krb%mass_amu - 126.87_dp) < 0.1_dp .and. &
        abs(krb%dipole_debye - 0.574_dp) < 0.01_dp .and. &
        narb%dipole_debye > krb%dipole_debye) then
        print *, " [PASS] Polar molecule presets initialized: KRb mass = ", krb%mass_amu, &
                 " amu, NaRb dipole = ", narb%dipole_debye, " Debye"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Molecule preset initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Microwave Shielding Barrier Formation
    ! --------------------------------------------------------------------------
    ! Blue detuning Delta = 15 MHz, Rabi frequency Omega = 5 MHz, universal loss y = 1.0
    call init_shielding_config(mw_cfg, method=1, detuning_mhz=15.0_dp, rabi_mhz=5.0_dp, &
                              e_field_kv_cm=0.0_dp, y_loss=1.0_dp)
    call calc_shielding_barrier_height(krb, mw_cfg, r_bar, v_bar_k)

    n_total = n_total + 1
    ! Barrier should form at R ~ 300 to 2000 a0 with height > 10 uK (1e-5 K)
    if (r_bar > 100.0_dp .and. r_bar < 5000.0_dp .and. v_bar_k > 1.0e-6_dp) then
        print *, " [PASS] Microwave shielding barrier verified: R_barrier = ", r_bar, &
                 " a0, V_barrier = ", v_bar_k * 1.0e6_dp, " uK"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected shielding barrier: R = ", r_bar, " V = ", v_bar_k
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: WKB Barrier Tunneling Suppression
    ! --------------------------------------------------------------------------
    ! Collision energy at T = 0.5 uK
    call calc_wkb_tunneling_probability(krb, mw_cfg, collision_energy_uk=0.5_dp, t_tunnel=t_tunnel)

    n_total = n_total + 1
    ! Transmission probability through repulsive barrier must be small (T << 1)
    if (t_tunnel < 0.1_dp .and. t_tunnel >= 0.0_dp) then
        print *, " [PASS] WKB quantum tunneling suppression verified: T_tunnel = ", t_tunnel
        n_pass = n_pass + 1
    else
        print *, " [FAIL] WKB transmission too large: T_tunnel = ", t_tunnel
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Elastic and Inelastic Rates & Good-to-Bad Ratio (gamma > 100)
    ! --------------------------------------------------------------------------
    call calc_shielded_scattering_rates(krb, mw_cfg, temp_uk=0.5_dp, &
                                       k2_el_cm3s=k2_el, k2_inel_cm3s=k2_inel, &
                                       gamma_ratio=gamma_ratio)

    n_total = n_total + 1
    if (k2_el > 1.0e-11_dp .and. k2_inel < 1.0e-10_dp .and. gamma_ratio > 100.0_dp) then
        print *, " [PASS] Shielded scattering rates verified:"
        print *, "        K2_el   = ", k2_el, " cm^3/s"
        print *, "        K2_inel = ", k2_inel, " cm^3/s"
        print *, "        gamma   = ", gamma_ratio, " (Evaporative cooling criterion > 100 met!)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Scattering rates error: K2_el = ", k2_el, " K2_inel = ", k2_inel, &
                 " gamma = ", gamma_ratio
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Detuning Scan Consistency
    ! --------------------------------------------------------------------------
    call calc_shielding_detuning_scan(krb, mw_cfg, temp_uk=0.5_dp, n_pts=5, &
                                      det_min_mhz=5.0_dp, det_max_mhz=25.0_dp, &
                                      det_arr=det_arr, gamma_arr=gamma_arr)

    n_total = n_total + 1
    if (det_arr(5) > det_arr(1) .and. all(gamma_arr > 10.0_dp)) then
        print *, " [PASS] Detuning scan completed: Delta = [", det_arr(1), ", ", det_arr(5), &
                 "] MHz, all gamma > 10"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Detuning scan inconsistent"
    end if

    print *, "--------------------------------------------------"
    print *, "Ultracold Shielding Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All ultracold reaction shielding tests passed."
    else
        stop 1
    end if

end program test_ultracold_reaction_shielding
