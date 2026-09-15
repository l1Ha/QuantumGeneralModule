! ==============================================================================
! GeneralModule: test_surface_electronic_friction.f90
!
! Unit Test Suite: Surface Electronic Friction & Non-Adiabatic EHP Dissipation
!
! Standard: Fortran 2008
! ==============================================================================

program test_surface_electronic_friction
    use mod_constants, only: dp
    use mod_surface_electronic_friction
    implicit none

    integer :: n_pass, n_total
    type(metal_surface_t) :: au_surf
    type(scattering_loss_result_t) :: loss
    real(dp) :: eta_near, eta_far
    real(dp) :: v_pot, force_eq, force_rep
    real(dp) :: t_arr(200), z_arr(200), v_arr(200)
    real(dp) :: gamma_vib

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Electronic Friction   "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Metal Surface Preset Initialization
    ! --------------------------------------------------------------------------
    call init_metal_surface(au_surf, "Au(111)", temp_k=300.0_dp)

    n_total = n_total + 1
    if (abs(au_surf%fermi_energy_ev - 5.53_dp) < 0.01_dp .and. &
        au_surf%eta0_au > 1.0e-3_dp .and. au_surf%decay_gamma_au > 1.0_dp) then
        print *, " [PASS] Au(111) initialized: E_F = ", au_surf%fermi_energy_ev, &
                 " eV, Peak friction eta0 = ", au_surf%eta0_au, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Metal surface initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Position-Dependent Electronic Friction Profile
    ! --------------------------------------------------------------------------
    ! z = 1.0 Bohr (inside surface selvedge) vs z = 8.0 Bohr (vacuum)
    eta_near = calc_electronic_friction_coeff(au_surf, z_bohr=1.0_dp)
    eta_far  = calc_electronic_friction_coeff(au_surf, z_bohr=8.0_dp)

    n_total = n_total + 1
    if (abs(eta_near - au_surf%eta0_au) < 1.0e-12_dp .and. &
        eta_far < 1.0e-5_dp .and. eta_far >= 0.0_dp) then
        print *, " [PASS] Friction profile verified: eta(z=1 Bohr) = ", eta_near, &
                 " a.u., eta(z=8 Bohr) = ", eta_far, " a.u. (decayed!)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Electronic friction profile error: near = ", eta_near, " far = ", eta_far
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Surface Morse Force & Equilibrium Check
    ! --------------------------------------------------------------------------
    ! Equilibrium distance ze = 3.5 Bohr
    call calc_surface_morse_force(z_bohr=3.5_dp, well_au=0.01_dp, alpha_au=0.5_dp, &
                                  ze_bohr=3.5_dp, v_pot_au=v_pot, force_au=force_eq)
    ! Repulsive wall at z = 2.0 Bohr < ze
    call calc_surface_morse_force(z_bohr=2.0_dp, well_au=0.01_dp, alpha_au=0.5_dp, &
                                  ze_bohr=3.5_dp, v_pot_au=v_pot, force_au=force_rep)

    n_total = n_total + 1
    if (abs(force_eq) < 1.0e-12_dp .and. force_rep > 0.0_dp) then
        print *, " [PASS] Conservative Morse force verified: F(ze) = ", force_eq, &
                 ", F(repulsive wall) = ", force_rep, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Morse force error: force_eq = ", force_eq, " force_rep = ", force_rep
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: GLE Trajectory & EHP Non-Adiabatic Energy Loss
    ! --------------------------------------------------------------------------
    ! NO molecule (M = 30 amu) colliding with Au(111) at E_i = 0.50 eV
    call integrate_gle_scattering_trajectory(au_surf, mass_amu=30.0_dp, e_incident_ev=0.50_dp, &
                                            dt_fs=0.5_dp, n_steps=200, t_arr_fs=t_arr, &
                                            z_arr_ang=z_arr, v_arr_ms=v_arr, loss_res=loss)

    n_total = n_total + 1
    ! EHP energy loss must be positive, less than total incident energy, and turnaround must happen
    if (loss%e_lost_ev > 1.0e-8_dp .and. loss%e_lost_ev < 0.50_dp .and. &
        loss%z_turnaround_ang > 0.5_dp .and. loss%e_final_ev < loss%e_incident_ev) then
        print *, " [PASS] GLE trajectory & EHP dissipation verified:"
        print *, "        E_incident = ", loss%e_incident_ev, " eV"
        print *, "        E_lost     = ", loss%e_lost_ev, " eV (dissipated to electron-hole pairs)"
        print *, "        E_final    = ", loss%e_final_ev, " eV, Turnaround = ", loss%z_turnaround_ang, " A"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] GLE trajectory error: E_lost = ", loss%e_lost_ev, " E_final = ", loss%e_final_ev
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Molecular Vibrational Relaxation Rate
    ! --------------------------------------------------------------------------
    gamma_vib = calc_vibrational_relaxation_rate(au_surf, z_ang=1.0_dp, mass_amu=30.0_dp)

    n_total = n_total + 1
    ! Vibrational quenching rate on metal surface should be finite (~ 0.01 - 10 ps^-1)
    if (gamma_vib > 1.0e-4_dp .and. gamma_vib < 100.0_dp) then
        print *, " [PASS] Molecular vibrational relaxation rate verified: Gamma_vib = ", gamma_vib, " ps^-1"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Vibrational relaxation rate error: ", gamma_vib
    end if

    print *, "--------------------------------------------------"
    print *, "Electronic Friction Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All surface electronic friction tests passed."
    else
        stop 1
    end if

end program test_surface_electronic_friction
