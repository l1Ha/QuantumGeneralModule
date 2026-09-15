! ==============================================================================
! GeneralModule: test_surface_reaction_er.f90
!
! Unit Test Suite: Surface Reactions & Eley-Rideal Catalytic Abstraction
!
! Standard: Fortran 2008
! ==============================================================================

program test_surface_reaction_er
    use mod_constants, only: dp
    use mod_surface_reaction_er
    implicit none

    integer :: n_pass, n_total
    type(er_reaction_system_t)  :: h_cu
    type(er_energy_partition_t) :: part
    real(dp) :: v_reac, v_prod
    real(dp) :: v_dist(0:8), sum_pv
    real(dp) :: sigma_er, k_er

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Surface ER Reaction    "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Preset System Initialization & Exothermicity
    ! --------------------------------------------------------------------------
    call init_er_reaction_system(h_cu, "H+H/Cu(111)")

    n_total = n_total + 1
    if (abs(h_cu%d_mol_ev - 4.75_dp) < 0.01_dp .and. &
        abs(h_cu%delta_e_exo_ev - 2.30_dp) < 0.01_dp .and. &
        h_cu%delta_e_exo_ev > 0.0_dp) then
        print *, " [PASS] H+H/Cu(111) initialized: D_mol = ", h_cu%d_mol_ev, &
                 " eV, D_chem = ", h_cu%d_chem_ev, " eV, Exothermicity Delta E = ", h_cu%delta_e_exo_ev, " eV"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] ER preset initialization error: Delta E = ", h_cu%delta_e_exo_ev
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: 2D Reactive Potential Energy Surface
    ! --------------------------------------------------------------------------
    ! Asymptote 1: Reactant channel (separated gas atom r=2.5 A, adatom at well z_ad ~ 1.0 A)
    call calc_er_potential_2d(h_cu, r_bond_ang=2.5_dp, z_cm_ang=2.0_dp, v_pot_ev=v_reac)
    ! Asymptote 2: Desorbed product channel (H2 bound r=0.741 A, desorbed far Z_cm=5.0 A)
    call calc_er_potential_2d(h_cu, r_bond_ang=0.741_dp, z_cm_ang=5.0_dp, v_pot_ev=v_prod)

    n_total = n_total + 1
    ! Desorbed H2 channel should be lower in energy by approx the exothermicity (-2.3 eV)
    if (v_prod < v_reac .and. v_prod < -1.5_dp) then
        print *, " [PASS] 2D Reactive PES verified: V_react = ", v_reac, " eV, V_prod = ", v_prod, " eV"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] 2D PES calculation error: V_reac = ", v_reac, " V_prod = ", v_prod
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Exothermicity Energy Budget Partitioning
    ! --------------------------------------------------------------------------
    ! Projectile incident energy E_i = 0.1 eV
    call calc_er_energy_partitioning(h_cu, e_incident_ev=0.10_dp, partition=part)

    n_total = n_total + 1
    ! Vibrational energy should receive ~50% (part%f_vib ~ 0.50), sum of parts = total
    if (abs(part%e_total_avail_ev - 2.40_dp) < 0.01_dp .and. &
        part%e_vib_ev > part%e_trans_ev .and. &
        abs(part%e_vib_ev + part%e_trans_ev + part%e_diss_ev + part%e_rot_ev - part%e_total_avail_ev) < 1.0e-12_dp) then
        print *, " [PASS] Energy partitioning verified: E_avail = ", part%e_total_avail_ev, " eV"
        print *, "        E_vib = ", part%e_vib_ev, " eV (~50%), E_trans = ", part%e_trans_ev, " eV (~35%)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Energy partitioning error: E_avail = ", part%e_total_avail_ev
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Nascent Vibrational Inversion P(v=2) > P(v=0)
    ! --------------------------------------------------------------------------
    call calc_er_vibrational_populations(h_cu, e_incident_ev=0.10_dp, max_v=8, v_dist=v_dist)
    sum_pv = sum(v_dist)

    n_total = n_total + 1
    ! Eley-Rideal signature: vibrational inversion! P(v=2) significantly larger than ground state P(v=0)
    if (abs(sum_pv - 1.0_dp) < 1.0e-4_dp .and. v_dist(2) > v_dist(0)) then
        print *, " [PASS] Vibrational population inversion confirmed:"
        print *, "        P(v=0) = ", v_dist(0), ", P(v=1) = ", v_dist(1), ", P(v=2) = ", v_dist(2), " (Inverted peak!)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Vibrational distribution unexpected: sum = ", sum_pv, " P(0)=", v_dist(0), " P(2)=", v_dist(2)
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Cross Section & Thermal Rate Constant
    ! --------------------------------------------------------------------------
    sigma_er = calc_er_reaction_cross_section(h_cu, e_incident_ev=0.15_dp)
    k_er = calc_er_thermal_rate_constant(h_cu, temp_k=300.0_dp)

    n_total = n_total + 1
    if (sigma_er > 0.1_dp .and. sigma_er < 1.0_dp .and. k_er > 1.0e-14_dp) then
        print *, " [PASS] ER cross section & rate constant verified:"
        print *, "        sigma_ER(E=0.15 eV) = ", sigma_er, " Angstrom^2"
        print *, "        k_ER(T=300 K)       = ", k_er, " cm^3/s"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Reaction cross section error: sigma = ", sigma_er, " k = ", k_er
    end if

    print *, "--------------------------------------------------"
    print *, "Surface ER Reaction Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All surface ER reaction tests passed."
    else
        stop 1
    end if

end program test_surface_reaction_er
