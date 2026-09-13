! ==============================================================================
! GeneralModule: test_dipolar_droplets_lhy.f90
!
! Unit Test Suite: Dipolar Quantum Droplets & Lee-Huang-Yang (LHY) Corrections
!
! Standard: Fortran 2008
! ==============================================================================

program test_dipolar_droplets_lhy
    use mod_constants, only: dp
    use mod_dipolar_droplets_lhy
    implicit none

    integer :: n_pass, n_total, stat
    type(dipolar_droplet_param_t) :: dy_droplet, dy_gas
    real(dp) :: q5_zero, q5_val, n0_au, n0_cm3, mu_eq, n_crit, e_dens

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Dipolar Droplets & LHY "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Dysprosium-162 Dipole Length a_dd
    ! --------------------------------------------------------------------------
    call init_dipolar_droplet_param("162Dy", a_scat_bohr=70.0_dp, param=dy_droplet, stat=stat)

    n_total = n_total + 1
    ! Experimental a_dd for 162Dy is approximately 130 - 132 Bohr
    if (stat == 0 .and. dy_droplet%a_dd_au > 120.0_dp .and. dy_droplet%a_dd_au < 140.0_dp) then
        print *, " [PASS] 162Dy dipole length a_dd = ", dy_droplet%a_dd_au, " a0 (Exp: ~131 a0)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected dipole length for Dy: ", dy_droplet%a_dd_au
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Droplet Regime Condition (epsilon_dd > 1.0)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    if (dy_droplet%is_droplet_regime .and. dy_droplet%epsilon_dd > 1.0_dp) then
        print *, " [PASS] Droplet regime active for a_s = 70 a0: epsilon_dd = ", dy_droplet%epsilon_dd
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Droplet regime not identified for a_s = 70 a0"
    end if

    ! For large scattering length a_s = 250 a0 > a_dd, droplet cannot form (gas phase)
    call init_dipolar_droplet_param("162Dy", a_scat_bohr=250.0_dp, param=dy_gas, stat=stat)
    n_total = n_total + 1
    if (.not. dy_gas%is_droplet_regime .and. dy_gas%epsilon_dd < 1.0_dp) then
        print *, " [PASS] Dilute gas regime for a_s = 250 a0: epsilon_dd = ", dy_gas%epsilon_dd, " < 1"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Gas regime detection error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Pelster-Lima Auxiliary Function Q_5(epsilon_dd)
    ! --------------------------------------------------------------------------
    q5_zero = calc_pelster_lima_q5(0.0_dp)
    n_total = n_total + 1
    if (abs(q5_zero - 1.0_dp) < 1.0e-12_dp) then
        print *, " [PASS] Non-dipolar limit Q_5(0) = 1.00000 (Exact LHY)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Q_5(0) mismatch: ", q5_zero
    end if

    q5_val = calc_pelster_lima_q5(dy_droplet%epsilon_dd)
    n_total = n_total + 1
    if (q5_val > 1.0_dp) then
        print *, " [PASS] Dipolar quantum fluctuation enhancement: Q_5(eps) = ", q5_val, " > 1.0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Q_5 dipolar enhancement missing"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Equilibrium Self-Bound Droplet Density n_0
    ! --------------------------------------------------------------------------
    n0_au = calc_equilibrium_droplet_density(dy_droplet)
    n0_cm3 = n0_au / ((5.29177210903e-9_dp)**3)

    n_total = n_total + 1
    ! Typical dipolar droplet density is ~ 10^14 - 10^17 cm^-3
    if (n0_cm3 > 1.0e13_dp .and. n0_cm3 < 1.0e18_dp) then
        print *, " [PASS] Self-bound droplet equilibrium density n_0 = ", n0_cm3, " cm^-3"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Droplet equilibrium density out of physical range: ", n0_cm3
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Negative Chemical Potential mu_0 < 0 (Self-Bound Criterion)
    ! --------------------------------------------------------------------------
    mu_eq = calc_droplet_chemical_potential(dy_droplet, n0_au)

    n_total = n_total + 1
    if (mu_eq < 0.0_dp) then
        print *, " [PASS] Self-bound negative chemical potential: mu_0 = ", mu_eq, " a.u. < 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Non-negative chemical potential: ", mu_eq
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: Critical Atom Number N_crit & Energy Density
    ! --------------------------------------------------------------------------
    n_crit = calc_critical_atom_number(dy_droplet)
    n_total = n_total + 1
    if (n_crit > 10.0_dp .and. n_crit < 1.0e6_dp) then
        print *, " [PASS] Critical droplet atom number N_crit = ", n_crit
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Critical atom number out of bounds: ", n_crit
    end if

    e_dens = calc_egpe_energy_density(dy_droplet, n0_au)
    n_total = n_total + 1
    if (e_dens < 0.0_dp) then
        print *, " [PASS] Droplet bound state energy density E/V = ", e_dens, " a.u. < 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Energy density not negative: ", e_dens
    end if

    print *, "--------------------------------------------------"
    print *, "Dipolar Droplets & LHY Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All dipolar quantum droplet tests passed."
    else
        stop 1
    end if

end program test_dipolar_droplets_lhy
