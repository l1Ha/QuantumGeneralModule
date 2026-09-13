! ==============================================================================
! GeneralModule Example 18: Ultracold Dipolar Quantum Droplets & LHY Fluctuations
!
! Features:
!  1. Dysprosium-162 (162Dy) magnetic dipole length a_dd ~ 130 a0 and dipole strength epsilon_dd
!  2. Beyond-mean-field Pelster-Lima auxiliary function Q_5(epsilon_dd)
!  3. Self-bound flat-top equilibrium density n_0 across the droplet-to-gas transition
!  4. Extended Gross-Pitaevskii (eGPE) energy density E/V and negative chemical potential mu < 0
!  5. Petrov-Chomaz critical atom number N_crit for 3D self-binding stability
!
! Standard: Fortran 2008
! ==============================================================================
program ex18_dipolar_quantum_droplets
    use mod_constants, only: dp, PI
    use mod_dipolar_droplets_lhy
    implicit none

    type(dipolar_droplet_param_t) :: param
    real(dp) :: a_s, eps_dd, q5, n0_au, n0_cm3, mu_0, e_dens, n_crit
    integer  :: u_out, i, n_scan, stat

    print *, "================================================================"
    print *, " Example 18: Dipolar Quantum Droplets & Lee-Huang-Yang Physics  "
    print *, "================================================================"

    ! 1. Initialize reference 162Dy droplet at a_s = 70 a0
    call init_dipolar_droplet_param("162Dy", 70.0_dp, param, stat)

    write(*, '(A, F10.2, A)') " 162Dy magnetic dipole length a_dd    : ", param%a_dd_au, " a0"
    write(*, '(A, F10.4)')    " Relative dipole strength epsilon_dd  : ", param%epsilon_dd
    write(*, '(A, F10.4)')    " Pelster-Lima LHY factor Q_5(eps_dd)  : ", param%q5_factor
    n0_au = calc_equilibrium_droplet_density(param)
    n0_cm3 = n0_au / ((5.29177210903e-9_dp)**3)
    mu_0   = calc_droplet_chemical_potential(param, n0_au)
    n_crit = calc_critical_atom_number(param)
    write(*, '(A, ES12.4, A)') " Self-bound equilibrium density n_0   : ", n0_cm3, " cm^-3"
    write(*, '(A, ES12.4, A)') " Bound chemical potential mu(n_0)     : ", mu_0, " a.u. (< 0)"
    write(*, '(A, F10.1)')     " Critical atom number N_crit          : ", n_crit
    print *, "----------------------------------------------------------------"

    ! 2. Scan scattering length across transition: a_s = 60 to 140 a0
    open(newunit=u_out, file="ex18_dipolar_droplet.dat", status="replace", action="write")
    write(u_out, '(A)') "# GeneralModule Example 18: 162Dy Dipolar Quantum Droplets"
    write(u_out, '(A)') "# Col 1: Scattering length a_s (Bohr)"
    write(u_out, '(A)') "# Col 2: Relative dipole strength epsilon_dd = a_dd / a_s"
    write(u_out, '(A)') "# Col 3: Pelster-Lima factor Q_5(epsilon_dd)"
    write(u_out, '(A)') "# Col 4: Equilibrium density n_0 (cm^-3)"
    write(u_out, '(A)') "# Col 5: Chemical potential mu(n_0) (a.u.)"
    write(u_out, '(A)') "# Col 6: Energy density E/V (a.u.)"
    write(u_out, '(A)') "# Col 7: Critical atom number N_crit"

    n_scan = 41
    write(*, '(A)') "   a_s (a0) | epsilon_dd |   Q_5   |    n_0 (cm^-3)    |   mu_0 (a.u.)   |   N_crit"
    write(*, '(A)') "------------+------------+---------+-------------------+-----------------+----------"

    do i = 1, n_scan
        a_s = 60.0_dp + real(i - 1, dp) * 2.0_dp
        call init_dipolar_droplet_param("162Dy", a_s, param, stat)

        eps_dd = param%epsilon_dd
        q5     = param%q5_factor

        if (param%is_droplet_regime) then
            n0_au  = calc_equilibrium_droplet_density(param)
            n0_cm3 = n0_au / ((5.29177210903e-9_dp)**3)
            mu_0   = calc_droplet_chemical_potential(param, n0_au)
            e_dens = calc_egpe_energy_density(param, n0_au)
            n_crit = calc_critical_atom_number(param)
        else
            n0_au  = 0.0_dp
            n0_cm3 = 0.0_dp
            mu_0   = 0.0_dp
            e_dens = 0.0_dp
            n_crit = 1.0e10_dp
        end if

        if (mod(i, 5) == 1) then
            write(*, '(F11.1, " | ", F10.3, " | ", F7.3, " | ", ES17.4, " | ", ES15.4, " | ", F9.1)') &
                a_s, eps_dd, q5, n0_cm3, mu_0, min(n_crit, 99999.9_dp)
        end if

        write(u_out, '(7ES16.6)') a_s, eps_dd, q5, n0_cm3, mu_0, e_dens, n_crit
    end do

    close(u_out)
    print *, "----------------------------------------------------------------"
    print *, " Output data saved to ex18_dipolar_droplet.dat"
    print *, " Example 18 completed successfully."

end program ex18_dipolar_quantum_droplets
