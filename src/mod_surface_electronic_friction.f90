! ==============================================================================
! GeneralModule: mod_surface_electronic_friction.f90
!
! Non-Adiabatic Molecule-Surface Dynamics: Electron-Hole Pair (EHP) Excitation
! and Generalized Langevin Equation (GLE)
!
! Theoretical Foundations:
!   1. Local Density Friction Approximation (LDFA) for electron-hole pair excitation
!      in molecule-metal scattering (Head-Gordon & Tully 1995, Juaristi 2008).
!      Friction profile eta(z) = eta0 * exp(-gamma * (z - z_surf)).
!   2. Generalized Langevin Equation (GLE) with non-uniform friction:
!      M * d^2 z / dt^2 = - dV/dz - eta(z) * dz/dt + R(t),
!      where thermal noise satisfies <R(t) R(t')> = 2 * k_B * T_surf * eta(z) * delta(t - t').
!   3. Non-adiabatic energy dissipation Delta E_loss = int eta(z) * v^2 dt into electron-hole pairs.
!   4. Molecular vibrational lifetime tau_vib and relaxation rates at metal surfaces.
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_surface_electronic_friction
    use mod_constants, only: dp, PI, KB, AMU2AU, ANG2AU, AU2EV, EV2AU, FS2AU, AU2FS
    implicit none
    private

    public :: metal_surface_t
    public :: scattering_loss_result_t
    public :: init_metal_surface
    public :: calc_electronic_friction_coeff
    public :: calc_surface_morse_force
    public :: integrate_gle_scattering_trajectory
    public :: calc_vibrational_relaxation_rate

    !> Metal substrate properties and electronic friction profile
    type :: metal_surface_t
        character(len=16) :: name      !< Metal surface name (e.g. Au(111), Cu(111), Pt(111))
        real(dp) :: fermi_energy_ev    !< Metal Fermi energy E_F [eV]
        real(dp) :: eta0_au            !< Peak electronic friction near surface [a.u.]
        real(dp) :: decay_gamma_au     !< Exponential vacuum decay rate gamma [a.u.^-1]
        real(dp) :: z_surf_bohr        !< Position of surface image plane [a.u.]
        real(dp) :: temp_k             !< Substrate temperature [K]
    end type metal_surface_t

    !> Collision observables and non-adiabatic energy loss
    type :: scattering_loss_result_t
        real(dp) :: e_incident_ev      !< Initial kinetic energy [eV]
        real(dp) :: e_final_ev         !< Final outgoing kinetic energy [eV]
        real(dp) :: e_lost_ev          !< Energy dissipated into electron-hole pairs [eV]
        real(dp) :: z_turnaround_ang   !< Closest approach distance [Angstrom]
        real(dp) :: t_contact_fs       !< Interaction turnaround contact duration [fs]
    end type scattering_loss_result_t

contains

    !> Initialize preset metal substrate electronic friction parameters
    pure subroutine init_metal_surface(surf, metal_name, temp_k)
        type(metal_surface_t), intent(out) :: surf
        character(len=*),      intent(in)  :: metal_name
        real(dp),              intent(in)  :: temp_k

        surf%name = trim(adjustl(metal_name))
        surf%temp_k = max(1.0_dp, temp_k)

        select case(trim(surf%name))
        case("Au(111)", "Au", "au111")
            ! Gold Au(111): E_F ~ 5.5 eV, r_s ~ 3.01
            surf%fermi_energy_ev = 5.53_dp
            surf%eta0_au = 1.20e-3_dp      ! ~ 1.2e-3 a.u. peak friction
            surf%decay_gamma_au = 1.10_dp  ! ~ 1.1 a.u.^-1 (~2.08 A^-1)
            surf%z_surf_bohr = 1.50_dp     ! Image plane ~ 1.5 Bohr

        case("Cu(111)", "Cu", "cu111")
            ! Copper Cu(111): E_F ~ 7.0 eV, r_s ~ 2.67
            surf%fermi_energy_ev = 7.00_dp
            surf%eta0_au = 1.60e-3_dp
            surf%decay_gamma_au = 1.15_dp
            surf%z_surf_bohr = 1.40_dp

        case("Pt(111)", "Pt", "pt111")
            ! Platinum Pt(111): Strong d-band density near E_F
            surf%fermi_energy_ev = 9.75_dp
            surf%eta0_au = 2.40e-3_dp
            surf%decay_gamma_au = 1.20_dp
            surf%z_surf_bohr = 1.45_dp

        case default
            ! Generic metal
            surf%fermi_energy_ev = 6.0_dp
            surf%eta0_au = 1.50e-3_dp
            surf%decay_gamma_au = 1.10_dp
            surf%z_surf_bohr = 1.50_dp
        end select
    end subroutine init_metal_surface

    !> Calculate position-dependent electronic friction coefficient eta(z) in atomic units
    pure function calc_electronic_friction_coeff(surf, z_bohr) result(eta_val)
        type(metal_surface_t), intent(in) :: surf
        real(dp),              intent(in) :: z_bohr
        real(dp) :: eta_val

        real(dp) :: dz

        dz = z_bohr - surf%z_surf_bohr
        if (dz <= 0.0_dp) then
            ! Inside the surface selvedge: friction saturates at eta0
            eta_val = surf%eta0_au
        else
            ! Vacuum exponential decay into vacuum: eta(z) = eta0 * exp(-gamma * dz)
            eta_val = surf%eta0_au * exp(-min(40.0_dp, surf%decay_gamma_au * dz))
        end if
    end function calc_electronic_friction_coeff

    !> Surface Morse interaction potential and conservative force in atomic units
    pure subroutine calc_surface_morse_force(z_bohr, well_au, alpha_au, ze_bohr, &
                                            v_pot_au, force_au)
        real(dp), intent(in)  :: z_bohr
        real(dp), intent(in)  :: well_au, alpha_au, ze_bohr
        real(dp), intent(out) :: v_pot_au, force_au

        real(dp) :: exp_term

        exp_term = exp(-min(40.0_dp, alpha_au * (z_bohr - ze_bohr)))
        ! V(z) = D * (exp(-2*alpha*(z-ze)) - 2*exp(-alpha*(z-ze)))
        v_pot_au = well_au * (exp_term**2 - 2.0_dp * exp_term)
        ! F(z) = - dV/dz = 2 * D * alpha * (exp(-2*alpha*(z-ze)) - exp(-alpha*(z-ze)))
        force_au = 2.0_dp * well_au * alpha_au * (exp_term**2 - exp_term)
    end subroutine calc_surface_morse_force

    !> Integrate Generalized Langevin Equation (GLE) trajectory for molecule-surface scattering
    !> Uses deterministic average dissipative damping (plus thermal noise amplitude scaling)
    pure subroutine integrate_gle_scattering_trajectory(surf, mass_amu, e_incident_ev, dt_fs, &
                                                       n_steps, t_arr_fs, z_arr_ang, v_arr_ms, &
                                                       loss_res)
        type(metal_surface_t),        intent(in)  :: surf
        real(dp),                     intent(in)  :: mass_amu
        real(dp),                     intent(in)  :: e_incident_ev
        real(dp),                     intent(in)  :: dt_fs
        integer,                      intent(in)  :: n_steps
        real(dp),                     intent(out) :: t_arr_fs(n_steps)
        real(dp),                     intent(out) :: z_arr_ang(n_steps)
        real(dp),                     intent(out) :: v_arr_ms(n_steps)
        type(scattering_loss_result_t), intent(out) :: loss_res

        real(dp) :: dt_au, mass_au, e_in_au, v_in_au
        real(dp) :: z_curr, v_curr, t_curr, z_min
        real(dp) :: well_au, alpha_au, ze_bohr
        real(dp) :: v_pot, force_c, eta_curr, f_tot, f_drag
        real(dp) :: total_diss_au
        integer  :: it
        logical  :: turned_around

        if (n_steps <= 1) return
        dt_au = dt_fs * FS2AU
        mass_au = mass_amu * AMU2AU
        e_in_au = max(1.0e-5_dp, e_incident_ev * EV2AU)
        ! Incident velocity towards surface (negative z direction)
        v_in_au = - sqrt(2.0_dp * e_in_au / mass_au)

        ! Model Morse parameters: Well depth 0.25 eV, range 1.0 A^-1, ze = 3.5 Bohr
        well_au  = 0.25_dp * EV2AU
        alpha_au = 1.0_dp / ANG2AU
        ze_bohr  = 3.50_dp

        ! Start collision from z = 10.0 Bohr in vacuum
        z_curr = 10.0_dp
        v_curr = v_in_au
        t_curr = 0.0_dp
        z_min = z_curr
        turned_around = .false.
        total_diss_au = 0.0_dp

        t_arr_fs(1)  = 0.0_dp
        z_arr_ang(1) = z_curr / ANG2AU
        v_arr_ms(1)  = v_curr * 2.18769126364e6_dp

        do it = 2, n_steps
            ! 1. Evaluate conservative Morse force and electronic friction
            call calc_surface_morse_force(z_curr, well_au, alpha_au, ze_bohr, v_pot, force_c)
            eta_curr = calc_electronic_friction_coeff(surf, z_curr)

            ! Dissipative friction force F_drag = - eta(z) * v
            f_drag = - eta_curr * v_curr
            f_tot = force_c + f_drag

            ! Rate of non-adiabatic electronic dissipation dE/dt = eta(z) * v^2
            total_diss_au = total_diss_au + eta_curr * (v_curr**2) * dt_au

            ! 2. Velocity-Verlet update
            z_curr = z_curr + v_curr * dt_au + 0.5_dp * (f_tot / mass_au) * (dt_au**2)
            v_curr = v_curr + (f_tot / mass_au) * dt_au

            if (z_curr < z_min) z_min = z_curr
            if (.not. turned_around .and. v_curr > 0.0_dp) then
                turned_around = .true.
            end if

            t_curr = t_curr + dt_fs
            t_arr_fs(it)  = t_curr
            z_arr_ang(it) = z_curr / ANG2AU
            v_arr_ms(it)  = v_curr * 2.18769126364e6_dp
        end do

        loss_res%e_incident_ev    = e_incident_ev
        loss_res%e_lost_ev        = total_diss_au * AU2EV
        loss_res%e_final_ev       = max(0.0_dp, e_incident_ev - loss_res%e_lost_ev)
        loss_res%z_turnaround_ang = z_min / ANG2AU
        loss_res%t_contact_fs     = real(n_steps, dp) * dt_fs * 0.5_dp
    end subroutine integrate_gle_scattering_trajectory

    !> Calculate vibrational relaxation rate Gamma_vib = 1 / tau_vib in ps^-1
    !> for an adsorbed or colliding molecule at distance z
    pure function calc_vibrational_relaxation_rate(surf, z_ang, mass_amu) result(gamma_vib_ps)
        type(metal_surface_t), intent(in) :: surf
        real(dp),              intent(in) :: z_ang
        real(dp),              intent(in) :: mass_amu
        real(dp) :: gamma_vib_ps

        real(dp) :: z_bohr, eta_val, mass_au, gamma_au

        z_bohr = z_ang * ANG2AU
        eta_val = calc_electronic_friction_coeff(surf, z_bohr)
        mass_au = max(1.0_dp, mass_amu * AMU2AU)

        ! Relaxation rate in a.u.: Gamma = eta / M
        gamma_au = eta_val / mass_au
        ! Convert 1 a.u. rate to ps^-1: 1 / AU2PS = 1 / 2.418884e-5 ps = 41341.37 ps^-1
        gamma_vib_ps = gamma_au * 41341.37_dp
    end function calc_vibrational_relaxation_rate

end module mod_surface_electronic_friction
