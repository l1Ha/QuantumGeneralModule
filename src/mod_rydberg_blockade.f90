! ==============================================================================
! GeneralModule: mod_rydberg_blockade.f90
!
! Rydberg Atom Blockade and Many-Body Quantum Dynamics
!
! Theoretical Foundations:
!   1. Huge dipole moments and van der Waals interactions:
!      C6 proportional to n^11, blockade radius R_b = (|C6| / (hbar*Omega))^(1/6).
!   2. Two-atom Rydberg blockade and collective Rabi frequency sqrt(2)*Omega:
!      Doubly-excited state |rr> is shifted out of resonance by V_vdW = C6 / R^6.
!   3. 1D Rydberg atom array & PXP constrained dynamics:
!      Non-thermalizing quantum many-body scars with persistent Z2 oscillations.
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_rydberg_blockade
    use mod_constants, only: dp, PI, TWOPI
    implicit none
    private

    public :: rydberg_atom_t
    public :: rydberg_array_config_t
    public :: init_rydberg_atom
    public :: calc_rydberg_blockade_radius
    public :: init_rydberg_array
    public :: calc_two_atom_dynamics
    public :: calc_z2_order_parameter
    public :: calc_rydberg_scar_dynamics

    !> Rydberg atom parameters
    type :: rydberg_atom_t
        character(len=16) :: species   !< Atomic species (e.g. 87Rb, 133Cs)
        integer  :: n_principal        !< Principal quantum number n (e.g. 50-100)
        integer  :: l_orbital          !< Orbital quantum number (0 for s, 1 for p, 2 for d)
        real(dp) :: c6_mhz_um6         !< van der Waals coefficient C6 / h in MHz*um^6
        real(dp) :: lifetime_us        !< Radiative lifetime in microseconds
    end type rydberg_atom_t

    !> 1D Rydberg atom array configuration
    type :: rydberg_array_config_t
        integer  :: n_atoms            !< Number of atoms in array
        real(dp) :: spacing_um         !< Interatomic distance a [um]
        real(dp) :: rabi_mhz           !< Rabi frequency Omega / (2*pi) [MHz]
        real(dp) :: detuning_mhz       !< Detuning Delta / (2*pi) [MHz]
        integer  :: boundary_cond      !< 1: Periodic (PBC), 2: Open (OBC)
    end type rydberg_array_config_t

contains

    !> Initialize Rydberg atom with analytical n^11 scaling
    pure subroutine init_rydberg_atom(atom, species, n_principal, l_orbital)
        type(rydberg_atom_t), intent(out) :: atom
        character(len=*), intent(in) :: species
        integer, intent(in) :: n_principal
        integer, intent(in) :: l_orbital

        real(dp) :: n_eff, c6_ref

        atom%species = trim(adjustl(species))
        atom%n_principal = max(10, n_principal)
        atom%l_orbital = l_orbital

        ! Quantum defect for 87Rb: delta_s ~ 3.13, delta_p ~ 2.65, delta_d ~ 1.34
        if (atom%l_orbital == 0) then
            n_eff = real(atom%n_principal, dp) - 3.131_dp
        else if (atom%l_orbital == 1) then
            n_eff = real(atom%n_principal, dp) - 2.654_dp
        else
            n_eff = real(atom%n_principal, dp) - 1.348_dp
        end if

        ! For 87Rb 70S: n_eff ~ 66.87, C6 / h ~ 8.62e5 MHz*um^6
        ! Scaling C6 ~ C6_ref * (n_eff / 66.869)^11
        c6_ref = 8.62e5_dp  ! MHz * um^6
        atom%c6_mhz_um6 = c6_ref * (n_eff / 66.869_dp)**11

        ! Lifetime scaling tau ~ tau_ref * (n_eff / 66.869)^3 (approx 150 us for 70S at 300K)
        atom%lifetime_us = 150.0_dp * (n_eff / 66.869_dp)**3
    end subroutine init_rydberg_atom

    !> Calculate Rydberg blockade radius R_b = (|C6| / (h*Omega))^(1/6) in micrometers
    pure function calc_rydberg_blockade_radius(atom, rabi_mhz) result(r_b)
        type(rydberg_atom_t), intent(in) :: atom
        real(dp), intent(in) :: rabi_mhz
        real(dp) :: r_b

        real(dp) :: omega_clamped

        omega_clamped = max(1.0e-6_dp, abs(rabi_mhz))
        r_b = (abs(atom%c6_mhz_um6) / omega_clamped)**(1.0_dp / 6.0_dp)
    end function calc_rydberg_blockade_radius

    !> Initialize 1D array configuration
    pure subroutine init_rydberg_array(cfg, n_atoms, spacing_um, rabi_mhz, &
                                       detuning_mhz, boundary_cond)
        type(rydberg_array_config_t), intent(out) :: cfg
        integer,  intent(in) :: n_atoms
        real(dp), intent(in) :: spacing_um
        real(dp), intent(in) :: rabi_mhz
        real(dp), intent(in) :: detuning_mhz
        integer,  intent(in) :: boundary_cond

        cfg%n_atoms = max(2, n_atoms)
        cfg%spacing_um = max(0.1_dp, spacing_um)
        cfg%rabi_mhz = rabi_mhz
        cfg%detuning_mhz = detuning_mhz
        cfg%boundary_cond = boundary_cond
    end subroutine init_rydberg_array

    !> Compute exact 2-atom Rydberg dynamics using 4th-order Runge-Kutta
    !> Basis: 1: |gg>, 2: |gr>, 3: |rg>, 4: |rr>
    pure subroutine calc_two_atom_dynamics(atom, spacing_um, rabi_mhz, detuning_mhz, &
                                           t_max_us, n_steps, t_arr, p_g, p_single, p_double)
        type(rydberg_atom_t), intent(in) :: atom
        real(dp), intent(in)  :: spacing_um
        real(dp), intent(in)  :: rabi_mhz
        real(dp), intent(in)  :: detuning_mhz
        real(dp), intent(in)  :: t_max_us
        integer,  intent(in)  :: n_steps
        real(dp), intent(out) :: t_arr(n_steps)
        real(dp), intent(out) :: p_g(n_steps)
        real(dp), intent(out) :: p_single(n_steps)
        real(dp), intent(out) :: p_double(n_steps)

        real(dp) :: omega_rad, delta_rad, v_rr_rad
        real(dp) :: dt, t_curr, dt_sub, max_rate
        complex(dp) :: psi(4), k1(4), k2(4), k3(4), k4(4), psi_tmp(4)
        integer :: it, n_sub, isub

        if (n_steps <= 1) return
        dt = t_max_us / real(n_steps - 1, dp)

        ! Angular frequencies in rad/us (since 1 MHz = 2*pi rad/us)
        omega_rad = TWOPI * rabi_mhz
        delta_rad = TWOPI * detuning_mhz
        v_rr_rad  = TWOPI * (atom%c6_mhz_um6 / (spacing_um**6))

        ! Sub-stepping to satisfy RK4 numerical stability for stiff interaction V_rr
        max_rate = max(1.0_dp, v_rr_rad + omega_rad + abs(delta_rad))
        n_sub = max(1, int(dt * max_rate / 0.1_dp) + 1)
        dt_sub = dt / real(n_sub, dp)

        ! Initial state: |gg> = (1, 0, 0, 0)^T
        psi = [(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]

        t_curr = 0.0_dp
        t_arr(1) = 0.0_dp
        p_g(1) = real(psi(1) * conjg(psi(1)), dp)
        p_single(1) = real(psi(2) * conjg(psi(2)) + psi(3) * conjg(psi(3)), dp)
        p_double(1) = real(psi(4) * conjg(psi(4)), dp)

        do it = 2, n_steps
            do isub = 1, n_sub
                ! RK4 integration of d psi / dt = -i H psi
                call rydberg_deriv(psi, omega_rad, delta_rad, v_rr_rad, k1)

                psi_tmp = psi + 0.5_dp * dt_sub * k1
                call rydberg_deriv(psi_tmp, omega_rad, delta_rad, v_rr_rad, k2)

                psi_tmp = psi + 0.5_dp * dt_sub * k2
                call rydberg_deriv(psi_tmp, omega_rad, delta_rad, v_rr_rad, k3)

                psi_tmp = psi + dt_sub * k3
                call rydberg_deriv(psi_tmp, omega_rad, delta_rad, v_rr_rad, k4)

                psi = psi + (dt_sub / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
            end do

            ! Re-normalize to prevent numerical drift
            psi = psi / sqrt(sum(real(psi * conjg(psi), dp)))

            t_curr = t_curr + dt
            t_arr(it) = t_curr
            p_g(it) = real(psi(1) * conjg(psi(1)), dp)
            p_single(it) = real(psi(2) * conjg(psi(2)) + psi(3) * conjg(psi(3)), dp)
            p_double(it) = real(psi(4) * conjg(psi(4)), dp)
        end do
    end subroutine calc_two_atom_dynamics

    !> Internal pure derivative for 2-atom Rydberg Hamiltonian
    pure subroutine rydberg_deriv(psi, omega_rad, delta_rad, v_rr_rad, dpsi)
        complex(dp), intent(in)  :: psi(4)
        real(dp),    intent(in)  :: omega_rad, delta_rad, v_rr_rad
        complex(dp), intent(out) :: dpsi(4)

        complex(dp) :: h_psi(4)
        complex(dp), parameter :: I_UNIT = (0.0_dp, 1.0_dp)

        ! H |gg> = 0.5 * Omega (|gr> + |rg>)
        h_psi(1) = 0.5_dp * omega_rad * (psi(2) + psi(3))

        ! H |gr> = 0.5 * Omega (|gg> + |rr>) - Delta |gr>
        h_psi(2) = 0.5_dp * omega_rad * (psi(1) + psi(4)) - delta_rad * psi(2)

        ! H |rg> = 0.5 * Omega (|gg> + |rr>) - Delta |rg>
        h_psi(3) = 0.5_dp * omega_rad * (psi(1) + psi(4)) - delta_rad * psi(3)

        ! H |rr> = 0.5 * Omega (|gr> + |rg>) + (V_rr - 2*Delta) |rr>
        h_psi(4) = 0.5_dp * omega_rad * (psi(2) + psi(3)) + (v_rr_rad - 2.0_dp * delta_rad) * psi(4)

        ! d psi / dt = -i H psi
        dpsi = - I_UNIT * h_psi
    end subroutine rydberg_deriv

    !> Calculate Z2 staggered order parameter: O_Z2 = (1/N) * sum_{i=1}^N (-1)^i * n_i
    pure function calc_z2_order_parameter(n_atoms, occ_arr) result(o_z2)
        integer,  intent(in) :: n_atoms
        real(dp), intent(in) :: occ_arr(n_atoms)
        real(dp) :: o_z2

        integer :: i
        real(dp) :: s

        s = 0.0_dp
        do i = 1, n_atoms
            if (mod(i, 2) == 1) then
                s = s - occ_arr(i)
            else
                s = s + occ_arr(i)
            end if
        end do
        o_z2 = s / real(n_atoms, dp)
    end function calc_z2_order_parameter

    !> Compute many-body scar dynamics and Z2 oscillation on a blockaded chain
    !> Uses effective constrained PXP model dynamics:
    !> O_Z2(t) exhibits coherent revivals at T_scar ~ 2*pi / (1.33 * Omega).
    pure subroutine calc_rydberg_scar_dynamics(cfg, atom, t_max_us, n_steps, &
                                              t_arr, z2_order_arr)
        type(rydberg_array_config_t), intent(in) :: cfg
        type(rydberg_atom_t),         intent(in) :: atom
        real(dp), intent(in)  :: t_max_us
        integer,  intent(in)  :: n_steps
        real(dp), intent(out) :: t_arr(n_steps)
        real(dp), intent(out) :: z2_order_arr(n_steps)

        real(dp) :: dt, t_curr, omega_eff, decay_rate
        real(dp) :: r_blockade, spacing
        integer :: it

        if (n_steps <= 1) return
        dt = t_max_us / real(n_steps - 1, dp)

        r_blockade = calc_rydberg_blockade_radius(atom, cfg%rabi_mhz)
        spacing = cfg%spacing_um

        ! Effective PXP scar oscillation frequency omega_scar ~ 1.33 * Omega
        omega_eff = TWOPI * 1.33_dp * cfg%rabi_mhz

        ! Damping rate due to finite lifetime and high-order blockade leakage
        decay_rate = 1.0_dp / atom%lifetime_us + 0.15_dp * max(0.0_dp, (spacing / r_blockade)**6)

        t_curr = 0.0_dp
        do it = 1, n_steps
            t_arr(it) = t_curr
            ! Z2 order starts at initial Neel state (|r g r g ...>, order = 0.5)
            ! and oscillates between positive and negative staggered order
            z2_order_arr(it) = 0.5_dp * cos(omega_eff * t_curr) * exp(-decay_rate * t_curr)
            t_curr = t_curr + dt
        end do
    end subroutine calc_rydberg_scar_dynamics

end module mod_rydberg_blockade
