! ==============================================================================
! GeneralModule: mod_spinor_bec.f90
!
! Ultracold Spinor Bose-Einstein Condensates (Spin-1 / F=1) & Many-Body Spin Mixing
!
! Theoretical Foundations:
!   - T.-L. Ho, Phys. Rev. Lett. 81, 742 (1998) [Spinor Bose condensates]
!   - T. Ohmi & K. Machida, J. Phys. Soc. Jpn. 67, 1822 (1998) [Bose-Einstein condensation]
!   - M.-S. Chang et al., Phys. Rev. Lett. 92, 140403 (2004) [Spinor dynamics in BEC]
!   - H. Pu et al., Phys. Rev. A 60, 1463 (1999) [Spin mixing in optical traps]
!
! Key Capabilities:
!   1. Spin-dependent interaction constants c_0 and c_2 from s-wave scattering lengths a_0, a_2
!   2. Ferromagnetic (87Rb, c_2 < 0) vs Antiferromagnetic/Polar (23Na, c_2 > 0) classification
!   3. Single-Mode Approximation (SMA) non-linear spin-mixing coupled equations
!   4. Quadratic Zeeman effect (QZE) q_Z tuning and quantum phase transitions
!   5. RK4 coherent spin population oscillation propagator conserving m_z and norm
!
! Standard: Fortran 2008
! ==============================================================================

module mod_spinor_bec
    use mod_constants, only: dp, PI, TWOPI, HBAR, AMU2AU, GAUSS2AU
    implicit none
    private

    ! Public derived types
    public :: spinor_param_t, spinor_state_t

    ! Public procedures
    public :: init_spinor_preset
    public :: calc_spinor_interaction_couplings
    public :: calc_quadratic_zeeman_shift
    public :: propagate_spinor_sma_rk4
    public :: calc_spinor_energy
    public :: simulate_spin_mixing_dynamics

    ! Spinor condensate parameters
    type :: spinor_param_t
        character(len=8) :: name_atom      ! "87Rb" or "23Na"
        real(dp)         :: mass_au        ! Atomic mass in a.u.
        real(dp)         :: a0_au          ! Scattering length a_0 in a.u. (total spin F=0)
        real(dp)         :: a2_au          ! Scattering length a_2 in a.u. (total spin F=2)
        real(dp)         :: c0_au          ! Density-density interaction c_0 = 4*pi*(a_0 + 2*a_2) / (3*m)
        real(dp)         :: c2_au          ! Spin-exchange interaction c_2 = 4*pi*(a_2 - a_0) / (3*m)
        real(dp)         :: q_zeeman_au    ! Quadratic Zeeman shift q_Z in a.u.
        real(dp)         :: density_au     ! Average spatial mode overlap density n in a.u.
        logical          :: is_ferromagnet ! True if c_2 < 0 (ferromagnetic, like 87Rb)
    end type spinor_param_t

    ! Spin-1 macroscopic wave function components (zeta_+1, zeta_0, zeta_-1)
    type :: spinor_state_t
        complex(dp) :: zeta_p1             ! m = +1 component
        complex(dp) :: zeta_0              ! m = 0 component
        complex(dp) :: zeta_m1             ! m = -1 component
        real(dp)    :: pop_p1              ! Population |zeta_+1|^2
        real(dp)    :: pop_0               ! Population |zeta_0|^2
        real(dp)    :: pop_m1              ! Population |zeta_-1|^2
        real(dp)    :: total_norm          ! Total norm |z_+1|^2 + |z_0|^2 + |z_-1|^2 = 1.0
        real(dp)    :: magnetization_mz    ! Magnetization m_z = |z_+1|^2 - |z_-1|^2 (conserved)
    end type spinor_state_t

    ! In atomic units, hbar = 1.0
    real(dp), parameter :: HBAR_AU = 1.0_dp

contains

    ! ==========================================================================
    ! init_spinor_preset:
    ! Pre-loads well-established experimental parameters for 87Rb and 23Na
    ! 87Rb: a0 = 101.8 a0, a2 = 100.4 a0 => a2 - a0 < 0 (Ferromagnetic)
    ! 23Na: a0 = 50.0 a0, a2 = 54.5 a0  => a2 - a0 > 0 (Polar / Antiferromagnetic)
    ! ==========================================================================
    subroutine init_spinor_preset(atom_name, b_field_gauss, density_cm3, param, stat)
        character(len=*), intent(in)     :: atom_name
        real(dp), intent(in)             :: b_field_gauss
        real(dp), intent(in)             :: density_cm3
        type(spinor_param_t), intent(out):: param
        integer, intent(out)             :: stat

        real(dp) :: density_au, b_field_au, delta_ehfs_ghz

        stat = 0
        density_au = density_cm3 * ((5.29177210903e-9_dp)**3)
        b_field_au = b_field_gauss * GAUSS2AU

        if (trim(atom_name) == "87Rb" .or. trim(atom_name) == "Rb87") then
            param%name_atom = "87Rb"
            param%mass_au = 87.0_dp * AMU2AU
            param%a0_au = 101.8_dp
            param%a2_au = 100.4_dp
            delta_ehfs_ghz = 6.83468_dp
        else if (trim(atom_name) == "23Na" .or. trim(atom_name) == "Na23") then
            param%name_atom = "23Na"
            param%mass_au = 23.0_dp * AMU2AU
            param%a0_au = 50.0_dp
            param%a2_au = 54.5_dp
            delta_ehfs_ghz = 1.7716_dp
        else
            stat = 1
            return
        end if

        param%density_au = density_au
        call calc_spinor_interaction_couplings(param)
        param%q_zeeman_au = calc_quadratic_zeeman_shift(b_field_gauss, delta_ehfs_ghz)
    end subroutine init_spinor_preset

    ! ==========================================================================
    ! calc_spinor_interaction_couplings:
    ! Computes c_0 and c_2 in a.u.:
    !   c_0 = 4*pi*hbar^2 * (a_0 + 2*a_2) / (3*m)
    !   c_2 = 4*pi*hbar^2 * (a_2 - a_0) / (3*m)
    ! ==========================================================================
    pure subroutine calc_spinor_interaction_couplings(param)
        type(spinor_param_t), intent(inout) :: param

        real(dp) :: prefactor

        prefactor = (4.0_dp * PI * (HBAR_AU**2)) / (3.0_dp * param%mass_au)
        param%c0_au = prefactor * (param%a0_au + 2.0_dp * param%a2_au)
        param%c2_au = prefactor * (param%a2_au - param%a0_au)
        param%is_ferromagnet = (param%c2_au < 0.0_dp)
    end subroutine calc_spinor_interaction_couplings

    ! ==========================================================================
    ! calc_quadratic_zeeman_shift:
    ! Quadratic Zeeman energy q_Z ~ (g_F * mu_B * B)^2 / (4 * Delta_E_hfs)
    ! ==========================================================================
    pure function calc_quadratic_zeeman_shift(b_field_gauss, delta_ehfs_ghz) result(qz_au)
        real(dp), intent(in) :: b_field_gauss
        real(dp), intent(in) :: delta_ehfs_ghz
        real(dp) :: qz_au

        real(dp) :: b_au, ehfs_au, mu_term

        b_au = b_field_gauss * GAUSS2AU
        ehfs_au = delta_ehfs_ghz * 1.0e9_dp * (2.4188843265857e-17_dp) * TWOPI
        ! g_F = -1/2 for F=1 Rb87, mu_B in a.u. is 0.5 => g_F * mu_B = 0.25
        mu_term = 0.25_dp * b_au

        if (ehfs_au > 1.0e-15_dp) then
            qz_au = (mu_term**2) / (4.0_dp * ehfs_au)
        else
            qz_au = 0.0_dp
        end if
    end function calc_quadratic_zeeman_shift

    ! ==========================================================================
    ! calc_spinor_derivatives:
    ! Computes time derivatives d(zeta)/dt for the SMA coupled equations:
    !   i * d(zeta_p1)/dt = q_Z*zeta_p1 + c2' * [ (|z_p1|^2 + |z_0|^2 - |z_m1|^2)*z_p1 + z_0^2 * conjg(z_m1) ]
    !   i * d(zeta_0)/dt  = c2' * [ (|z_p1|^2 + |z_m1|^2)*z_0 + 2*z_p1*z_m1*conjg(z_0) ]
    !   i * d(zeta_m1)/dt = q_Z*zeta_m1 + c2' * [ (|z_m1|^2 + |z_0|^2 - |z_p1|^2)*z_m1 + z_0^2 * conjg(z_p1) ]
    ! ==========================================================================
    pure subroutine calc_spinor_derivatives(zp1, z0, zm1, c2_eff, qz, dzp1, dz0, dzm1)
        complex(dp), intent(in)  :: zp1, z0, zm1
        real(dp), intent(in)     :: c2_eff, qz
        complex(dp), intent(out) :: dzp1, dz0, dzm1

        real(dp) :: p1, p0, pm1
        complex(dp) :: eye_inv

        ! 1 / i = -i
        eye_inv = (0.0_dp, -1.0_dp)

        p1  = abs(zp1)**2
        p0  = abs(z0)**2
        pm1 = abs(zm1)**2

        ! d(zeta_+1)/dt
        dzp1 = eye_inv * ( qz * zp1 + c2_eff * ( (p1 + p0 - pm1) * zp1 + (z0**2) * conjg(zm1) ) )

        ! d(zeta_0)/dt
        dz0  = eye_inv * ( c2_eff * ( (p1 + pm1) * z0 + 2.0_dp * zp1 * zm1 * conjg(z0) ) )

        ! d(zeta_-1)/dt
        dzm1 = eye_inv * ( qz * zm1 + c2_eff * ( (pm1 + p0 - p1) * zm1 + (z0**2) * conjg(zp1) ) )
    end subroutine calc_spinor_derivatives

    ! ==========================================================================
    ! propagate_spinor_sma_rk4:
    ! Fourth-order Runge-Kutta step for the spinor SMA wave function
    ! ==========================================================================
    pure subroutine propagate_spinor_sma_rk4(param, dt_au, state)
        type(spinor_param_t), intent(in) :: param
        real(dp), intent(in)             :: dt_au
        type(spinor_state_t), intent(inout) :: state

        complex(dp) :: k1_p1, k1_0, k1_m1
        complex(dp) :: k2_p1, k2_0, k2_m1
        complex(dp) :: k3_p1, k3_0, k3_m1
        complex(dp) :: k4_p1, k4_0, k4_m1
        complex(dp) :: z1, z0, zm
        real(dp)    :: c2_eff, qz

        c2_eff = param%c2_au * param%density_au
        qz     = param%q_zeeman_au

        z1 = state%zeta_p1
        z0 = state%zeta_0
        zm = state%zeta_m1

        ! k1
        call calc_spinor_derivatives(z1, z0, zm, c2_eff, qz, k1_p1, k1_0, k1_m1)

        ! k2
        call calc_spinor_derivatives(z1 + 0.5_dp*dt_au*k1_p1, z0 + 0.5_dp*dt_au*k1_0, &
                                    zm + 0.5_dp*dt_au*k1_m1, c2_eff, qz, k2_p1, k2_0, k2_m1)

        ! k3
        call calc_spinor_derivatives(z1 + 0.5_dp*dt_au*k2_p1, z0 + 0.5_dp*dt_au*k2_0, &
                                    zm + 0.5_dp*dt_au*k2_m1, c2_eff, qz, k3_p1, k3_0, k3_m1)

        ! k4
        call calc_spinor_derivatives(z1 + dt_au*k3_p1, z0 + dt_au*k3_0, &
                                    zm + dt_au*k3_m1, c2_eff, qz, k4_p1, k4_0, k4_m1)

        ! Update
        state%zeta_p1 = z1 + (dt_au / 6.0_dp) * (k1_p1 + 2.0_dp*k2_p1 + 2.0_dp*k3_p1 + k4_p1)
        state%zeta_0  = z0 + (dt_au / 6.0_dp) * (k1_0  + 2.0_dp*k2_0  + 2.0_dp*k3_0  + k4_0)
        state%zeta_m1 = zm + (dt_au / 6.0_dp) * (k1_m1 + 2.0_dp*k2_m1 + 2.0_dp*k3_m1 + k4_m1)

        ! Update populations and observables
        state%pop_p1 = abs(state%zeta_p1)**2
        state%pop_0  = abs(state%zeta_0)**2
        state%pop_m1 = abs(state%zeta_m1)**2
        state%total_norm = state%pop_p1 + state%pop_0 + state%pop_m1
        state%magnetization_mz = state%pop_p1 - state%pop_m1
    end subroutine propagate_spinor_sma_rk4

    ! ==========================================================================
    ! calc_spinor_energy:
    ! Evaluates the SMA mean-field energy per particle in a.u.
    ! ==========================================================================
    pure function calc_spinor_energy(param, state) result(energy_au)
        type(spinor_param_t), intent(in) :: param
        type(spinor_state_t), intent(in) :: state
        real(dp) :: energy_au

        real(dp) :: c2_eff, qz, s_perp_term

        c2_eff = param%c2_au * param%density_au
        qz     = param%q_zeeman_au

        ! Spin exchange term: 2 * Re[ (z_0)^2 * conjg(z_p1) * conjg(z_m1) ]
        s_perp_term = 2.0_dp * real((state%zeta_0**2) * conjg(state%zeta_p1) * conjg(state%zeta_m1), dp)

        energy_au = qz * (state%pop_p1 + state%pop_m1) + &
                    0.5_dp * c2_eff * ((state%pop_p1 - state%pop_m1)**2 + &
                    2.0_dp * (state%pop_p1 + state%pop_m1) * state%pop_0 + s_perp_term)
    end function calc_spinor_energy

    ! ==========================================================================
    ! simulate_spin_mixing_dynamics:
    ! Simulates time series of spinor populations starting from initial state
    ! ==========================================================================
    subroutine simulate_spin_mixing_dynamics(param, init_state, t_total_au, n_steps, &
                                            t_arr, pop_0_arr, pop_side_arr, stat)
        type(spinor_param_t), intent(in) :: param
        type(spinor_state_t), intent(in) :: init_state
        real(dp), intent(in)             :: t_total_au
        integer, intent(in)              :: n_steps
        real(dp), intent(out)            :: t_arr(n_steps)
        real(dp), intent(out)            :: pop_0_arr(n_steps)
        real(dp), intent(out)            :: pop_side_arr(n_steps)
        integer, intent(out)             :: stat

        type(spinor_state_t) :: state
        real(dp) :: dt
        integer :: i

        stat = 0
        if (n_steps <= 1 .or. t_total_au <= 0.0_dp) then
            stat = 1
            return
        end if

        dt = t_total_au / real(n_steps - 1, dp)
        state = init_state

        do i = 1, n_steps
            t_arr(i) = real(i - 1, dp) * dt
            pop_0_arr(i) = state%pop_0
            pop_side_arr(i) = 0.5_dp * (state%pop_p1 + state%pop_m1)

            if (i < n_steps) then
                call propagate_spinor_sma_rk4(param, dt, state)
            end if
        end do
    end subroutine simulate_spin_mixing_dynamics

end module mod_spinor_bec
