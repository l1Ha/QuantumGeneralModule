! ==============================================================================
! GeneralModule: mod_dipolar_droplets_lhy.f90
!
! Ultracold Dipolar Quantum Droplets & Lee-Huang-Yang (LHY) Quantum Fluctuations
!
! Theoretical Foundations:
!   - T. D. Lee, K. Huang, & C. N. Yang, Phys. Rev. 106, 1135 (1957) [LHY quantum energy]
!   - D. S. Petrov, Phys. Rev. Lett. 115, 155301 (2015) [Droplet stabilization by LHY]
!   - H. Kadau et al., Nature 530, 194 (2016) [Observation of quantum ferrofluid droplets]
!   - I. Ferrier-Barbut et al., Phys. Rev. Lett. 116, 215301 (2016) [Self-bound dipolar droplets]
!   - A. R. P. Lima & A. Pelster, Phys. Rev. A 84, 041604(R) (2011) [Dipolar LHY corrections]
!   - F. Chomaz et al., Rep. Prog. Phys. 86, 026401 (2023) [Dipolar physics review]
!
! Key Capabilities:
!   1. Magnetic dipole length a_dd and relative dipolar strength epsilon_dd = a_dd / a_s
!   2. Pelster-Lima dipolar LHY auxiliary function Q_5(epsilon_dd)
!   3. Self-bound equilibrium droplet density n_0 and negative chemical potential mu_0 < 0
!   4. Critical atom number N_crit for droplet self-binding vs evaporation
!   5. Extended Gross-Pitaevskii (eGPE) local energy density and chemical potential
!
! Standard: Fortran 2008
! ==============================================================================

module mod_dipolar_droplets_lhy
    use mod_constants, only: dp, PI, TWOPI, HBAR, AMU2AU
    implicit none
    private

    ! Public derived types
    public :: dipolar_droplet_param_t, droplet_state_t

    ! Public procedures
    public :: init_dipolar_droplet_param
    public :: calc_pelster_lima_q5
    public :: calc_equilibrium_droplet_density
    public :: calc_droplet_chemical_potential
    public :: calc_critical_atom_number
    public :: calc_egpe_energy_density

    ! Dipolar droplet physical parameters
    type :: dipolar_droplet_param_t
        character(len=8) :: atom_name        ! e.g., "162Dy" or "166Er"
        real(dp)         :: mass_au          ! Atomic mass in a.u.
        real(dp)         :: dipole_mu_b      ! Magnetic dipole moment in Bohr magnetons (mu_B)
        real(dp)         :: a_dd_au          ! Dipole length a_dd in a.u.
        real(dp)         :: a_scat_au        ! s-wave contact scattering length a_s in a.u.
        real(dp)         :: epsilon_dd       ! Relative dipolar strength epsilon_dd = a_dd / a_s
        real(dp)         :: g_contact_au     ! Contact interaction g = 4*pi*hbar^2*a_s / m
        real(dp)         :: q5_factor        ! Pelster-Lima LHY correction factor Q_5(epsilon_dd)
        real(dp)         :: gamma_lhy_au     ! LHY coefficient gamma_LHY = (32 / 3*sqrt(pi)) * g * a_s^(3/2) * Q_5
        logical          :: is_droplet_regime! True if epsilon_dd > 1.0 (mean-field collapse prevented by LHY)
    end type dipolar_droplet_param_t

    ! Droplet equilibrium state properties
    type :: droplet_state_t
        real(dp) :: n0_density_au     ! Peak / saturation equilibrium density in a.u.
        real(dp) :: n0_density_cm3    ! Density in cm^-3
        real(dp) :: mu_chem_pot_au    ! Chemical potential mu in a.u. (negative for self-bound state)
        real(dp) :: energy_per_particle_au ! Bound energy E/N in a.u.
        real(dp) :: n_critical        ! Critical particle number N_crit
        real(dp) :: healing_length_au ! Quantum droplet healing length xi = hbar / sqrt(m * |mu|)
    end type droplet_state_t

    ! Fundamental constants in atomic units:
    ! hbar = 1.0, m_e = 1.0, mu_B = 0.5 a.u.
    real(dp), parameter :: HBAR_AU = 1.0_dp
    real(dp), parameter :: MU_B_AU = 0.5_dp
    ! Magnetic permeability in a.u.: mu_0 = 4*pi * alpha^2 ~ 4*pi * (1/137.035999)^2
    ! In atomic units, fine-structure constant alpha = 1 / 137.035999084
    real(dp), parameter :: ALPHA_FS = 1.0_dp / 137.035999084_dp
    real(dp), parameter :: MU0_AU   = 4.0_dp * PI * (ALPHA_FS**2)

    ! Length conversion: 1 a.u. = 5.29177210903e-9 cm
    real(dp), parameter :: AU_TO_CM = 5.29177210903e-9_dp

contains

    ! ==========================================================================
    ! init_dipolar_droplet_param:
    ! Initializes parameters for magnetic dipolar atoms (162Dy: 10 mu_B, 166Er: 7 mu_B)
    ! ==========================================================================
    subroutine init_dipolar_droplet_param(atom_name, a_scat_bohr, param, stat)
        character(len=*), intent(in)          :: atom_name
        real(dp), intent(in)                  :: a_scat_bohr
        type(dipolar_droplet_param_t), intent(out) :: param
        integer, intent(out)                  :: stat

        real(dp) :: mu_mag_au

        stat = 0
        if (a_scat_bohr <= 0.0_dp) then
            stat = 1
            return
        end if

        if (trim(atom_name) == "162Dy" .or. trim(atom_name) == "Dy162") then
            param%atom_name   = "162Dy"
            param%mass_au     = 161.926798_dp * AMU2AU
            param%dipole_mu_b = 9.93_dp  ! Dysprosium has the largest magnetic moment in periodic table (~ 10 mu_B)
        else if (trim(atom_name) == "166Er" .or. trim(atom_name) == "Er166") then
            param%atom_name   = "166Er"
            param%mass_au     = 165.930293_dp * AMU2AU
            param%dipole_mu_b = 7.0_dp   ! Erbium ~ 7 mu_B
        else
            ! Generic default
            param%atom_name   = trim(atom_name)
            param%mass_au     = 162.0_dp * AMU2AU
            param%dipole_mu_b = 10.0_dp
        end if

        param%a_scat_au = a_scat_bohr

        ! Magnetic dipole moment in a.u. = dipole_mu_b * MU_B_AU
        mu_mag_au = param%dipole_mu_b * MU_B_AU

        ! Dipolar length: a_dd = mu_0 * mu^2 * m / (12 * pi * hbar^2)
        param%a_dd_au = (MU0_AU * (mu_mag_au**2) * param%mass_au) / (12.0_dp * PI * (HBAR_AU**2))

        ! Relative dipolar strength epsilon_dd = a_dd / a_s
        param%epsilon_dd = param%a_dd_au / param%a_scat_au
        param%is_droplet_regime = (param%epsilon_dd > 1.0_dp)

        ! Contact interaction parameter g = 4 * pi * hbar^2 * a_s / m
        param%g_contact_au = (4.0_dp * PI * (HBAR_AU**2) * param%a_scat_au) / param%mass_au

        ! Pelster-Lima Q5 factor
        param%q5_factor = calc_pelster_lima_q5(param%epsilon_dd)

        ! LHY coefficient gamma_LHY:
        !   Delta E_LHY / V = (64 / (15*sqrt(pi))) * (hbar^2 / m) * a_s^(5/2) * Q_5 * n^(5/2)
        !   Chemical potential contribution: Delta mu_LHY = (32 / (3*sqrt(pi))) * g * a_s^(3/2) * Q_5 * n^(3/2)
        param%gamma_lhy_au = (32.0_dp / (3.0_dp * sqrt(PI))) * param%g_contact_au * &
                             (param%a_scat_au**1.5_dp) * param%q5_factor
    end subroutine init_dipolar_droplet_param

    ! ==========================================================================
    ! calc_pelster_lima_q5:
    ! Evaluates the Pelster-Lima beyond-mean-field dipolar integral:
    !   Q_5(eps) = (3/2) \int_0^1 du (1 - u^2) [ 1 + eps*(3*u^2 - 1) ]^(5/2)
    ! Accurate analytical expansion: Q_5(eps) ~ 1 + (3/2)*eps^2 + (3/8)*eps^3 ...
    ! ==========================================================================
    pure function calc_pelster_lima_q5(eps) result(q5)
        real(dp), intent(in) :: eps
        real(dp) :: q5

        integer, parameter :: N_QUAD = 200
        real(dp) :: du, u, f_u, sum_int, bracket
        integer :: i

        if (eps <= 0.0_dp) then
            q5 = 1.0_dp
            return
        end if

        du = 1.0_dp / real(N_QUAD, dp)
        sum_int = 0.0_dp

        do i = 1, N_QUAD
            u = (real(i, dp) - 0.5_dp) * du
            bracket = 1.0_dp + eps * (3.0_dp * u**2 - 1.0_dp)
            if (bracket > 0.0_dp) then
                f_u = (1.0_dp - u**2) * (bracket**2.5_dp)
                sum_int = sum_int + f_u * du
            end if
        end do

        q5 = 1.5_dp * sum_int
    end function calc_pelster_lima_q5

    ! ==========================================================================
    ! calc_equilibrium_droplet_density:
    ! Computes the self-bound flat-top equilibrium density n_0 of a dipolar droplet
    ! where mean-field attraction is balanced by repulsive LHY quantum fluctuations:
    !   n_0 = [ 15*sqrt(pi) / (128 * a_s^(5/2) * Q_5) * (m*g_eff / hbar^2) ]^2
    !       = (25 * pi / 1024) * [ (epsilon_dd - 1)^2 / (a_s^3 * Q_5^2) ]
    ! ==========================================================================
    pure function calc_equilibrium_droplet_density(param) result(n0_au)
        type(dipolar_droplet_param_t), intent(in) :: param
        real(dp) :: n0_au

        real(dp) :: g_eff

        if (param%epsilon_dd <= 1.0_dp .or. param%gamma_lhy_au <= 0.0_dp) then
            ! Below threshold, no self-bound droplet exists in free space
            n0_au = 0.0_dp
            return
        end if

        g_eff = param%g_contact_au * (1.0_dp - param%epsilon_dd)
        n0_au = (25.0_dp / 36.0_dp) * ((abs(g_eff) / param%gamma_lhy_au)**2)
    end function calc_equilibrium_droplet_density

    ! ==========================================================================
    ! calc_droplet_chemical_potential:
    ! Evaluates the chemical potential at density n:
    !   mu(n) = g * (1 - epsilon_dd) * n + gamma_LHY * n^(3/2)
    ! At equilibrium density n_0, mu(n_0) < 0 ensures self-binding.
    ! ==========================================================================
    pure function calc_droplet_chemical_potential(param, density_au) result(mu_au)
        type(dipolar_droplet_param_t), intent(in) :: param
        real(dp), intent(in)                      :: density_au
        real(dp)                                  :: mu_au

        real(dp) :: g_eff

        if (density_au <= 0.0_dp) then
            mu_au = 0.0_dp
            return
        end if

        ! Effective mean-field coupling (negative for epsilon_dd > 1 along dipole axis)
        g_eff = param%g_contact_au * (1.0_dp - param%epsilon_dd)

        mu_au = g_eff * density_au + param%gamma_lhy_au * (density_au**1.5_dp)
    end function calc_droplet_chemical_potential

    ! ==========================================================================
    ! calc_critical_atom_number:
    ! Estimates critical atom number N_crit for 3D self-bound droplet stability
    ! using Petrov-Chomaz scaling: N_crit ~ 18.6 * (a_s / (a_dd - a_s))^(5/2) / sqrt(Q_5)
    ! ==========================================================================
    pure function calc_critical_atom_number(param) result(n_crit)
        type(dipolar_droplet_param_t), intent(in) :: param
        real(dp) :: n_crit

        real(dp) :: delta_eps

        delta_eps = param%epsilon_dd - 1.0_dp
        if (delta_eps <= 0.0_dp) then
            n_crit = 1.0e30_dp  ! Cannot form stable droplet
            return
        end if

        n_crit = 18.6_dp * ((1.0_dp / delta_eps)**2.5_dp) / sqrt(param%q5_factor)
    end function calc_critical_atom_number

    ! ==========================================================================
    ! calc_egpe_energy_density:
    ! Evaluates the total local eGPE energy density:
    !   E/V = 0.5 * g*(1 - epsilon_dd)*n^2 + (2/5) * gamma_LHY * n^(5/2)
    ! ==========================================================================
    pure function calc_egpe_energy_density(param, density_au) result(energy_density_au)
        type(dipolar_droplet_param_t), intent(in) :: param
        real(dp), intent(in)                      :: density_au
        real(dp)                                  :: energy_density_au

        real(dp) :: g_eff

        if (density_au <= 0.0_dp) then
            energy_density_au = 0.0_dp
            return
        end if

        g_eff = param%g_contact_au * (1.0_dp - param%epsilon_dd)
        energy_density_au = 0.5_dp * g_eff * (density_au**2) + &
                            0.4_dp * param%gamma_lhy_au * (density_au**2.5_dp)
    end function calc_egpe_energy_density

end module mod_dipolar_droplets_lhy
