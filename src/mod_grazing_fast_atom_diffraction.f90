! ==============================================================================
! GeneralModule: mod_grazing_fast_atom_diffraction.f90
!
! Grazing Incidence Fast Atom Diffraction (GIFAD) and
! Surface Rainbow Scattering
!
! Theoretical Foundations:
!   1. Fast atom diffraction (keV beams, grazing angles theta_in < 1-2 degrees).
!      Decoupling of fast longitudinal channeling motion along crystallographic axes:
!      Transverse quantum energy E_perp = E_tot * sin^2(theta_in) ~ meV to eV.
!      Transverse de Broglie wavelength lambda_perp = h / sqrt(2*M*E_perp) ~ 0.5-2 Angstrom.
!   2. Transverse 1D diffraction: discrete Bragg peaks at sin(theta_m) = m * lambda_perp / a_x.
!   3. Classical surface rainbow scattering and quantum Airy supernumerary peaks:
!      Rainbow angle theta_R = 2*pi * zeta / a_x.
!   4. Sub-picometer topological surface corrugation reconstruction.
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_grazing_fast_atom_diffraction
    use mod_constants, only: dp, PI, TWOPI, H_PLANCK, AMU2AU, ANG2AU, AU2EV, EV2AU
    implicit none
    private

    public :: gifad_experiment_t
    public :: gifad_spectrum_t
    public :: init_gifad_experiment
    public :: calc_gifad_transverse_kinematics
    public :: calc_gifad_rainbow_angle
    public :: calc_gifad_diffraction_spectrum
    public :: calc_surface_corrugation_from_rainbow

    !> GIFAD grazing beam and channel configuration
    type :: gifad_experiment_t
        character(len=16) :: projectile_name  !< Projectile species (e.g. He, Ne, H)
        real(dp) :: mass_amu                  !< Projectile mass [amu]
        real(dp) :: e_beam_kev                !< Total incident beam energy [keV]
        real(dp) :: theta_in_deg              !< Grazing incidence angle theta_in [deg]
        real(dp) :: ax_channel_ang            !< Transverse channel spacing a_x [Angstrom]
        real(dp) :: corrugation_ang           !< Corrugation amplitude zeta [Angstrom]
        real(dp) :: e_perp_ev                 !< Transverse kinetic energy E_perp [eV]
        real(dp) :: lambda_perp_ang           !< Transverse de Broglie wavelength [Angstrom]
    end type gifad_experiment_t

    !> Angular diffraction spectrum and rainbow parameters
    type :: gifad_spectrum_t
        integer  :: max_order                 !< Maximum diffraction order considered
        integer  :: n_open_orders             !< Number of allowed diffraction peaks
        real(dp) :: theta_rainbow_deg         !< Classical rainbow deflection angle [deg]
        real(dp) :: angles_deg(61)            !< Outgoing deflection angles theta_m [deg]
        real(dp) :: intensities(61)           !< Normalized intensities I(m)
        integer  :: orders(61)                !< Diffraction orders m
    end type gifad_spectrum_t

contains

    !> Initialize GIFAD experimental parameters
    pure subroutine init_gifad_experiment(exp_cfg, projectile, mass_amu, e_kev, &
                                         theta_deg, ax_ang, corrugation_ang)
        type(gifad_experiment_t), intent(out) :: exp_cfg
        character(len=*),         intent(in)  :: projectile
        real(dp),                 intent(in)  :: mass_amu
        real(dp),                 intent(in)  :: e_kev
        real(dp),                 intent(in)  :: theta_deg
        real(dp),                 intent(in)  :: ax_ang
        real(dp),                 intent(in)  :: corrugation_ang

        real(dp) :: e_tot_ev, th_rad, mass_kg

        exp_cfg%projectile_name = trim(adjustl(projectile))
        exp_cfg%mass_amu = max(1.0_dp, mass_amu)
        exp_cfg%e_beam_kev = max(0.01_dp, e_kev)
        exp_cfg%theta_in_deg = max(0.01_dp, theta_deg)
        exp_cfg%ax_channel_ang = max(0.5_dp, ax_ang)
        exp_cfg%corrugation_ang = max(0.001_dp, corrugation_ang)

        ! Total energy in eV: E_tot = E_kev * 1000 eV
        e_tot_ev = exp_cfg%e_beam_kev * 1000.0_dp
        th_rad = exp_cfg%theta_in_deg * (PI / 180.0_dp)

        ! Transverse energy: E_perp = E_tot * sin^2(theta_in)
        exp_cfg%e_perp_ev = e_tot_ev * (sin(th_rad)**2)

        ! Transverse de Broglie wavelength: lambda_perp = h / sqrt(2 * M * E_perp)
        mass_kg = exp_cfg%mass_amu * 1.66053906660e-27_dp
        exp_cfg%lambda_perp_ang = (H_PLANCK / sqrt(2.0_dp * mass_kg * &
            max(1.0e-25_dp, exp_cfg%e_perp_ev * 1.602176634e-19_dp))) * 1.0e10_dp
    end subroutine init_gifad_experiment

    !> Calculate transverse kinematics (E_perp and lambda_perp)
    pure subroutine calc_gifad_transverse_kinematics(exp_cfg, e_perp_ev, lambda_perp_ang)
        type(gifad_experiment_t), intent(in)  :: exp_cfg
        real(dp),                 intent(out) :: e_perp_ev
        real(dp),                 intent(out) :: lambda_perp_ang

        e_perp_ev = exp_cfg%e_perp_ev
        lambda_perp_ang = exp_cfg%lambda_perp_ang
    end subroutine calc_gifad_transverse_kinematics

    !> Compute classical surface rainbow deflection angle theta_R in degrees
    !> For sinusoidal channel corrugation: theta_R = 2*pi * (zeta / a_x) * (180 / pi)
    pure function calc_gifad_rainbow_angle(exp_cfg) result(theta_rainbow_deg)
        type(gifad_experiment_t), intent(in) :: exp_cfg
        real(dp) :: theta_rainbow_deg

        real(dp) :: slope_max

        ! Maximum corrugation slope: d xi / dx = (2*pi*zeta / ax)
        slope_max = TWOPI * (exp_cfg%corrugation_ang / exp_cfg%ax_channel_ang)
        theta_rainbow_deg = atan(slope_max) * (180.0_dp / PI)
    end function calc_gifad_rainbow_angle

    !> Compute 1D transverse diffraction spectrum with quantum rainbow envelope
    pure subroutine calc_gifad_diffraction_spectrum(exp_cfg, max_order, spec)
        type(gifad_experiment_t), intent(in)  :: exp_cfg
        integer,                  intent(in)  :: max_order
        type(gifad_spectrum_t),   intent(out) :: spec

        real(dp) :: th_r, sin_th_m, th_m_deg, arg_bessel
        real(dp) :: j_m, intensity, sum_int
        integer  :: m, idx, m_clamped

        m_clamped = min(30, max(1, max_order))
        spec%max_order = m_clamped
        th_r = calc_gifad_rainbow_angle(exp_cfg)
        spec%theta_rainbow_deg = th_r

        idx = 0
        sum_int = 0.0_dp

        ! Argument for Bessel envelope representing rainbow scattering:
        ! u = k_perp * zeta = (2*pi / lambda_perp) * zeta
        arg_bessel = TWOPI * (exp_cfg%corrugation_ang / exp_cfg%lambda_perp_ang)

        do m = -m_clamped, m_clamped
            idx = idx + 1
            spec%orders(idx) = m

            ! Transverse grating equation: sin(theta_m) = m * lambda_perp / a_x
            sin_th_m = real(m, dp) * exp_cfg%lambda_perp_ang / exp_cfg%ax_channel_ang

            if (abs(sin_th_m) <= 1.0_dp) then
                th_m_deg = asin(sin_th_m) * (180.0_dp / PI)
                spec%angles_deg(idx) = th_m_deg

                ! Hard-wall Eikonal approximation for 1D corrugation: I(m) ~ J_m(k_perp * zeta)^2
                j_m = bessel_jn_gifad(m, arg_bessel)
                intensity = j_m**2

                spec%intensities(idx) = intensity
                sum_int = sum_int + intensity
            else
                spec%angles_deg(idx) = 90.0_dp * sign(1.0_dp, sin_th_m)
                spec%intensities(idx) = 0.0_dp
            end if
        end do
        spec%n_open_orders = idx

        ! Normalize total diffraction intensity to 1.0
        if (sum_int > 1.0e-14_dp) then
            do m = 1, idx
                spec%intensities(m) = spec%intensities(m) / sum_int
            end do
        end if
    end subroutine calc_gifad_diffraction_spectrum

    !> Simple robust pure Bessel function of integer order J_m(x)
    pure function bessel_jn_gifad(m_order, x_val) result(res)
        integer,  intent(in) :: m_order
        real(dp), intent(in) :: x_val
        real(dp) :: res

        integer :: m_abs, k
        real(dp) :: term, s, x_half

        m_abs = abs(m_order)
        x_half = 0.5_dp * x_val

        term = 1.0_dp
        do k = 1, m_abs
            term = term * x_half / real(k, dp)
        end do

        s = term
        do k = 1, 24
            term = - term * (x_half**2) / (real(k, dp) * real(m_abs + k, dp))
            s = s + term
            if (abs(term) < 1.0e-16_dp * abs(s)) exit
        end do

        if (m_order < 0 .and. mod(m_abs, 2) == 1) then
            res = - s
        else
            res = s
        end if
    end function bessel_jn_gifad

    !> Inverse reconstruction of surface corrugation amplitude zeta from rainbow angle
    pure function calc_surface_corrugation_from_rainbow(ax_channel_ang, theta_rainbow_deg) &
        result(corrugation_ang)
        real(dp), intent(in) :: ax_channel_ang
        real(dp), intent(in) :: theta_rainbow_deg
        real(dp) :: corrugation_ang

        real(dp) :: th_r_rad

        th_r_rad = theta_rainbow_deg * (PI / 180.0_dp)
        ! zeta = a_x * tan(theta_R) / (2*pi)
        corrugation_ang = (ax_channel_ang * tan(th_r_rad)) / TWOPI
    end function calc_surface_corrugation_from_rainbow

end module mod_grazing_fast_atom_diffraction
