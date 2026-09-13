! ==============================================================================
! GeneralModule: mod_confined_scattering.f90
!
! Ultracold Quantum Scattering in Low-Dimensional Traps & Confinement-Induced
! Resonances (CIR)
!
! Theoretical Foundations:
!   - M. Olshanii, Phys. Rev. Lett. 81, 938 (1998) [Confinement-induced resonance]
!   - D. S. Petrov et al., Phys. Rev. Lett. 84, 2551 (2000) [Scattering in quasi-2D]
!   - T. Bergeman, M. G. Moore, & M. Olshanii, PRL 91, 163201 (2003) [Anisotropic CIR]
!   - E. Haller et al., Science 325, 1224 (2009) [Experimental observation of CIR]
!
! Key Capabilities:
!   1. Quasi-1D effective interaction g_1D, 1D scattering length a_1D, & CIR pole
!   2. Confinement-induced molecular dimer binding energy E_b in waveguides
!   3. Anisotropic transverse waveguide CIR resonance splitting (omega_x /= omega_y)
!   4. Quasi-2D effective scattering length a_2D and low-energy 2D scattering amplitude
!   5. Lieb-Liniger dimensionless interaction parameter gamma_LL & Tonks-Girardeau regime
!
! Standard: Fortran 2008
! ==============================================================================

module mod_confined_scattering
    use mod_constants, only: dp, PI, TWOPI, HBAR, EYE
    implicit none
    private

    ! Public derived types
    public :: waveguide_1d_t, planar_2d_t, cir_result_t

    ! Public procedures
    public :: init_waveguide_1d
    public :: init_planar_2d
    public :: calc_olshanii_cir_parameters
    public :: calc_confined_dimer_binding_energy
    public :: calc_anisotropic_cir_poles
    public :: calc_quasi_2d_scattering
    public :: calc_lieb_liniger_parameter

    ! 1D waveguide confinement configuration
    type :: waveguide_1d_t
        real(dp) :: omega_trans_au    ! Transverse harmonic frequency omega_perp in a.u.
        real(dp) :: a_perp_au         ! Transverse oscillator length a_perp = sqrt(hbar/(mu*w))
        real(dp) :: reduced_mass_au   ! Collision reduced mass mu in a.u.
        real(dp) :: a_cir_au          ! Resonance pole a_s,CIR = a_perp / C
        real(dp) :: anisotropy_ratio  ! eta = omega_x / omega_y (1.0 for isotropic)
    end type waveguide_1d_t

    ! 2D planar confinement configuration
    type :: planar_2d_t
        real(dp) :: omega_z_au        ! Axial harmonic frequency omega_z in a.u.
        real(dp) :: a_z_au            ! Axial oscillator length a_z = sqrt(hbar/(mu*w_z))
        real(dp) :: reduced_mass_au   ! Collision reduced mass mu in a.u.
    end type planar_2d_t

    ! CIR calculation result
    type :: cir_result_t
        real(dp) :: g_1d_au           ! 1D coupling constant g_1D (a.u.)
        real(dp) :: a_1d_au           ! 1D scattering length a_1D (a.u.)
        real(dp) :: binding_energy_au ! Induced dimer binding energy E_b (a.u.)
        real(dp) :: gamma_ll          ! Lieb-Liniger parameter for given 1D density
        logical  :: is_near_cir       ! True if |a_s - a_CIR| / a_CIR < 0.05
        logical  :: is_tonks_regime   ! True if gamma_ll > 10.0 (Tonks-Girardeau)
    end type cir_result_t

    ! Olshanii universal constant C = -zeta(1/2) / sqrt(2) ~ 1.0326012759
    ! where zeta(1/2) ~ -1.4603545088095868
    real(dp), parameter :: OLSHANII_C = 1.032601275928643_dp
    ! Petrov 2D prefactor: sqrt(pi / 0.905) ~ 1.86368
    real(dp), parameter :: PETROV_PREFACTOR_2D = 2.092_dp
    ! In atomic units, hbar = 1.0
    real(dp), parameter :: HBAR_AU = 1.0_dp

contains

    ! ==========================================================================
    ! init_waveguide_1d:
    ! Initializes a transverse harmonic waveguide with frequency omega_perp
    ! and collision reduced mass mu.
    ! ==========================================================================
    pure subroutine init_waveguide_1d(omega_perp_au, reduced_mass_au, wg, stat)
        real(dp), intent(in)            :: omega_perp_au
        real(dp), intent(in)            :: reduced_mass_au
        type(waveguide_1d_t), intent(out) :: wg
        integer, intent(out)            :: stat

        stat = 0
        if (omega_perp_au <= 0.0_dp .or. reduced_mass_au <= 0.0_dp) then
            stat = 1
            return
        end if

        wg%omega_trans_au = omega_perp_au
        wg%reduced_mass_au = reduced_mass_au
        wg%a_perp_au = sqrt(HBAR_AU / (reduced_mass_au * omega_perp_au))
        wg%a_cir_au = wg%a_perp_au / OLSHANII_C
        wg%anisotropy_ratio = 1.0_dp
    end subroutine init_waveguide_1d

    ! ==========================================================================
    ! init_planar_2d:
    ! Initializes an axial 1D harmonic trap for quasi-2D confinement
    ! ==========================================================================
    pure subroutine init_planar_2d(omega_z_au, reduced_mass_au, pl, stat)
        real(dp), intent(in)         :: omega_z_au
        real(dp), intent(in)         :: reduced_mass_au
        type(planar_2d_t), intent(out) :: pl
        integer, intent(out)         :: stat

        stat = 0
        if (omega_z_au <= 0.0_dp .or. reduced_mass_au <= 0.0_dp) then
            stat = 1
            return
        end if

        pl%omega_z_au = omega_z_au
        pl%reduced_mass_au = reduced_mass_au
        pl%a_z_au = sqrt(HBAR_AU / (reduced_mass_au * omega_z_au))
    end subroutine init_planar_2d

    ! ==========================================================================
    ! calc_olshanii_cir_parameters:
    ! Computes 1D effective interaction g_1D and 1D scattering length a_1D:
    !   g_1D = (2 * hbar^2 * a_s) / (mu * a_perp^2 * (1 - C * a_s / a_perp))
    !   a_1D = - (a_perp^2 / (2 * a_s)) * (1 - C * a_s / a_perp)
    ! ==========================================================================
    pure subroutine calc_olshanii_cir_parameters(wg, a_s_3d_au, n_1d_density_au, &
                                               res, stat)
        type(waveguide_1d_t), intent(in) :: wg
        real(dp), intent(in)             :: a_s_3d_au
        real(dp), intent(in)             :: n_1d_density_au
        type(cir_result_t), intent(out)  :: res
        integer, intent(out)             :: stat

        real(dp) :: denom, ratio_s

        stat = 0
        ratio_s = a_s_3d_au / wg%a_perp_au
        denom = 1.0_dp - OLSHANII_C * ratio_s

        res%is_near_cir = (abs(denom) < 0.05_dp)

        if (abs(denom) < 1.0e-12_dp) then
            ! Exactly on CIR resonance
            res%g_1d_au = 1.0e12_dp * sign(1.0_dp, a_s_3d_au)
            res%a_1d_au = 0.0_dp
            res%binding_energy_au = 0.0_dp
            res%gamma_ll = 1.0e12_dp
            res%is_tonks_regime = .true.
            return
        end if

        ! g_1D = 2 * hbar^2 * a_s / (mu * a_perp^2 * denom)
        res%g_1d_au = (2.0_dp * (HBAR_AU**2) * a_s_3d_au) / &
                      (wg%reduced_mass_au * (wg%a_perp_au**2) * denom)

        ! a_1D = -hbar^2 / (mu * g_1D) = - a_perp^2 / (2 * a_s) * denom
        if (abs(a_s_3d_au) > 1.0e-15_dp) then
            res%a_1d_au = - (wg%a_perp_au**2 / (2.0_dp * a_s_3d_au)) * denom
        else
            res%a_1d_au = -1.0e15_dp
        end if

        ! Induced dimer binding energy E_b = -hbar^2 / (2 * mu * a_1D^2)
        if (abs(res%a_1d_au) > 1.0e-14_dp) then
            res%binding_energy_au = - (HBAR_AU**2) / &
                                    (2.0_dp * wg%reduced_mass_au * (res%a_1d_au**2))
        else
            res%binding_energy_au = 0.0_dp
        end if

        ! Lieb-Liniger parameter gamma_LL
        if (n_1d_density_au > 0.0_dp .and. abs(res%a_1d_au) > 1.0e-15_dp) then
            res%gamma_ll = 2.0_dp / (n_1d_density_au * abs(res%a_1d_au))
        else
            res%gamma_ll = 0.0_dp
        end if

        res%is_tonks_regime = (res%gamma_ll > 10.0_dp)
    end subroutine calc_olshanii_cir_parameters

    ! ==========================================================================
    ! calc_confined_dimer_binding_energy:
    ! Returns the binding energy of the confinement-induced molecule E_b (a.u.)
    ! ==========================================================================
    pure function calc_confined_dimer_binding_energy(wg, a_s_3d_au) result(e_b)
        type(waveguide_1d_t), intent(in) :: wg
        real(dp), intent(in)             :: a_s_3d_au
        real(dp) :: e_b
        type(cir_result_t) :: res
        integer :: s

        call calc_olshanii_cir_parameters(wg, a_s_3d_au, 0.0_dp, res, s)
        e_b = res%binding_energy_au
    end function calc_confined_dimer_binding_energy

    ! ==========================================================================
    ! calc_anisotropic_cir_poles:
    ! Computes the split CIR poles for an anisotropic waveguide with
    ! frequencies omega_x /= omega_y (anisotropy ratio eta = omega_x / omega_y).
    ! The CIR splits into two distinct poles: a_CIR,x and a_CIR,y
    ! ==========================================================================
    subroutine calc_anisotropic_cir_poles(omega_x_au, omega_y_au, reduced_mass_au, &
                                         a_cir_x, a_cir_y, stat)
        real(dp), intent(in)  :: omega_x_au
        real(dp), intent(in)  :: omega_y_au
        real(dp), intent(in)  :: reduced_mass_au
        real(dp), intent(out) :: a_cir_x
        real(dp), intent(out) :: a_cir_y
        integer, intent(out)  :: stat

        real(dp) :: a_px, a_py, eta

        stat = 0
        if (omega_x_au <= 0.0_dp .or. omega_y_au <= 0.0_dp .or. reduced_mass_au <= 0.0_dp) then
            stat = 1
            a_cir_x = 0.0_dp
            a_cir_y = 0.0_dp
            return
        end if

        a_px = sqrt(HBAR_AU / (reduced_mass_au * omega_x_au))
        a_py = sqrt(HBAR_AU / (reduced_mass_au * omega_y_au))
        eta = omega_x_au / omega_y_au

        ! Bergeman et al. PRL 91, 163201 (2003):
        ! Shifted CIR poles along x and y principal axes
        a_cir_x = a_px / (OLSHANII_C + 0.15_dp * (eta - 1.0_dp))
        a_cir_y = a_py / (OLSHANII_C - 0.15_dp * (eta - 1.0_dp))
    end subroutine calc_anisotropic_cir_poles

    ! ==========================================================================
    ! calc_quasi_2d_scattering:
    ! In quasi-2D confinement with harmonic oscillator length a_z:
    !   a_2D = C_2D * a_z * exp(-sqrt(pi / 2) * a_z / a_s)
    ! Low-energy 2D scattering amplitude:
    !   f_2D(k) = 2*pi / [ ln(2 / (pi * k * a_2D)) + i * pi / 2 ]
    ! ==========================================================================
    pure subroutine calc_quasi_2d_scattering(pl, a_s_3d_au, k_wave_au, &
                                            a_2d_au, f_2d_amp, sigma_2d_au, stat)
        type(planar_2d_t), intent(in)   :: pl
        real(dp), intent(in)            :: a_s_3d_au
        real(dp), intent(in)            :: k_wave_au
        real(dp), intent(out)           :: a_2d_au
        complex(dp), intent(out)        :: f_2d_amp
        real(dp), intent(out)           :: sigma_2d_au
        integer, intent(out)            :: stat

        real(dp) :: exponent, log_term
        complex(dp) :: denom

        stat = 0
        if (a_s_3d_au <= 0.0_dp .or. k_wave_au <= 0.0_dp) then
            stat = 1
            a_2d_au = 0.0_dp
            f_2d_amp = (0.0_dp, 0.0_dp)
            sigma_2d_au = 0.0_dp
            return
        end if

        exponent = - sqrt(PI * 0.5_dp) * (pl%a_z_au / a_s_3d_au)
        ! Protect against numerical underflow
        if (exponent < -60.0_dp) exponent = -60.0_dp

        a_2d_au = PETROV_PREFACTOR_2D * pl%a_z_au * exp(exponent)

        log_term = log(2.0_dp / (PI * k_wave_au * a_2d_au))
        denom = cmplx(log_term, PI * 0.5_dp, kind=dp)

        f_2d_amp = (TWOPI) / denom
        sigma_2d_au = (abs(f_2d_amp)**2) / k_wave_au
    end subroutine calc_quasi_2d_scattering

    ! ==========================================================================
    ! calc_lieb_liniger_parameter:
    ! Evaluates dimensionless interaction parameter gamma_LL = m * g_1D / (hbar^2 * n_1D)
    ! ==========================================================================
    pure function calc_lieb_liniger_parameter(g_1d_au, n_1d_au, mass_au) result(gamma_ll)
        real(dp), intent(in) :: g_1d_au
        real(dp), intent(in) :: n_1d_au
        real(dp), intent(in) :: mass_au
        real(dp) :: gamma_ll

        if (n_1d_au > 0.0_dp .and. mass_au > 0.0_dp) then
            gamma_ll = (mass_au * abs(g_1d_au)) / ((HBAR_AU**2) * n_1d_au)
        else
            gamma_ll = 0.0_dp
        end if
    end function calc_lieb_liniger_parameter

end module mod_confined_scattering
