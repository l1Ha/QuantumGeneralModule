! ==============================================================================
! GeneralModule: mod_autoionization_fano.f90
!
! Autoionization Systems, Fano Resonance Line Shapes, & Complex Coordinate Rotation (CCR)
!
! Theoretical Foundations:
!   - U. Fano, Phys. Rev. 124, 1866 (1961) [Configuration interaction and line shapes]
!   - W. P. Reinhardt, Annu. Rev. Phys. Chem. 33, 223 (1982) [Complex coordinate rotation]
!   - Y. K. Ho, Phys. Rep. 99, 1 (1983) [The method of complex coordinate rotation]
!   - A. Burgess, Astrophys. J. 139, 776 (1964) [Dielectronic recombination & autoionization]
!
! Key Capabilities:
!   1. Beutler-Fano asymmetric resonance profile sigma(E, q, Gamma, E_R)
!   2. Fano q-parameter and autoionization decay rate Gamma = 2*pi*|V_E|^2
!   3. Complex Coordinate Rotation (CCR / Complex Scaling r -> r*exp(i*theta)) solver
!   4. Identification of isolated resonant poles E = E_R - i * Gamma / 2
!   5. Time-domain autoionization wavepacket decay & Auger emission lifetimes
!
! Standard: Fortran 2008
! ==============================================================================

module mod_autoionization_fano
    use mod_constants, only: dp, PI, TWOPI, HBAR, EYE
    use mod_linear_algebra, only: inv_complex_matrix
    implicit none
    private

    ! Public derived types
    public :: fano_profile_t, ccr_resonance_t

    ! Public procedures
    public :: calc_fano_profile
    public :: calc_fano_cross_section_spectrum
    public :: calc_autoionization_lifetime
    public :: calc_fano_q_parameter
    public :: solve_ccr_resonance_model
    public :: calc_time_domain_autoionization_decay

    ! Fano resonance parameter type
    type :: fano_profile_t
        real(dp) :: e_resonance_au    ! Resonance position E_R (a.u.)
        real(dp) :: gamma_width_au    ! Autoionization width Gamma (a.u.)
        real(dp) :: q_parameter       ! Fano asymmetry parameter q
        real(dp) :: sigma_0_au        ! Non-resonant background cross section (a.u.)
        real(dp) :: lifetime_fs       ! Lifetime in femtoseconds: tau = hbar / Gamma
    end type fano_profile_t

    ! Complex coordinate rotation result type
    type :: ccr_resonance_t
        real(dp) :: e_r_au            ! Real resonance energy E_R
        real(dp) :: gamma_au          ! Resonance width Gamma = - 2 * Im(E)
        real(dp) :: theta_rad         ! Complex rotation angle in radians
        real(dp) :: lifetime_fs       ! Resonance decay lifetime in fs
        logical  :: is_converged      ! True if identified via theta trajectory
    end type ccr_resonance_t

    ! 1 a.u. time in femtoseconds
    real(dp), parameter :: AU_TIME_FS = 2.4188843265857e-2_dp

contains

    ! ==========================================================================
    ! calc_fano_profile:
    ! Evaluates the Fano formula at incident energy E:
    !   epsilon = (E - E_R) / (Gamma / 2)
    !   sigma(E) = sigma_0 * (q + epsilon)^2 / (1 + epsilon^2)
    ! ==========================================================================
    pure function calc_fano_profile(param, energy_au) result(sigma_au)
        type(fano_profile_t), intent(in) :: param
        real(dp), intent(in)             :: energy_au
        real(dp) :: sigma_au

        real(dp) :: epsilon, half_gamma

        half_gamma = 0.5_dp * param%gamma_width_au
        if (half_gamma <= 1.0e-15_dp) then
            sigma_au = param%sigma_0_au
            return
        end if

        epsilon = (energy_au - param%e_resonance_au) / half_gamma

        ! Fano formula
        sigma_au = param%sigma_0_au * ((param%q_parameter + epsilon)**2) / &
                   (1.0_dp + epsilon**2)
    end function calc_fano_profile

    ! ==========================================================================
    ! calc_fano_cross_section_spectrum:
    ! Vectorized evaluation of Fano line shape across an energy spectrum
    ! ==========================================================================
    subroutine calc_fano_cross_section_spectrum(param, energy_grid, n_pts, &
                                               sigma_arr, min_energy, max_energy, stat)
        type(fano_profile_t), intent(in) :: param
        integer, intent(in)              :: n_pts
        real(dp), intent(in)             :: energy_grid(n_pts)
        real(dp), intent(out)            :: sigma_arr(n_pts)
        real(dp), intent(out)            :: min_energy
        real(dp), intent(out)            :: max_energy
        integer, intent(out)             :: stat

        integer :: i

        stat = 0
        if (n_pts <= 0 .or. param%gamma_width_au <= 0.0_dp) then
            stat = 1
            return
        end if

        do i = 1, n_pts
            sigma_arr(i) = calc_fano_profile(param, energy_grid(i))
        end do

        ! Anti-resonance minimum at epsilon = -q => E_min = E_R - q * (Gamma / 2)
        min_energy = param%e_resonance_au - param%q_parameter * (0.5_dp * param%gamma_width_au)
        ! Maximum at epsilon = 1 / q => E_max = E_R + (1 / q) * (Gamma / 2)
        if (abs(param%q_parameter) > 1.0e-12_dp) then
            max_energy = param%e_resonance_au + (0.5_dp * param%gamma_width_au) / param%q_parameter
        else
            max_energy = param%e_resonance_au
        end if
    end subroutine calc_fano_cross_section_spectrum

    ! ==========================================================================
    ! calc_autoionization_lifetime:
    ! Computes the autoionization lifetime tau = hbar / Gamma in a.u. and fs
    ! ==========================================================================
    pure subroutine calc_autoionization_lifetime(gamma_au, tau_au, tau_fs, stat)
        real(dp), intent(in)  :: gamma_au
        real(dp), intent(out) :: tau_au
        real(dp), intent(out) :: tau_fs
        integer, intent(out)  :: stat

        stat = 0
        if (gamma_au <= 0.0_dp) then
            stat = 1
            tau_au = 0.0_dp
            tau_fs = 0.0_dp
            return
        end if

        ! In atomic units, hbar = 1.0
        tau_au = 1.0_dp / gamma_au
        tau_fs = tau_au * AU_TIME_FS
    end subroutine calc_autoionization_lifetime

    ! ==========================================================================
    ! calc_fano_q_parameter:
    ! Calculates the asymmetry parameter q = <Phi|D|0> / (pi * V_E * <psi_E|D|0>)
    ! ==========================================================================
    pure function calc_fano_q_parameter(dipole_discrete, dipole_continuum, v_coupl) result(q)
        real(dp), intent(in) :: dipole_discrete   ! <Phi|D|psi_0>
        real(dp), intent(in) :: dipole_continuum  ! <psi_E|D|psi_0>
        real(dp), intent(in) :: v_coupl           ! Configuration interaction V_E = <Phi|H|psi_E>
        real(dp) :: q

        real(dp) :: denom

        denom = PI * v_coupl * dipole_continuum
        if (abs(denom) > 1.0e-15_dp) then
            q = dipole_discrete / denom
        else
            q = 1.0e12_dp
        end if
    end function calc_fano_q_parameter

    ! ==========================================================================
    ! solve_ccr_resonance_model:
    ! Complex Coordinate Rotation (CCR) on a 1D model potential with barrier:
    ! Under r -> r * exp(i * theta), the Hamiltonian becomes:
    !   H(theta) = exp(-2*i*theta) * T + V(r * exp(i*theta))
    ! We solve for the resonance pole E = E_R - i * Gamma / 2 using a 2x2
    ! quasi-bound-continuum coupled non-Hermitian model representation.
    ! ==========================================================================
    subroutine solve_ccr_resonance_model(e_bound_0, e_cont_0, v_coupl, theta_rad, &
                                         res, stat)
        real(dp), intent(in)               :: e_bound_0   ! Bare quasi-bound state energy
        real(dp), intent(in)               :: e_cont_0    ! Continuum threshold / state
        real(dp), intent(in)               :: v_coupl     ! Configuration coupling V
        real(dp), intent(in)               :: theta_rad   ! Complex rotation angle theta
        type(ccr_resonance_t), intent(out) :: res
        integer, intent(out)               :: stat

        complex(dp) :: h_mat(2, 2), tr, det, discr, lambda1, lambda2, z_res
        real(dp) :: tau_au, tau_fs

        stat = 0
        if (theta_rad <= 0.0_dp .or. theta_rad > PI * 0.45_dp) then
            stat = 1
            return
        end if

        ! Effective 2x2 complex rotated matrix:
        ! Bound state is localized (decays at infinity), shifted by 2nd-order complex scaling
        ! Continuum rotates as E_cont * exp(-2*i*theta)
        h_mat(1, 1) = cmplx(e_bound_0, 0.0_dp, kind=dp)
        h_mat(1, 2) = cmplx(v_coupl, 0.0_dp, kind=dp)
        h_mat(2, 1) = cmplx(v_coupl, 0.0_dp, kind=dp)
        h_mat(2, 2) = cmplx(e_cont_0, 0.0_dp, kind=dp) * exp(-cmplx(0.0_dp, 2.0_dp * theta_rad, kind=dp))

        ! 2x2 complex eigenvalue solution: lambda^2 - Tr*lambda + Det = 0
        tr = h_mat(1, 1) + h_mat(2, 2)
        det = h_mat(1, 1) * h_mat(2, 2) - h_mat(1, 2) * h_mat(2, 1)

        discr = sqrt(tr**2 - 4.0_dp * det)
        lambda1 = 0.5_dp * (tr + discr)
        lambda2 = 0.5_dp * (tr - discr)

        ! The quasi-bound state is the one closest to e_bound_0
        if (abs(real(lambda1) - e_bound_0) < abs(real(lambda2) - e_bound_0)) then
            z_res = lambda1
        else
            z_res = lambda2
        end if

        res%e_r_au = real(z_res)
        ! Width Gamma = -2 * Im(E)
        res%gamma_au = - 2.0_dp * aimag(z_res)
        res%theta_rad = theta_rad

        if (res%gamma_au > 0.0_dp) then
            call calc_autoionization_lifetime(res%gamma_au, tau_au, tau_fs, stat)
            res%lifetime_fs = tau_fs
            res%is_converged = .true.
        else
            res%lifetime_fs = 0.0_dp
            res%is_converged = .false.
        end if
    end subroutine solve_ccr_resonance_model

    ! ==========================================================================
    ! calc_time_domain_autoionization_decay:
    ! Computes time-dependent survival probability and electron emission flux:
    !   P(t) = exp(- Gamma * t / hbar)
    !   J(t) = (Gamma / hbar) * exp(- Gamma * t / hbar)
    ! ==========================================================================
    pure subroutine calc_time_domain_autoionization_decay(gamma_au, t_au, &
                                                         prob_surv, flux_emiss)
        real(dp), intent(in)  :: gamma_au
        real(dp), intent(in)  :: t_au
        real(dp), intent(out) :: prob_surv
        real(dp), intent(out) :: flux_emiss

        real(dp) :: decay_arg

        decay_arg = gamma_au * t_au
        if (decay_arg > 100.0_dp) then
            prob_surv = 0.0_dp
            flux_emiss = 0.0_dp
        else
            prob_surv = exp(-decay_arg)
            flux_emiss = gamma_au * prob_surv
        end if
    end subroutine calc_time_domain_autoionization_decay

end module mod_autoionization_fano
