! ==============================================================================
! GeneralModule: mod_three_body_recombination.f90
!
! Ultracold Three-Body Recombination and Efimov Few-Body Physics
!
! Theoretical Foundations:
!   - V. Efimov, Phys. Lett. B 33, 563 (1970) [Energy levels of 3-boson systems]
!   - E. Braaten & H.-W. Hammer, Phys. Rep. 428, 259 (2006) [Universality in few-body]
!   - B. D. Esry, C. H. Greene, & J. P. Burke, Jr., PRL 83, 1751 (1999) [3-body loss]
!   - P. O. Kraemer et al., Nature 440, 315 (2006) [Observation of Efimov resonances]
!
! Key Capabilities:
!   1. Transcendental root solver for Efimov hyper-radial scale parameter s_0
!   2. Braaten-Hammer universal three-body recombination rate K_3(a) (a > 0 & a < 0)
!   3. Efimov trimer resonance poles, interference minima, and discrete scale factor
!   4. Finite-temperature unitary limit K_3(T) with T^(-2) power-law scaling
!   5. Heteronuclear AAB three-body transcendental root solver with mass ratio beta
!
! Standard: Fortran 2008
! ==============================================================================

module mod_three_body_recombination
    use mod_constants, only: dp, PI, TWOPI, HBAR, KB
    implicit none
    private

    ! Public derived types
    public :: efimov_param_t, three_body_loss_t

    ! Public procedures
    public :: solve_efimov_s0_identical_bosons
    public :: solve_efimov_s0_heteronuclear
    public :: calc_efimov_scale_factor
    public :: calc_three_body_recombination_a_positive
    public :: calc_three_body_recombination_a_negative
    public :: calc_three_body_recombination_universal
    public :: calc_unitary_three_body_loss_temperature
    public :: scan_efimov_recombination_spectrum

    ! Efimov universal parameters type
    type :: efimov_param_t
        real(dp) :: s0             ! Hyper-radial transcendental root (s0 ~ 1.00624)
        real(dp) :: scale_factor   ! Discrete scaling factor lambda = exp(pi / s0) ~ 22.69
        real(dp) :: a_star         ! Three-body parameter a_* (a > 0 trimer at threshold)
        real(dp) :: a_minus        ! Efimov resonance pole position a_- (a < 0 trimer)
        real(dp) :: eta_star       ! Inelasticity parameter eta_* (decay to deep states)
    end type efimov_param_t

    ! Recombination result type
    type :: three_body_loss_t
        real(dp) :: k3_au          ! Three-body rate coefficient in a.u. (a.u.^4 / a.u.)
        real(dp) :: k3_si          ! Three-body rate coefficient in cm^6 / s
        real(dp) :: phase_arg      ! Efimov phase argument s_0 * ln(|a| / a_*)
        logical  :: is_resonance   ! True if near resonance pole
        logical  :: is_minimum     ! True if near interference minimum
    end type three_body_loss_t

    ! Universal dimensionless prefactors (Braaten & Hammer 2006)
    ! C_+ = 128 * pi^2 * (4*pi - 3*sqrt(3)) ~ 67.1160358
    real(dp), parameter :: C_PLUS_PREFACTOR = 67.11603581787_dp
    ! C_- ~ 4590.0
    real(dp), parameter :: C_MINUS_PREFACTOR = 4590.0_dp
    ! Atomic unit of three-body recombination rate (a.u. length^4 * time / mass)
    ! 1 a.u. of K_3 = a_0^4 * (hbar / m_e) ~ 6.57968e-10 cm^6 / s for electron mass
    real(dp), parameter :: AU_LENGTH_CM = 5.29177210903e-9_dp
    real(dp), parameter :: AU_TIME_SEC = 2.4188843265857e-17_dp

contains

    ! ==========================================================================
    ! solve_efimov_s0_identical_bosons:
    ! Solves the transcendental Efimov root equation for three identical bosons:
    !   f(s) = s * cosh(s * pi / 2) - (8 / sqrt(3)) * sinh(s * pi / 6) = 0
    ! Expected solution: s_0 = 1.006237825...
    ! ==========================================================================
    pure function solve_efimov_s0_identical_bosons() result(s0)
        real(dp) :: s0
        real(dp) :: s_curr, f_val, df_val, ds
        integer  :: iter

        ! Initial guess
        s_curr = 1.0_dp

        ! Newton-Raphson iteration
        do iter = 1, 30
            f_val = s_curr * cosh(s_curr * PI * 0.5_dp) - &
                    (8.0_dp / sqrt(3.0_dp)) * sinh(s_curr * PI / 6.0_dp)

            df_val = cosh(s_curr * PI * 0.5_dp) + &
                     s_curr * (PI * 0.5_dp) * sinh(s_curr * PI * 0.5_dp) - &
                     (8.0_dp / sqrt(3.0_dp)) * (PI / 6.0_dp) * cosh(s_curr * PI / 6.0_dp)

            ds = -f_val / df_val
            s_curr = s_curr + ds
            if (abs(ds) < 1.0e-13_dp) exit
        end do

        s0 = s_curr
    end function solve_efimov_s0_identical_bosons

    ! ==========================================================================
    ! solve_efimov_s0_heteronuclear:
    ! Solves the transcendental root for AAB heteronuclear three-body system
    ! with mass ratio beta = m_A / m_B (two identical bosons A and one particle B).
    ! Transcendental equation:
    !   s * cosh(s * pi / 2) - (2 / sin(2*alpha)) * sinh(s * alpha) = 0
    ! where alpha = arcsin(1 / (1 + beta))
    ! ==========================================================================
    pure function solve_efimov_s0_heteronuclear(mass_ratio_beta) result(s0)
        real(dp), intent(in) :: mass_ratio_beta
        real(dp) :: s0
        real(dp) :: alpha, sin_2alpha, s_curr, f_val, df_val, ds
        integer  :: iter

        if (mass_ratio_beta <= 0.0_dp) then
            s0 = 0.0_dp
            return
        end if

        alpha = asin(1.0_dp / (1.0_dp + mass_ratio_beta))
        sin_2alpha = sin(2.0_dp * alpha)

        if (mass_ratio_beta < 0.5_dp) then
            s_curr = 2.0_dp
        else
            s_curr = 1.0_dp
        end if

        do iter = 1, 40
            f_val = s_curr * cosh(s_curr * PI * 0.5_dp) - &
                    (2.0_dp / sin_2alpha) * sinh(s_curr * alpha)

            df_val = cosh(s_curr * PI * 0.5_dp) + &
                     s_curr * (PI * 0.5_dp) * sinh(s_curr * PI * 0.5_dp) - &
                     (2.0_dp / sin_2alpha) * alpha * cosh(s_curr * alpha)

            ds = -f_val / df_val
            if (abs(ds) > 1.0_dp) ds = sign(1.0_dp, ds)
            s_curr = s_curr + ds
            if (s_curr < 0.1_dp) s_curr = 0.5_dp
            if (abs(ds) < 1.0e-13_dp) exit
        end do

        s0 = s_curr
    end function solve_efimov_s0_heteronuclear

    ! ==========================================================================
    ! calc_efimov_scale_factor:
    ! Discrete scaling factor lambda = exp(pi / s_0).
    ! For 3 identical bosons: lambda ~ exp(pi / 1.00624) ~ 22.694
    ! ==========================================================================
    pure function calc_efimov_scale_factor(s0) result(lambda_scale)
        real(dp), intent(in) :: s0
        real(dp) :: lambda_scale

        if (s0 > 0.0_dp) then
            lambda_scale = exp(PI / s0)
        else
            lambda_scale = 0.0_dp
        end if
    end function calc_efimov_scale_factor

    ! ==========================================================================
    ! calc_three_body_recombination_a_positive:
    ! Braaten-Hammer universal formula for a > 0:
    !   K_3 = (C_+ * hbar / m) * a^4 * [ sin^2(s_0 * ln(a / a_*)) + sinh^2(eta_*) ]
    ! ==========================================================================
    pure subroutine calc_three_body_recombination_a_positive(a_scat, mass_atom, &
                                                           param, res, stat)
        real(dp), intent(in)              :: a_scat       ! Scattering length (a.u., > 0)
        real(dp), intent(in)              :: mass_atom    ! Atom mass (a.u.)
        type(efimov_param_t), intent(in)  :: param        ! Efimov parameters
        type(three_body_loss_t), intent(out) :: res       ! Loss results
        integer, intent(out)              :: stat

        real(dp) :: phase, sin_val, sinh_val, factor, k3_au_val

        stat = 0
        if (a_scat <= 0.0_dp .or. mass_atom <= 0.0_dp .or. param%a_star <= 0.0_dp) then
            stat = 1
            res%k3_au = 0.0_dp
            res%k3_si = 0.0_dp
            return
        end if

        phase = param%s0 * log(a_scat / param%a_star)
        res%phase_arg = phase

        sin_val = sin(phase)
        sinh_val = sinh(param%eta_star)

        factor = sin_val**2 + sinh_val**2

        ! K_3 in atomic units: [length]^4 * [hbar / m]
        k3_au_val = (C_PLUS_PREFACTOR * HBAR / mass_atom) * (a_scat**4) * factor
        res%k3_au = k3_au_val

        ! Convert to cm^6 / s
        res%k3_si = k3_au_val * ((AU_LENGTH_CM)**6 / AU_TIME_SEC)

        res%is_resonance = .false.
        res%is_minimum = (abs(sin_val) < 0.05_dp)
    end subroutine calc_three_body_recombination_a_positive

    ! ==========================================================================
    ! calc_three_body_recombination_a_negative:
    ! Braaten-Hammer universal formula for a < 0:
    !   K_3 = (C_- * hbar / m) * a^4 * sinh(2*eta_*) / [ sin^2(s_0*ln(|a|/a_-)) + sinh^2(eta_*) ]
    ! ==========================================================================
    pure subroutine calc_three_body_recombination_a_negative(a_scat, mass_atom, &
                                                           param, res, stat)
        real(dp), intent(in)              :: a_scat       ! Scattering length (a.u., < 0)
        real(dp), intent(in)              :: mass_atom    ! Atom mass (a.u.)
        type(efimov_param_t), intent(in)  :: param        ! Efimov parameters
        type(three_body_loss_t), intent(out) :: res       ! Loss results
        integer, intent(out)              :: stat

        real(dp) :: abs_a, a_minus_val, phase, sin_val, denom, k3_au_val

        stat = 0
        abs_a = abs(a_scat)
        if (abs_a <= 0.0_dp .or. mass_atom <= 0.0_dp) then
            stat = 1
            res%k3_au = 0.0_dp
            res%k3_si = 0.0_dp
            return
        end if

        ! If a_minus not given, determine from a_star by universal ratio
        if (param%a_minus > 0.0_dp) then
            a_minus_val = param%a_minus
        else
            a_minus_val = param%a_star * exp(-PI / (2.0_dp * param%s0)) / sqrt(2.0_dp)
        end if

        phase = param%s0 * log(abs_a / a_minus_val)
        res%phase_arg = phase

        sin_val = sin(phase)
        denom = sin_val**2 + sinh(param%eta_star)**2

        if (denom < 1.0e-30_dp) denom = 1.0e-30_dp

        k3_au_val = (C_MINUS_PREFACTOR * HBAR / mass_atom) * (abs_a**4) * &
                    sinh(2.0_dp * param%eta_star) / denom

        res%k3_au = k3_au_val
        res%k3_si = k3_au_val * ((AU_LENGTH_CM)**6 / AU_TIME_SEC)

        res%is_resonance = (abs(sin_val) < 0.05_dp)
        res%is_minimum = .false.
    end subroutine calc_three_body_recombination_a_negative

    ! ==========================================================================
    ! calc_three_body_recombination_universal:
    ! Dispatcher for arbitrary scattering length sign a > 0 or a < 0
    ! ==========================================================================
    subroutine calc_three_body_recombination_universal(a_scat, mass_atom, &
                                                       param, res, stat)
        real(dp), intent(in)              :: a_scat
        real(dp), intent(in)              :: mass_atom
        type(efimov_param_t), intent(in)  :: param
        type(three_body_loss_t), intent(out) :: res
        integer, intent(out)              :: stat

        if (a_scat > 0.0_dp) then
            call calc_three_body_recombination_a_positive(a_scat, mass_atom, param, res, stat)
        else if (a_scat < 0.0_dp) then
            call calc_three_body_recombination_a_negative(a_scat, mass_atom, param, res, stat)
        else
            res%k3_au = 0.0_dp
            res%k3_si = 0.0_dp
            res%phase_arg = 0.0_dp
            res%is_resonance = .false.
            res%is_minimum = .false.
            stat = 0
        end if
    end subroutine calc_three_body_recombination_universal

    ! ==========================================================================
    ! calc_unitary_three_body_loss_temperature:
    ! At the unitary limit |a| -> infinity, K_3 is bounded by temperature:
    !   K_3(T) = C_unit * (hbar^5 / (m^3 * (k_B * T)^2))
    ! Theoretical prefactor C_unit ~ 65.9 * (1 - exp(-4*eta_*))
    ! ==========================================================================
    pure subroutine calc_unitary_three_body_loss_temperature(temp_kelvin, mass_atom, &
                                                            eta_star, k3_au, k3_si, stat)
        real(dp), intent(in)  :: temp_kelvin   ! Temperature in Kelvin
        real(dp), intent(in)  :: mass_atom     ! Atom mass in a.u.
        real(dp), intent(in)  :: eta_star      ! Inelasticity parameter
        real(dp), intent(out) :: k3_au         ! Recombination rate in a.u.
        real(dp), intent(out) :: k3_si         ! Recombination rate in cm^6 / s
        integer, intent(out)  :: stat

        real(dp) :: e_th, c_unit

        stat = 0
        if (temp_kelvin <= 0.0_dp .or. mass_atom <= 0.0_dp) then
            stat = 1
            k3_au = 0.0_dp
            k3_si = 0.0_dp
            return
        end if

        e_th = KB * temp_kelvin  ! Thermal energy in a.u.
        c_unit = 65.9_dp * (1.0_dp - exp(-4.0_dp * eta_star))

        k3_au = c_unit * (HBAR**5) / ((mass_atom**3) * (e_th**2))
        k3_si = k3_au * ((AU_LENGTH_CM)**6 / AU_TIME_SEC)
    end subroutine calc_unitary_three_body_loss_temperature

    ! ==========================================================================
    ! scan_efimov_recombination_spectrum:
    ! Scans three-body recombination rate K_3 over an array of scattering lengths
    ! ==========================================================================
    subroutine scan_efimov_recombination_spectrum(a_grid, n_pts, mass_atom, &
                                                 param, k3_au_arr, k3_si_arr, stat)
        integer, intent(in)               :: n_pts
        real(dp), intent(in)              :: a_grid(n_pts)
        real(dp), intent(in)              :: mass_atom
        type(efimov_param_t), intent(in)  :: param
        real(dp), intent(out)             :: k3_au_arr(n_pts)
        real(dp), intent(out)             :: k3_si_arr(n_pts)
        integer, intent(out)              :: stat

        integer :: i
        type(three_body_loss_t) :: res
        integer :: s

        stat = 0
        do i = 1, n_pts
            call calc_three_body_recombination_universal(a_grid(i), mass_atom, param, res, s)
            if (s /= 0) stat = s
            k3_au_arr(i) = res%k3_au
            k3_si_arr(i) = res%k3_si
        end do
    end subroutine scan_efimov_recombination_spectrum

end module mod_three_body_recombination
