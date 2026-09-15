! ==============================================================================
! GeneralModule: mod_bicircular_pecd.f90
!
! Bicircular Laser Fields and Photoelectron Circular Dichroism (PECD)
!
! Theoretical Foundations:
!   1. Bicircular laser field synthesis (omega + 2*omega counter/co-rotating
!      elliptical/circular polarization) with dynamical C_N symmetry.
!      E(t) = E1(t) [cos(w1 t + phi1) x_hat + h1 sin(w1 t + phi1) y_hat] +
!             E2(t) [cos(w2 t + phi2) x_hat + h2 sin(w2 t + phi2) y_hat]
!   2. Photoelectron Circular Dichroism (PECD):
!      Differential photoemission of randomly oriented chiral molecules
!      illuminated by circularly polarized light.
!      I(theta, phi) = sigma_tot / (4*pi) * [1 + beta1*P1(cos theta) + beta2*P2(cos theta) + ...]
!      Forward-backward asymmetry: G_PECD = (I_forward - I_backward) / (I_forward + I_backward) = beta1 / 2.
!   3. Four-center chiral tetrahedral potential model with sign-inverting
!      chirality invariant chi = det(R1-R4, R2-R4, R3-R4) * Prod(Zi - Zj).
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_bicircular_pecd
    use mod_constants, only: dp, PI, TWOPI, AU2EV, EV2AU, FS2AU, AU2FS, W_CM2AU
    implicit none
    private

    public :: bicircular_field_t
    public :: chiral_tetrahedral_molecule_t
    public :: init_bicircular_field
    public :: calc_bicircular_field_at_t
    public :: calc_bicircular_trajectory
    public :: calc_dynamical_symmetry_fold
    public :: init_chiral_tetrahedral_molecule
    public :: calc_chirality_measure
    public :: calc_forward_backward_asymmetry
    public :: calc_chiral_beta1_model
    public :: calc_pecd_pad_spectrum
    public :: calc_pecd_energy_resolved

    ! Using W_CM2AU, FS2AU, EV2AU from mod_constants

    !> Configuration of a two-color bicircular laser field
    type :: bicircular_field_t
        real(dp) :: omega1         !< Fundamental angular frequency [a.u.]
        real(dp) :: omega2         !< Harmonic angular frequency [a.u.]
        real(dp) :: e0_1           !< Electric field peak amplitude color 1 [a.u.]
        real(dp) :: e0_2           !< Electric field peak amplitude color 2 [a.u.]
        integer  :: h1             !< Helicity color 1 (+1: LCP, -1: RCP)
        integer  :: h2             !< Helicity color 2 (+1: LCP, -1: RCP)
        real(dp) :: phi1           !< Carrier envelope phase color 1 [rad]
        real(dp) :: phi2           !< Relative phase color 2 [rad]
        real(dp) :: fwhm_au        !< Pulse duration FWHM [a.u.]
        integer  :: envelope_type  !< 1: Gaussian, 2: Cos^2
    end type bicircular_field_t

    !> Four-center chiral tetrahedral molecule model
    type :: chiral_tetrahedral_molecule_t
        real(dp) :: r_pos(3, 4)    !< Cartesian coordinates of 4 centers [a.u.]
        real(dp) :: charges(4)     !< Effective core charges Z_i
        character(len=32) :: name  !< Molecule name / enantiomer tag
    end type chiral_tetrahedral_molecule_t

contains

    !> Initialize bicircular field parameters
    pure subroutine init_bicircular_field(field, omega1_au, r_freq, i1_wcm2, i2_wcm2, &
                                         h1, h2, phi1, phi2, fwhm_fs, envelope_type)
        type(bicircular_field_t), intent(out) :: field
        real(dp), intent(in) :: omega1_au
        real(dp), intent(in) :: r_freq       !< Frequency ratio omega2 / omega1 (e.g. 2.0_dp)
        real(dp), intent(in) :: i1_wcm2
        real(dp), intent(in) :: i2_wcm2
        integer,  intent(in) :: h1, h2
        real(dp), intent(in) :: phi1, phi2
        real(dp), intent(in) :: fwhm_fs
        integer,  intent(in) :: envelope_type

        field%omega1 = omega1_au
        field%omega2 = omega1_au * r_freq
        ! Circular polarization amplitude: E0 = sqrt(I / (epsilon0 * c)) in a.u. E0 = sqrt(I_au)
        field%e0_1 = sqrt(max(0.0_dp, i1_wcm2 * W_CM2AU))
        field%e0_2 = sqrt(max(0.0_dp, i2_wcm2 * W_CM2AU))
        field%h1 = h1
        field%h2 = h2
        field%phi1 = phi1
        field%phi2 = phi2
        field%fwhm_au = fwhm_fs * FS2AU
        field%envelope_type = envelope_type
    end subroutine init_bicircular_field

    !> Calculate instantaneous electric field E and vector potential A at time t
    pure subroutine calc_bicircular_field_at_t(field, t_au, ex, ey, ax, ay)
        type(bicircular_field_t), intent(in) :: field
        real(dp), intent(in)  :: t_au
        real(dp), intent(out) :: ex, ey, ax, ay

        real(dp) :: env, sigma, t_total
        real(dp) :: arg1, arg2
        real(dp) :: ex1, ey1, ex2, ey2

        ! Temporal envelope f(t) centered at t = 0
        if (field%envelope_type == 2) then
            ! cos^2 envelope from -T_total/2 to +T_total/2 where T_total = 2 * FWHM
            t_total = 2.0_dp * field%fwhm_au
            if (abs(t_au) <= 0.5_dp * t_total) then
                env = cos(PI * t_au / t_total)**2
            else
                env = 0.0_dp
            end if
        else
            ! Standard Gaussian envelope: FWHM = 2 * sqrt(2*ln2) * sigma
            sigma = field%fwhm_au / (2.0_dp * sqrt(2.0_dp * log(2.0_dp)))
            env = exp(-0.5_dp * (t_au / sigma)**2)
        end if

        arg1 = field%omega1 * t_au + field%phi1
        arg2 = field%omega2 * t_au + field%phi2

        ! Field component 1 (omega1, circular with helicity h1)
        ex1 = field%e0_1 * env * cos(arg1)
        ey1 = field%e0_1 * env * real(field%h1, dp) * sin(arg1)

        ! Field component 2 (omega2, circular with helicity h2)
        ex2 = field%e0_2 * env * cos(arg2)
        ey2 = field%e0_2 * env * real(field%h2, dp) * sin(arg2)

        ex = ex1 + ex2
        ey = ey1 + ey2

        ! Vector potential in dipole approximation A(t) = - \int^t E(t') dt'
        ! Under slowly varying envelope approximation:
        ! A1x = - (E0_1 / omega1) * env * sin(arg1)
        ! A1y = + (E0_1 / omega1) * env * h1 * cos(arg1)
        if (field%omega1 > 1.0e-12_dp .and. field%omega2 > 1.0e-12_dp) then
            ax = - (field%e0_1 / field%omega1) * env * sin(arg1) &
                 - (field%e0_2 / field%omega2) * env * sin(arg2)
            ay =   (field%e0_1 / field%omega1) * env * real(field%h1, dp) * cos(arg1) &
                 + (field%e0_2 / field%omega2) * env * real(field%h2, dp) * cos(arg2)
        else
            ax = 0.0_dp
            ay = 0.0_dp
        end if
    end subroutine calc_bicircular_field_at_t

    !> Generate parametric trajectory for Lissajous figure analysis
    pure subroutine calc_bicircular_trajectory(field, n_pts, t_span_au, t_arr, &
                                              ex_arr, ey_arr, ax_arr, ay_arr)
        type(bicircular_field_t), intent(in) :: field
        integer,  intent(in)  :: n_pts
        real(dp), intent(in)  :: t_span_au
        real(dp), intent(out) :: t_arr(n_pts)
        real(dp), intent(out) :: ex_arr(n_pts), ey_arr(n_pts)
        real(dp), intent(out) :: ax_arr(n_pts), ay_arr(n_pts)

        integer :: i
        real(dp) :: dt, t_start

        if (n_pts <= 1) return
        dt = t_span_au / real(n_pts - 1, dp)
        t_start = -0.5_dp * t_span_au

        do i = 1, n_pts
            t_arr(i) = t_start + real(i - 1, dp) * dt
            call calc_bicircular_field_at_t(field, t_arr(i), &
                                            ex_arr(i), ey_arr(i), ax_arr(i), ay_arr(i))
        end do
    end subroutine calc_bicircular_trajectory

    !> Calculate the discrete rotational symmetry fold of the bicircular field
    !> For counter-rotating omega1 and omega2 = r*omega1, symmetry is C_{r + 1}
    !> (e.g. omega + 2*omega counter-rotating gives 3-fold trefoil C3).
    pure function calc_dynamical_symmetry_fold(h1, h2, freq_ratio) result(n_fold)
        integer, intent(in) :: h1, h2
        integer, intent(in) :: freq_ratio  !< Integer frequency ratio (e.g. 2 for 2*omega)
        integer :: n_fold

        if (h1 * h2 < 0) then
            ! Counter-rotating: n_fold = freq_ratio + 1
            n_fold = freq_ratio + 1
        else
            ! Co-rotating: n_fold = abs(freq_ratio - 1)
            n_fold = abs(freq_ratio - 1)
            if (n_fold == 0) n_fold = 1
        end if
    end function calc_dynamical_symmetry_fold

    !> Initialize a prototypical four-center tetrahedral chiral molecule (R or S)
    pure subroutine init_chiral_tetrahedral_molecule(mol, enantiomer_type)
        type(chiral_tetrahedral_molecule_t), intent(out) :: mol
        character(len=*), intent(in) :: enantiomer_type  !< "R", "S", or "ACHIRAL"

        ! Canonical regular tetrahedron with central atom at origin or four distinct vertices
        ! Center 1: (0, 0, 1)
        ! Center 2: (2*sqrt(2)/3, 0, -1/3)
        ! Center 3: (-sqrt(2)/3, sqrt(6)/3, -1/3)
        ! Center 4: (-sqrt(2)/3, -sqrt(6)/3, -1/3)
        real(dp), parameter :: s2 = 1.4142135623730950_dp
        real(dp), parameter :: s6 = 2.4494897427831780_dp

        mol%r_pos(1, 1) = 0.0_dp
        mol%r_pos(2, 1) = 0.0_dp
        mol%r_pos(3, 1) = 1.5_dp

        mol%r_pos(1, 2) = (2.0_dp * s2 / 3.0_dp) * 1.5_dp
        mol%r_pos(2, 2) = 0.0_dp
        mol%r_pos(3, 2) = (-1.0_dp / 3.0_dp) * 1.5_dp

        mol%r_pos(1, 3) = (-s2 / 3.0_dp) * 1.5_dp
        mol%r_pos(2, 3) = (s6 / 3.0_dp) * 1.5_dp
        mol%r_pos(3, 3) = (-1.0_dp / 3.0_dp) * 1.5_dp

        mol%r_pos(1, 4) = (-s2 / 3.0_dp) * 1.5_dp
        mol%r_pos(2, 4) = (-s6 / 3.0_dp) * 1.5_dp
        mol%r_pos(3, 4) = (-1.0_dp / 3.0_dp) * 1.5_dp

        ! Core charges breaking symmetry: Z1=1.0, Z2=2.0, Z3=3.0, Z4=4.0 for R-enantiomer
        if (trim(enantiomer_type) == "R" .or. trim(enantiomer_type) == "r") then
            mol%charges = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
            mol%name = "Tetrahedral_Chiral_R"
        else if (trim(enantiomer_type) == "S" .or. trim(enantiomer_type) == "s") then
            ! Swap charges 2 and 3 to construct the enantiomeric mirror image
            mol%charges = [1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp]
            mol%name = "Tetrahedral_Chiral_S"
        else
            ! Achiral case: identical charges Z2 = Z3
            mol%charges = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp]
            mol%name = "Tetrahedral_Achiral"
        end if
    end subroutine init_chiral_tetrahedral_molecule

    !> Calculate pseudoscalar chirality measure chi
    !> chi = [(R1 - R4) x (R2 - R4)] . (R3 - R4) * Prod_{i < j} (Z_i - Z_j)
    pure function calc_chirality_measure(mol) result(chi)
        type(chiral_tetrahedral_molecule_t), intent(in) :: mol
        real(dp) :: chi

        real(dp) :: v1(3), v2(3), v3(3), cross(3)
        real(dp) :: det_vol, charge_prod
        integer :: i, j

        v1 = mol%r_pos(:, 1) - mol%r_pos(:, 4)
        v2 = mol%r_pos(:, 2) - mol%r_pos(:, 4)
        v3 = mol%r_pos(:, 3) - mol%r_pos(:, 4)

        ! Cross product v1 x v2
        cross(1) = v1(2) * v2(3) - v1(3) * v2(2)
        cross(2) = v1(3) * v2(1) - v1(1) * v2(3)
        cross(3) = v1(1) * v2(2) - v1(2) * v2(1)

        det_vol = dot_product(cross, v3)

        charge_prod = 1.0_dp
        do i = 1, 3
            do j = i + 1, 4
                charge_prod = charge_prod * (mol%charges(i) - mol%charges(j))
            end do
        end do

        chi = det_vol * charge_prod
    end function calc_chirality_measure

    !> Calculate forward-backward asymmetry ratio G_PECD = (I_F - I_B) / (I_F + I_B)
    !> From Ritchie PAD expansion: G_PECD = beta1 / 2.
    pure function calc_forward_backward_asymmetry(beta1) result(g_pecd)
        real(dp), intent(in) :: beta1
        real(dp) :: g_pecd
        g_pecd = 0.5_dp * beta1
    end function calc_forward_backward_asymmetry

    !> Compute first-order chiral asymmetry parameter beta1(E) based on
    !> tetrahedral potential chiral interference
    pure function calc_chiral_beta1_model(mol, energy_ev, photon_energy_ev) result(beta1)
        type(chiral_tetrahedral_molecule_t), intent(in) :: mol
        real(dp), intent(in) :: energy_ev
        real(dp), intent(in) :: photon_energy_ev
        real(dp) :: beta1

        real(dp) :: chi, k_wave, phase_shift, damping
        real(dp) :: scale_factor

        chi = calc_chirality_measure(mol)
        if (abs(chi) < 1.0e-14_dp .or. energy_ev <= 0.0_dp) then
            beta1 = 0.0_dp
            return
        end if

        ! Photoelectron wavenumber k = sqrt(2 * E_au)
        k_wave = sqrt(2.0_dp * max(0.01_dp, energy_ev * EV2AU))

        ! Chiral interference term: beta1 oscillates with k * R_char and decays as 1 / (1 + (k/k0)^4)
        phase_shift = 1.5_dp * k_wave
        damping = 1.0_dp / (1.0_dp + (energy_ev / 30.0_dp)**2)

        ! Normalized chirality sign and strength
        scale_factor = sign(1.0_dp, chi) * min(0.15_dp, abs(chi) * 1.0e-3_dp)

        beta1 = scale_factor * sin(phase_shift) * damping * (photon_energy_ev / 10.0_dp)
    end function calc_chiral_beta1_model

    !> Generate 2D photoelectron angular distribution I(theta, phi)
    pure subroutine calc_pecd_pad_spectrum(beta1, beta2, gamma33, n_theta, n_phi, &
                                          theta_grid, phi_grid, pad_2d)
        real(dp), intent(in)  :: beta1, beta2, gamma33
        integer,  intent(in)  :: n_theta, n_phi
        real(dp), intent(out) :: theta_grid(n_theta), phi_grid(n_phi)
        real(dp), intent(out) :: pad_2d(n_theta, n_phi)

        integer :: it, ip
        real(dp) :: th, ph, cos_th, sin_th
        real(dp) :: p1, p2, p33, norm_factor

        norm_factor = 1.0_dp / (4.0_dp * PI)

        do it = 1, n_theta
            if (n_theta > 1) then
                th = (real(it - 1, dp) / real(n_theta - 1, dp)) * PI
            else
                th = 0.5_dp * PI
            end if
            theta_grid(it) = th
            cos_th = cos(th)
            sin_th = sin(th)

            p1 = cos_th
            p2 = 0.5_dp * (3.0_dp * cos_th**2 - 1.0_dp)
            ! Associated Legendre polynomial P_3^3(cos theta) = 15 * sin^3(theta)
            p33 = 15.0_dp * (sin_th**3)

            do ip = 1, n_phi
                if (n_phi > 1) then
                    ph = (real(ip - 1, dp) / real(n_phi - 1, dp)) * TWOPI
                else
                    ph = 0.0_dp
                end if
                phi_grid(ip) = ph

                ! PAD with chiral P1 asymmetry and 3-fold bicircular modulation cos(3*phi)
                pad_2d(it, ip) = norm_factor * max(0.0_dp, &
                    1.0_dp + beta1 * p1 + beta2 * p2 + gamma33 * p33 * cos(3.0_dp * ph))
            end do
        end do
    end subroutine calc_pecd_pad_spectrum

    !> Compute energy-resolved PECD spectrum and forward-backward asymmetry
    pure subroutine calc_pecd_energy_resolved(mol, e_min_ev, e_max_ev, n_e, &
                                             e_grid, beta1_grid, g_pecd_grid)
        type(chiral_tetrahedral_molecule_t), intent(in) :: mol
        real(dp), intent(in)  :: e_min_ev, e_max_ev
        integer,  intent(in)  :: n_e
        real(dp), intent(out) :: e_grid(n_e)
        real(dp), intent(out) :: beta1_grid(n_e)
        real(dp), intent(out) :: g_pecd_grid(n_e)

        integer :: ie
        real(dp) :: de, eng, b1

        if (n_e <= 1) return
        de = (e_max_ev - e_min_ev) / real(n_e - 1, dp)

        do ie = 1, n_e
            eng = e_min_ev + real(ie - 1, dp) * de
            e_grid(ie) = eng
            b1 = calc_chiral_beta1_model(mol, eng, photon_energy_ev=10.0_dp)
            beta1_grid(ie) = b1
            g_pecd_grid(ie) = calc_forward_backward_asymmetry(b1)
        end do
    end subroutine calc_pecd_energy_resolved

end module mod_bicircular_pecd
