! ==============================================================================
! GeneralModule: mod_triatomic_geometry.f90
!
! Triatomic Molecular Geometry, Jacobi Coordinates, LEPS Reactive Surfaces, &
! Conical Intersection Geometric Phase
!
! Theoretical Foundations:
!   - J. Z. H. Zhang, Theory and Application of Quantum Molecular Dynamics (1999)
!   - S. Sato, J. Chem. Phys. 23, 592 (1955) [London-Eyring-Polanyi-Sato potential]
!   - C. A. Mead & D. G. Truhlar, J. Chem. Phys. 70, 2284 (1979) [Conical intersections]
!   - M. V. Berry, Proc. R. Soc. Lond. A 392, 45 (1984) [Geometric phase]
!
! Key Capabilities:
!   1. Transformation between internuclear distances (r12, r23, r31) and Jacobi (r, R, gamma)
!   2. Analytical London-Eyring-Polanyi-Sato (LEPS) reactive potential energy surface
!   3. Linear E x e / C2v conical intersection model with diabatic coupling
!   4. Berry geometric phase line integral around conical intersection singularity (Phi = pi)
!   5. Minimum energy reaction path (MEP) gradient and transition state barrier
!
! Standard: Fortran 2008
! ==============================================================================

module mod_triatomic_geometry
    use mod_constants, only: dp, PI, TWOPI, HBAR
    implicit none
    private

    ! Public derived types
    public :: jacobi_coord_t, internuclear_dist_t, leps_param_t, conical_intersection_t

    ! Public procedures
    public :: jacobi_to_internuclear
    public :: internuclear_to_jacobi
    public :: calc_leps_potential
    public :: init_default_h3_leps
    public :: calc_conical_intersection_adiabats
    public :: calc_berry_phase_around_ci

    ! Jacobi coordinates for A + BC system
    type :: jacobi_coord_t
        real(dp) :: r_diatom       ! BC internuclear distance r
        real(dp) :: r_atom_diatom  ! Distance R from A to center-of-mass of BC
        real(dp) :: gamma_rad      ! Jacobi angle gamma = arccos(r . R) in radians
        real(dp) :: mass_a         ! Mass of atom A in a.u.
        real(dp) :: mass_b         ! Mass of atom B in a.u.
        real(dp) :: mass_c         ! Mass of atom C in a.u.
    end type jacobi_coord_t

    ! Internuclear distances
    type :: internuclear_dist_t
        real(dp) :: r12            ! Distance between atom 1 and 2
        real(dp) :: r23            ! Distance between atom 2 and 3
        real(dp) :: r31            ! Distance between atom 3 and 1
    end type internuclear_dist_t

    ! LEPS surface parameters (for each of the 3 pairs: 12, 23, 31)
    type :: leps_param_t
        real(dp) :: d_e(3)         ! Well depths D_e in a.u.
        real(dp) :: r_e(3)         ! Equilibrium distances r_e in a.u.
        real(dp) :: beta(3)        ! Morse range parameters beta in a.u.^-1
        real(dp) :: sato_delta(3)  ! Sato empirical parameters Delta (typically 0.05 ~ 0.20)
    end type leps_param_t

    ! 2D Conical Intersection parameters
    type :: conical_intersection_t
        real(dp) :: x_ci           ! Intersection location x_0
        real(dp) :: y_ci           ! Intersection location y_0
        real(dp) :: kappa_tuning   ! Tuning coupling slope kappa_x
        real(dp) :: lambda_coupl   ! Non-adiabatic coupling slope lambda_y
        real(dp) :: e_ci           ! Degeneracy energy E_CI
    end type conical_intersection_t

contains

    ! ==========================================================================
    ! jacobi_to_internuclear:
    ! Converts Jacobi coordinates (r, R, gamma) to internuclear distances (r12, r23, r31)
    ! Atom 1 is A, atom 2 is B, atom 3 is C.
    ! r = distance(B, C) = r23.
    ! R = distance from A to COM(BC).
    ! ==========================================================================
    pure subroutine jacobi_to_internuclear(jac, dist)
        type(jacobi_coord_t), intent(in)     :: jac
        type(internuclear_dist_t), intent(out) :: dist

        real(dp) :: d_b, d_c, r, r_com, cos_g

        r = jac%r_diatom
        r_com = jac%r_atom_diatom
        cos_g = cos(jac%gamma_rad)

        ! Distance from COM(BC) to atom B and C:
        ! d_b = m_c / (m_b + m_c) * r
        ! d_c = m_b / (m_b + m_c) * r
        d_b = (jac%mass_c / (jac%mass_b + jac%mass_c)) * r
        d_c = (jac%mass_b / (jac%mass_b + jac%mass_c)) * r

        ! r23 is diatom distance r
        dist%r23 = r

        ! r12 = distance(A, B) = sqrt(R^2 + d_b^2 - 2*R*d_b*cos(gamma))
        dist%r12 = sqrt(max(1.0e-12_dp, r_com**2 + d_b**2 - 2.0_dp * r_com * d_b * cos_g))

        ! r31 = distance(A, C) = sqrt(R^2 + d_c^2 + 2*R*d_c*cos(gamma))
        dist%r31 = sqrt(max(1.0e-12_dp, r_com**2 + d_c**2 + 2.0_dp * r_com * d_c * cos_g))
    end subroutine jacobi_to_internuclear

    ! ==========================================================================
    ! internuclear_to_jacobi:
    ! Converts internuclear distances (r12, r23, r31) back to Jacobi coordinates
    ! ==========================================================================
    pure subroutine internuclear_to_jacobi(dist, mass_a, mass_b, mass_c, jac)
        type(internuclear_dist_t), intent(in) :: dist
        real(dp), intent(in)                  :: mass_a, mass_b, mass_c
        type(jacobi_coord_t), intent(out)     :: jac

        real(dp) :: d_b, d_c, r23_sq, r12_sq, r31_sq, r_com_sq, cos_g

        jac%mass_a = mass_a
        jac%mass_b = mass_b
        jac%mass_c = mass_c
        jac%r_diatom = dist%r23

        d_b = (mass_c / (mass_b + mass_c)) * dist%r23
        d_c = (mass_b / (mass_b + mass_c)) * dist%r23

        r23_sq = dist%r23**2
        r12_sq = dist%r12**2
        r31_sq = dist%r31**2

        ! By Stewart's theorem: R^2 = (m_b*r12^2 + m_c*r31^2)/(m_b+m_c) - d_b*d_c
        r_com_sq = (mass_b * r12_sq + mass_c * r31_sq) / (mass_b + mass_c) - d_b * d_c
        jac%r_atom_diatom = sqrt(max(1.0e-12_dp, r_com_sq))

        ! Jacobi angle gamma
        if (jac%r_atom_diatom > 1.0e-8_dp .and. dist%r23 > 1.0e-8_dp) then
            cos_g = (r_com_sq + d_b**2 - r12_sq) / (2.0_dp * jac%r_atom_diatom * d_b)
            if (cos_g > 1.0_dp) cos_g = 1.0_dp
            if (cos_g < -1.0_dp) cos_g = -1.0_dp
            jac%gamma_rad = acos(cos_g)
        else
            jac%gamma_rad = 0.0_dp
        end if
    end subroutine internuclear_to_jacobi

    ! ==========================================================================
    ! calc_leps_potential:
    ! Computes the London-Eyring-Polanyi-Sato (LEPS) triatomic potential energy:
    !   V = Q1 + Q2 + Q3 - sqrt(1/2 * [ (J1 - J2)^2 + (J2 - J3)^2 + (J3 - J1)^2 ])
    ! where Q_i and J_i are Coulomb and exchange integrals from Morse curves.
    ! ==========================================================================
    pure function calc_leps_potential(dist, leps) result(v_pot)
        type(internuclear_dist_t), intent(in) :: dist
        type(leps_param_t), intent(in)        :: leps
        real(dp) :: v_pot

        real(dp) :: r(3), d_e, r_e, b, s, x
        real(dp) :: q(3), j(3), diff1, diff2, diff3
        integer  :: k

        r(1) = dist%r12
        r(2) = dist%r23
        r(3) = dist%r31

        do k = 1, 3
            d_e = leps%d_e(k)
            r_e = leps%r_e(k)
            b   = leps%beta(k)
            s   = leps%sato_delta(k)

            x = exp(-b * (r(k) - r_e))

            ! Standard London-Eyring-Polanyi-Sato (LEPS) formulation:
            ! Coulomb integral Q_k and exchange integral J_k
            q(k) = (d_e / (4.0_dp * (1.0_dp + s))) * &
                   ((3.0_dp + s) * x**2 - (2.0_dp + 6.0_dp * s) * x)
            j(k) = (d_e / (4.0_dp * (1.0_dp + s))) * &
                   ((1.0_dp + 3.0_dp * s) * x**2 - (6.0_dp + 2.0_dp * s) * x)
        end do

        diff1 = j(1) - j(2)
        diff2 = j(2) - j(3)
        diff3 = j(3) - j(1)

        v_pot = q(1) + q(2) + q(3) - sqrt(0.5_dp * (diff1**2 + diff2**2 + diff3**2))
    end function calc_leps_potential

    ! ==========================================================================
    ! init_default_h3_leps:
    ! Initializes parameters for the standard benchmark H + H_2 -> H_3 -> H_2 + H LEPS surface
    ! ==========================================================================
    pure subroutine init_default_h3_leps(leps)
        type(leps_param_t), intent(out) :: leps
        integer :: k

        do k = 1, 3
            leps%d_e(k) = 0.1744_dp        ! H2 dissociation energy ~ 4.75 eV = 0.1744 a.u.
            leps%r_e(k) = 1.401_dp         ! H2 equilibrium bond length ~ 0.741 A = 1.401 a.u.
            leps%beta(k) = 1.044_dp        ! Morse parameter in a.u.^-1
            leps%sato_delta(k) = 0.10_dp   ! Standard Sato parameter
        end do
    end subroutine init_default_h3_leps

    ! ==========================================================================
    ! calc_conical_intersection_adiabats:
    ! Solves adiabatic energies for a conical intersection (linear E x e model):
    !   H_diab = [ E_CI + kappa*x,   lambda*y ]
    !            [ lambda*y,        E_CI - kappa*x ]
    !   E_+/- = E_CI +/- sqrt( (kappa*x)^2 + (lambda*y)^2 )
    ! ==========================================================================
    pure subroutine calc_conical_intersection_adiabats(ci, x, y, e_lower, e_upper, gap)
        type(conical_intersection_t), intent(in) :: ci
        real(dp), intent(in)                     :: x, y
        real(dp), intent(out)                    :: e_lower, e_upper, gap

        real(dp) :: dx, dy, delta_rad

        dx = x - ci%x_ci
        dy = y - ci%y_ci

        delta_rad = sqrt((ci%kappa_tuning * dx)**2 + (ci%lambda_coupl * dy)**2)

        e_lower = ci%e_ci - delta_rad
        e_upper = ci%e_ci + delta_rad
        gap = 2.0_dp * delta_rad
    end subroutine calc_conical_intersection_adiabats

    ! ==========================================================================
    ! calc_berry_phase_around_ci:
    ! Integrates the non-adiabatic vector potential A(R) = <1|grad|2> around a closed
    ! circle of radius R centered at (x_ci, y_ci) to compute the geometric Berry phase:
    !   Phi_B = \oint A . dR = pi (Exact topological result for conical intersection)
    ! ==========================================================================
    function calc_berry_phase_around_ci(ci, radius, n_steps) result(berry_phase)
        type(conical_intersection_t), intent(in) :: ci
        real(dp), intent(in)                     :: radius
        integer, intent(in)                      :: n_steps
        real(dp) :: berry_phase

        integer :: i
        real(dp) :: d_phi, phi, dx, dy, mixing_angle, last_angle, d_angle, total_angle

        d_phi = TWOPI / real(n_steps, dp)
        total_angle = 0.0_dp
        last_angle = atan2(ci%lambda_coupl * sin(0.0_dp), ci%kappa_tuning * cos(0.0_dp))

        ! Line integral of polar angle around CI singularity
        ! For linear CI, wave function mixing angle theta_mix = theta / 2
        do i = 1, n_steps
            phi = real(i, dp) * d_phi
            dx = radius * cos(phi)
            dy = radius * sin(phi)

            mixing_angle = atan2(ci%lambda_coupl * dy, ci%kappa_tuning * dx)
            d_angle = mixing_angle - last_angle
            ! Unwind 2pi branches
            if (d_angle > PI) d_angle = d_angle - TWOPI
            if (d_angle < -PI) d_angle = d_angle + TWOPI
            total_angle = total_angle + d_angle
            last_angle = mixing_angle
        end do

        ! Closed loop around singularity accumulates Berry phase |Delta theta| / 2 = pi
        berry_phase = 0.5_dp * abs(total_angle)
    end function calc_berry_phase_around_ci

end module mod_triatomic_geometry
