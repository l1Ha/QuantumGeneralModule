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
!   - K. Fukui, J. Phys. Chem. 74, 4161 (1970) [intrinsic reaction coordinate]
!   - K. Ishida, K. Morokuma, A. Komornicki, J. Chem. Phys. 66, 2153 (1977) [IRC numerics]
!
! Key Capabilities:
!   1. Transformation between internuclear distances (r12, r23, r31) and Jacobi (r, R, gamma)
!   2. Analytical London-Eyring-Polanyi-Sato (LEPS) reactive potential energy surface
!      with analytic gradient dV/dr_k
!   3. Linear E x e / C2v conical intersection model with diabatic coupling
!   4. Berry geometric phase line integral around conical intersection singularity (Phi = pi)
!   5. Minimum energy reaction path (MEP) gradient and transition state barrier
!   6. Intrinsic reaction coordinate (IRC): mass-weighted steepest descent in
!      the Wilson G-matrix metric with Gonzalez-Schlegel constrained steps
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
    public :: calc_leps_gradient
    public :: intrinsic_reaction_coordinate
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
    ! calc_leps_gradient:
    ! Analytic gradient of the LEPS surface with respect to the three
    ! internuclear distances, treated as independent: grad(k) = dV/dr_k for
    ! k = 1,2,3 mapping to (r12, r23, r31). With x_k = exp(-b_k (r_k - re_k)),
    ! dx_k/dr_k = -b_k x_k and c_k = de_k / (4 (1 + s_k)):
    !   Q_k = c_k[(3+s_k) x^2 - (2+6 s_k) x],  Q'_k = c_k b_k[(2+6 s_k) x - 2(3+s_k) x^2]
    !   J_k = c_k[(1+3 s_k) x^2 - (6+2 s_k) x], J'_k = c_k b_k[(6+2 s_k) x - 2(1+3 s_k) x^2]
    !   V = sum_k Q_k - sqrt(Delta),  Delta = 1/2 sum (J_i - J_j)^2
    !   dV/dr_k = Q'_k - (3 J_k - J_sum) J'_k / (2 sqrt(Delta)),  J_sum = J1+J2+J3
    ! (d Delta/d J_k = 2 J_k - J_j - J_j' = 3 J_k - J_sum.)
    ! Collinear constraint r31 = r12 + r23: chain rule gives
    !   dV/dr12|collinear = grad(1) + grad(3),  dV/dr23|collinear = grad(2) + grad(3).
    ! ==========================================================================
    pure function calc_leps_gradient(dist, leps) result(grad)
        type(internuclear_dist_t), intent(in) :: dist
        type(leps_param_t), intent(in)        :: leps
        real(dp) :: grad(3)

        real(dp) :: r(3), x, b, s, c, de
        real(dp) :: q_p(3), j(3), j_p(3), j_sum, sq
        integer  :: k

        r(1) = dist%r12
        r(2) = dist%r23
        r(3) = dist%r31

        do k = 1, 3
            de = leps%d_e(k)
            b  = leps%beta(k)
            s  = leps%sato_delta(k)
            c  = de / (4.0_dp * (1.0_dp + s))
            x  = exp(-b * (r(k) - leps%r_e(k)))
            q_p(k) = c * b * ((2.0_dp + 6.0_dp * s) * x - 2.0_dp * (3.0_dp + s) * x**2)
            j(k)   = c * ((1.0_dp + 3.0_dp * s) * x**2 - (6.0_dp + 2.0_dp * s) * x)
            j_p(k) = c * b * ((6.0_dp + 2.0_dp * s) * x - 2.0_dp * (1.0_dp + 3.0_dp * s) * x**2)
        end do

        j_sum = j(1) + j(2) + j(3)
        sq = sqrt(max(0.5_dp * ((j(1) - j(2))**2 + (j(2) - j(3))**2 &
                                + (j(3) - j(1))**2), 1.0e-300_dp))
        do k = 1, 3
            grad(k) = q_p(k) - (3.0_dp * j(k) - j_sum) * j_p(k) / (2.0_dp * sq)
        end do
    end function calc_leps_gradient

    ! ==========================================================================
    ! intrinsic_reaction_coordinate:
    ! Intrinsic reaction coordinate (Fukui 1970) for a COLLINEAR triatomic,
    ! integrated in the mass-weighted metric with Gonzalez-Schlegel (1990)
    ! constrained second-order steps. In bond coordinates (r12, r23) the
    ! velocity-space kinetic metric of the collinear A1-A2-A3 chain is
    !   T = 1/2 v^T G_vel v,
    !   G_vel = 1/M [[ m1(m2+m3),  m1 m3 ],  M = m1+m2+m3,
    !                [ m1 m3,  m3(m1+m2) ]]
    ! whose inverse is exactly the Wilson G-matrix
    !   G = [[ 1/m1+1/m2,   -1/m2  ],
    !        [  -1/m2,   1/m2+1/m3 ]].
    ! The IRC obeys the unit-speed gradient ODE (mass-weighted arc length s):
    !   d r/ds = -G g / sqrt(g^T G g),  so dV/ds = -sqrt(g^T G g) <= 0:
    ! the energy decreases monotonically and the path ends at a valley floor
    ! (reactant/product asymptote), where sqrt(g^T G g) vanishes.
    ! Each step uses the GS construction: pivot r_p = r_n + ds*v_n with the
    ! unit-speed direction v_n = -G g_n / sqrt(g_n^T G g_n); the new point is
    ! searched on the metric-orthogonal bisector r = r_p + xi*t through the
    ! scalar parallelism equation  sin(angle( G g(r), r_{n+1} - r_n )) = 0
    ! (bracketed bisection + secant, ~10 gradient evaluations). Plain
    ! unit-speed Euler steps r_{n+1} = r_n + ds*v_n are only conditionally
    ! stable: with transverse valley curvature k_t the stability window is
    ! ds < 2 sqrt(g^T G g)/k_t, which collapses on the flat valley floor and
    ! produces a growing transverse zigzag (observed at |g| ~ 1.5e-3 for the
    ! H3 LEPS valley, ds = 0.01). The GS constrained step has no such
    ! instability and keeps mirror symmetry to roundoff.
    ! The caller supplies r_start = r_saddle +/- eps * u_react, with u_react
    ! the unstable (imaginary-frequency) normal mode in mass-weighted
    ! coordinates (for symmetric H3: u = (1,-1)/sqrt(2)).
    ! stop_reason: 0 = n_max steps, 1 = gradient converged (valley floor),
    !              2 = left the domain [r_min, r_max], 3 = energy ascent
    !              detected (numerical; decrease ds_step), -1 = invalid input.
    ! ==========================================================================
    subroutine intrinsic_reaction_coordinate(leps, masses, r_start, n_max, ds_step, grad_tol, &
                                             r_min, r_max, s_path, r12_path, r23_path, v_path, &
                                             n_out, stop_reason)
        type(leps_param_t), intent(in) :: leps
        real(dp), intent(in)  :: masses(3)
        real(dp), intent(in)  :: r_start(2)
        integer, intent(in)   :: n_max
        real(dp), intent(in)  :: ds_step, grad_tol, r_min, r_max
        real(dp), intent(out), allocatable :: s_path(:), r12_path(:), r23_path(:), v_path(:)
        integer, intent(out)  :: n_out, stop_reason

        type(internuclear_dist_t) :: dist
        real(dp) :: m1, m2, m3, m_tot
        real(dp) :: g11, g12, g22                 ! Wilson G (inverse velocity metric)
        real(dp) :: h11, h12, h22                 ! velocity metric G_vel
        real(dp) :: r(2), g(2), gg(2), w(2), v_n(2), t_vec(2), u_vec(2)
        real(dp) :: r_p(2), r_new(2), r_eul(2)
        real(dp) :: chord(2), gnrm, gwnrm, xi, xl, xr, fl, fr, fm, tmp
        real(dp) :: s_acc, v_old, v_new, v_eul
        real(dp) :: g3(3)
        real(dp), allocatable :: s_tmp(:), r1_tmp(:), r2_tmp(:), v_tmp(:)
        integer  :: istep, ibis, ibrack

        stop_reason = -1
        n_out = 0
        m1 = masses(1); m2 = masses(2); m3 = masses(3)
        if (min(m1, m2, m3) <= 0.0_dp) return
        if (n_max < 2 .or. ds_step <= 0.0_dp .or. grad_tol <= 0.0_dp) return
        if (r_start(1) <= r_min .or. r_start(1) >= r_max) return
        if (r_start(2) <= r_min .or. r_start(2) >= r_max) return

        m_tot = m1 + m2 + m3
        ! Wilson G-matrix (inverse of the velocity metric), collinear bond coords
        g11 = 1.0_dp / m1 + 1.0_dp / m2
        g12 = -1.0_dp / m2
        g22 = 1.0_dp / m2 + 1.0_dp / m3
        ! velocity metric G_vel (for the unit-speed normalization)
        h11 = m1 * (m2 + m3) / m_tot
        h12 = m1 * m3 / m_tot
        h22 = m3 * (m1 + m2) / m_tot

        allocate(s_tmp(n_max + 1), r1_tmp(n_max + 1), r2_tmp(n_max + 1), v_tmp(n_max + 1))

        r = r_start
        dist%r12 = r(1); dist%r23 = r(2); dist%r31 = r(1) + r(2)
        v_old = calc_leps_potential(dist, leps)
        n_out = 1
        s_tmp(1) = 0.0_dp
        r1_tmp(1) = r(1)
        r2_tmp(1) = r(2)
        v_tmp(1) = v_old
        s_acc = 0.0_dp
        stop_reason = 0

        do istep = 1, n_max
            ! --- gradient and unit-speed direction at the current point --------
            dist%r12 = r(1); dist%r23 = r(2); dist%r31 = r(1) + r(2)
            g3 = calc_leps_gradient(dist, leps)
            ! collinear chain rule: dV/dr12 = grad(1)+grad(3), dV/dr23 = grad(2)+grad(3)
            g(1) = g3(1) + g3(3)
            g(2) = g3(2) + g3(3)

            ! mass-weighted gradient norm sqrt(g^T G g): convergence at valley floor
            gnrm = sqrt(max(g(1)**2 * g11 + 2.0_dp * g(1) * g(2) * g12 &
                            + g(2)**2 * g22, 0.0_dp))
            if (gnrm < grad_tol) then
                stop_reason = 1
                exit
            end if

            ! v_n = -G g / sqrt(g^T G g): unit metric speed since |G g|_G = sqrt(g^T G g)
            w(1) = -(g11 * g(1) + g12 * g(2))
            w(2) = -(g12 * g(1) + g22 * g(2))
            v_n(1) = w(1) / gnrm
            v_n(2) = w(2) / gnrm

            ! --- Gonzalez-Schlegel constrained step -----------------------------
            r_p(1) = r(1) + ds_step * v_n(1)
            r_p(2) = r(2) + ds_step * v_n(2)
            ! metric-orthogonal unit tangent t: t^T G_vel v_n = 0, |t|_G = 1
            u_vec(1) = h11 * v_n(1) + h12 * v_n(2)
            u_vec(2) = h12 * v_n(1) + h22 * v_n(2)
            t_vec(1) = -u_vec(2)
            t_vec(2) = u_vec(1)
            tmp = sqrt(max(t_vec(1)**2 * h11 + 2.0_dp * t_vec(1) * t_vec(2) * h12 &
                           + t_vec(2)**2 * h22, 1.0e-300_dp))
            t_vec = t_vec / tmp

            ! psi(xi) = sin(angle( G g(r_p + xi t), chord )) on the bisector,
            ! chord = ds*v_n + xi*t; root => gradient parallel to the chord.
            ! Bracket the root on [-ds, +ds] (fallback: unit-speed Euler step,
            ! which is stable away from the flat floor).
            ibrack = 0
            xl = 0.0_dp; xr = 0.0_dp
            fl = 0.0_dp; fr = 0.0_dp
            do ibis = -1, 1, 2
                xi = real(ibis, dp) * ds_step
                r_new(1) = r_p(1) + xi * t_vec(1)
                r_new(2) = r_p(2) + xi * t_vec(2)
                if (r_new(1) <= r_min .or. r_new(1) >= r_max .or. &
                    r_new(2) <= r_min .or. r_new(2) >= r_max) cycle
                chord(1) = r_new(1) - r(1)
                chord(2) = r_new(2) - r(2)
                dist%r12 = r_new(1); dist%r23 = r_new(2); dist%r31 = r_new(1) + r_new(2)
                g3 = calc_leps_gradient(dist, leps)
                gg(1) = g3(1) + g3(3)
                gg(2) = g3(2) + g3(3)
                w(1) = g11 * gg(1) + g12 * gg(2)
                w(2) = g12 * gg(1) + g22 * gg(2)
                fm = (w(1) * chord(2) - w(2) * chord(1)) &
                     / max(sqrt(w(1)**2 + w(2)**2) * sqrt(chord(1)**2 + chord(2)**2), 1.0e-300_dp)
                if (abs(fm) <= 1.0e-14_dp) then
                    xl = xi
                    fl = 0.0_dp
                    ibrack = 1
                    exit
                end if
                if (ibis == -1) then
                    xl = xi; fl = fm
                else
                    xr = xi; fr = fm
                    if (fl * fr <= 0.0_dp) ibrack = 1
                end if
            end do

            if (ibrack == 1 .and. fl /= 0.0_dp) then
                ! bisection + secant on [xl, xr]
                do ibis = 1, 40
                    xi = 0.5_dp * (xl + xr)
                    r_new(1) = r_p(1) + xi * t_vec(1)
                    r_new(2) = r_p(2) + xi * t_vec(2)
                    chord(1) = r_new(1) - r(1)
                    chord(2) = r_new(2) - r(2)
                    dist%r12 = r_new(1); dist%r23 = r_new(2); dist%r31 = r_new(1) + r_new(2)
                    g3 = calc_leps_gradient(dist, leps)
                    gg(1) = g3(1) + g3(3)
                    gg(2) = g3(2) + g3(3)
                    w(1) = g11 * gg(1) + g12 * gg(2)
                    w(2) = g12 * gg(1) + g22 * gg(2)
                    fm = (w(1) * chord(2) - w(2) * chord(1)) &
                         / max(sqrt(w(1)**2 + w(2)**2) * sqrt(chord(1)**2 + chord(2)**2), 1.0e-300_dp)
                    if (abs(fm) < 1.0e-14_dp .or. 0.5_dp * (xr - xl) < 1.0e-13_dp) exit
                    if (fl * fm <= 0.0_dp) then
                        xr = xi; fr = fm
                    else
                        xl = xi; fl = fm
                    end if
                    ! secant acceleration inside the bracket
                    if (abs(fr - fl) > 1.0e-300_dp) then
                        xi = xl - fl * (xr - xl) / (fr - fl)
                        if (xi > min(xl, xr) .and. xi < max(xl, xr)) then
                            r_new(1) = r_p(1) + xi * t_vec(1)
                            r_new(2) = r_p(2) + xi * t_vec(2)
                            chord(1) = r_new(1) - r(1)
                            chord(2) = r_new(2) - r(2)
                            dist%r12 = r_new(1); dist%r23 = r_new(2); dist%r31 = r_new(1) + r_new(2)
                            g3 = calc_leps_gradient(dist, leps)
                            gg(1) = g3(1) + g3(3)
                            gg(2) = g3(2) + g3(3)
                            w(1) = g11 * gg(1) + g12 * gg(2)
                            w(2) = g12 * gg(1) + g22 * gg(2)
                            fm = (w(1) * chord(2) - w(2) * chord(1)) &
                                 / max(sqrt(w(1)**2 + w(2)**2) * sqrt(chord(1)**2 + chord(2)**2), 1.0e-300_dp)
                            if (fl * fm <= 0.0_dp) then
                                xr = xi; fr = fm
                            else
                                xl = xi; fl = fm
                            end if
                        end if
                    end if
                end do
                xi = 0.5_dp * (xl + xr)
                r_new(1) = r_p(1) + xi * t_vec(1)
                r_new(2) = r_p(2) + xi * t_vec(2)
            else
                ! fallback: unit-speed Euler step (stable far from the flat floor)
                r_eul(1) = r(1) + ds_step * v_n(1)
                r_eul(2) = r(2) + ds_step * v_n(2)
                r_new = r_eul
            end if

            if (r_new(1) <= r_min .or. r_new(1) >= r_max .or. &
                r_new(2) <= r_min .or. r_new(2) >= r_max) then
                stop_reason = 2
                exit
            end if

            dist%r12 = r_new(1); dist%r23 = r_new(2); dist%r31 = r_new(1) + r_new(2)
            v_new = calc_leps_potential(dist, leps)
            if (v_new > v_old + 1.0e-9_dp) then
                stop_reason = 3     ! ascent: ds too large or non-descent direction
                exit
            end if

            s_acc = s_acc + ds_step
            n_out = n_out + 1
            s_tmp(n_out) = s_acc
            r1_tmp(n_out) = r_new(1)
            r2_tmp(n_out) = r_new(2)
            v_tmp(n_out) = v_new
            v_old = v_new
            r = r_new
        end do

        s_path = s_tmp(1:n_out)
        r12_path = r1_tmp(1:n_out)
        r23_path = r2_tmp(1:n_out)
        v_path = v_tmp(1:n_out)
        deallocate(s_tmp, r1_tmp, r2_tmp, v_tmp)
    end subroutine intrinsic_reaction_coordinate

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
