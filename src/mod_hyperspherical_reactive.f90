! ==============================================================================
! GeneralModule: mod_hyperspherical_reactive.f90
!
! Reactive Triatomic Quantum Scattering Dynamics in Hyperspherical Coordinates
!
! Theoretical Foundations:
!   - B. R. Johnson, J. Chem. Phys. 73, 5051 (1980) [Hyperspherical reactive scattering]
!   - R. T. Pack & G. A. Parker, J. Chem. Phys. 87, 3888 (1987) [Quantum reactive scattering]
!   - W. H. Miller, J. Chem. Phys. 62, 1899 (1975) [Semiclassical transition state theory]
!   - D. G. Truhlar & B. C. Garrett, Acc. Chem. Res. 13, 440 (1980) [Variational TST]
!
! Key Capabilities:
!   1. Delves mass-scaled hyperspherical coordinates (rho, alpha) and Jacobi transformation
!   2. Reaction skew angle beta_skew for arbitrary A + BC mass combinations (H+H2 -> 60 deg)
!   3. Eckart and parabolic quantum tunneling transmission probability P(E) through barrier
!   4. Cumulative reaction probability N(E) summing over quantized transition state channels
!   5. Canonical thermal reaction rate constant k(T) via exact energy integration and TST
!
! Standard: Fortran 2008
! ==============================================================================

module mod_hyperspherical_reactive
    use mod_constants, only: dp, PI, TWOPI, HBAR, AMU2AU, K2AU
    implicit none
    private

    ! Public derived types
    public :: reaction_mass_t, transition_state_t

    ! Public procedures
    public :: init_reaction_mass
    public :: jacobi_to_hyperspherical
    public :: hyperspherical_to_jacobi
    public :: calc_eckart_transmission
    public :: calc_cumulative_reaction_probability
    public :: calc_canonical_rate_constant
    public :: calc_tst_wigner_rate

    ! Reaction kinematic mass parameters
    type :: reaction_mass_t
        real(dp) :: mass_a_au      ! Mass of incident atom A
        real(dp) :: mass_b_au      ! Mass of central atom B
        real(dp) :: mass_c_au      ! Mass of target atom C
        real(dp) :: total_mass_au  ! Total mass M = mA + mB + mC
        real(dp) :: mu_bc_au       ! Diatom reduced mass mu_BC = mB*mC / (mB + mC)
        real(dp) :: mu_a_bc_au     ! Atom-diatom reduced mass mu_{A,BC}
        real(dp) :: mu_three_au    ! Three-body reduced mass mu = sqrt(mA*mB*mC / M)
        real(dp) :: scale_factor_d ! Delves mass-scaling factor d = (mu_{A,BC} / mu_BC)^(1/4)
        real(dp) :: skew_angle_rad ! Reaction skew angle beta = arctan(sqrt(mB*M / (mA*mC)))
        real(dp) :: skew_angle_deg ! Skew angle in degrees
    end type reaction_mass_t

    ! Transition state parameters
    type :: transition_state_t
        real(dp) :: v_barrier_au    ! Classical saddle point barrier height V^# in a.u.
        real(dp) :: omega_im_au     ! Imaginary barrier frequency omega^# in a.u. (reaction coordinate)
        real(dp) :: omega_bend_au   ! Bending vibrational frequency at TS in a.u.
        real(dp) :: omega_symm_au   ! Symmetric stretch frequency at TS in a.u.
        integer  :: n_trans_states  ! Number of quantized transverse levels to sum
    end type transition_state_t

    ! In atomic units, hbar = 1.0, k_B = 1.0 (when temperature is in a.u.)
    real(dp), parameter :: HBAR_AU = 1.0_dp
    ! Conversion from atomic unit of bimolecular rate (a_0^3 / (a.u. time)) to cm^3 / (molecule * s)
    ! 1 a.u. volume = (5.291772109e-9 cm)^3 = 1.4818471e-25 cm^3
    ! 1 a.u. time = 2.4188843265857e-17 s
    ! 1 a.u. rate = 1.4818471e-25 / 2.4188843e-17 ~ 6.1260468e-9 cm^3 / s
    real(dp), parameter :: AU2CM3_S = 6.126046835e-9_dp

contains

    ! ==========================================================================
    ! init_reaction_mass:
    ! Initializes mass-scaling parameters and reaction skew angle
    ! ==========================================================================
    pure subroutine init_reaction_mass(ma_au, mb_au, mc_au, rmass)
        real(dp), intent(in)            :: ma_au, mb_au, mc_au
        type(reaction_mass_t), intent(out) :: rmass

        real(dp) :: m_tot, mu_bc, mu_a_bc, sin_beta, cos_beta

        rmass%mass_a_au = ma_au
        rmass%mass_b_au = mb_au
        rmass%mass_c_au = mc_au

        m_tot = ma_au + mb_au + mc_au
        rmass%total_mass_au = m_tot

        mu_bc = (mb_au * mc_au) / (mb_au + mc_au)
        rmass%mu_bc_au = mu_bc

        mu_a_bc = (ma_au * (mb_au + mc_au)) / m_tot
        rmass%mu_a_bc_au = mu_a_bc

        rmass%mu_three_au = sqrt((ma_au * mb_au * mc_au) / m_tot)
        rmass%scale_factor_d = (mu_a_bc / mu_bc)**0.25_dp

        ! Reaction skew angle between reactants (A + BC) and products (AB + C):
        !   sin(beta) = sqrt( mB * M / ( (mA + mB) * (mB + mc) ) )
        !   cos(beta) = sqrt( mA * mC / ( (mA + mB) * (mB + mc) ) )
        sin_beta = sqrt((mb_au * m_tot) / ((ma_au + mb_au) * (mb_au + mc_au)))
        cos_beta = sqrt((ma_au * mc_au) / ((ma_au + mb_au) * (mb_au + mc_au)))

        rmass%skew_angle_rad = atan2(sin_beta, cos_beta)
        rmass%skew_angle_deg = rmass%skew_angle_rad * (180.0_dp / PI)
    end subroutine init_reaction_mass

    ! ==========================================================================
    ! jacobi_to_hyperspherical:
    ! Converts atom-diatom Jacobi coordinates (r, R) to mass-scaled Delves (rho, alpha)
    !   S = d * R,   s = r / d
    !   rho = sqrt(S^2 + s^2),   alpha = arctan(s / S)
    ! ==========================================================================
    pure subroutine jacobi_to_hyperspherical(r_diatom, r_atom_diatom, d_scale, rho, alpha)
        real(dp), intent(in)  :: r_diatom        ! Diatom distance r (a.u.)
        real(dp), intent(in)  :: r_atom_diatom   ! Atom-diatom distance R (a.u.)
        real(dp), intent(in)  :: d_scale         ! Kinematic scaling factor d
        real(dp), intent(out) :: rho             ! Hyperradius rho (a.u.)
        real(dp), intent(out) :: alpha           ! Delves hyperangle alpha (rad)

        real(dp) :: cap_s, small_s

        cap_s   = d_scale * r_atom_diatom
        small_s = r_diatom / d_scale

        rho   = sqrt(cap_s**2 + small_s**2)
        alpha = atan2(small_s, cap_s)
    end subroutine jacobi_to_hyperspherical

    ! ==========================================================================
    ! hyperspherical_to_jacobi:
    ! Converts mass-scaled Delves (rho, alpha) back to Jacobi coordinates (r, R)
    ! ==========================================================================
    pure subroutine hyperspherical_to_jacobi(rho, alpha, d_scale, r_diatom, r_atom_diatom)
        real(dp), intent(in)  :: rho, alpha, d_scale
        real(dp), intent(out) :: r_diatom, r_atom_diatom

        real(dp) :: cap_s, small_s

        cap_s   = rho * cos(alpha)
        small_s = rho * sin(alpha)

        r_atom_diatom = cap_s / d_scale
        r_diatom      = small_s * d_scale
    end subroutine hyperspherical_to_jacobi

    ! ==========================================================================
    ! calc_eckart_transmission:
    ! Calculates quantum transmission probability P(E) through an inverted parabolic
    ! or Eckart barrier with imaginary frequency omega^#:
    !   P(E) = 1 / [ 1 + exp( - 2*pi*(E - V_b) / (hbar * omega^#) ) ]
    ! ==========================================================================
    pure function calc_eckart_transmission(energy_au, v_barrier_au, omega_im_au) result(prob)
        real(dp), intent(in) :: energy_au, v_barrier_au, omega_im_au
        real(dp) :: prob

        real(dp) :: arg

        if (omega_im_au <= 0.0_dp) then
            if (energy_au >= v_barrier_au) then
                prob = 1.0_dp
            else
                prob = 0.0_dp
            end if
            return
        end if

        arg = - TWOPI * (energy_au - v_barrier_au) / (HBAR_AU * omega_im_au)

        ! Avoid numerical overflow in exp(arg)
        if (arg > 60.0_dp) then
            prob = 0.0_dp
        else if (arg < -60.0_dp) then
            prob = 1.0_dp
        else
            prob = 1.0_dp / (1.0_dp + exp(arg))
        end if
    end function calc_eckart_transmission

    ! ==========================================================================
    ! calc_cumulative_reaction_probability:
    ! Calculates cumulative reaction probability N(E) = sum_n P_n(E)
    ! summing over quantized transverse modes (symmetric stretch + doubly degenerate bend)
    ! ==========================================================================
    pure function calc_cumulative_reaction_probability(ts, energy_au) result(n_prob)
        type(transition_state_t), intent(in) :: ts
        real(dp), intent(in)                 :: energy_au
        real(dp)                             :: n_prob

        integer  :: v_sym, v_bend, d_bend
        real(dp) :: e_trans, v_eff_barrier, p_trans

        n_prob = 0.0_dp

        ! Sum over symmetric stretch modes and bending modes
        do v_sym = 0, ts%n_trans_states
            do v_bend = 0, ts%n_trans_states
                ! Degeneracy for 2D isotropic bending mode: g_b = v_bend + 1
                d_bend = v_bend + 1

                ! Zero-point and vibrational excitation at TS
                e_trans = (real(v_sym, dp) + 0.5_dp) * ts%omega_symm_au + &
                          (real(v_bend, dp) + 1.0_dp) * ts%omega_bend_au

                v_eff_barrier = ts%v_barrier_au + e_trans
                p_trans = calc_eckart_transmission(energy_au, v_eff_barrier, ts%omega_im_au)

                n_prob = n_prob + real(d_bend, dp) * p_trans

                ! Break if transmission for this mode is negligibly small
                if (p_trans < 1.0e-8_dp .and. energy_au < v_eff_barrier) exit
            end do
        end do
    end function calc_cumulative_reaction_probability

    ! ==========================================================================
    ! calc_canonical_rate_constant:
    ! Computes canonical thermal rate constant k(T) by numerical quadrature of N(E):
    !   k(T) = [ 1 / (h * Q_reac(T)) ] * \int_0^\infty N(E) * exp(-E / (k_B * T)) dE
    ! Returns k(T) in cm^3 / (molecule * s)
    ! ==========================================================================
    function calc_canonical_rate_constant(ts, rmass, temp_kelvin, n_e_steps) result(k_rate_cm3_s)
        type(transition_state_t), intent(in) :: ts
        type(reaction_mass_t), intent(in)    :: rmass
        real(dp), intent(in)                 :: temp_kelvin
        integer, intent(in)                  :: n_e_steps
        real(dp)                             :: k_rate_cm3_s

        real(dp) :: k_b_t_au, e_max_au, de_au, integral_sum, energy, n_val, weight
        real(dp) :: q_trans_reac, q_rot_reac, q_vib_reac, q_tot_reac
        real(dp) :: rate_au
        integer  :: i

        k_b_t_au = temp_kelvin * K2AU

        ! Energy integration range up to V_b + 12 * k_B * T
        e_max_au = ts%v_barrier_au + 12.0_dp * k_b_t_au
        de_au = e_max_au / real(n_e_steps, dp)

        integral_sum = 0.0_dp
        do i = 1, n_e_steps
            energy = (real(i, dp) - 0.5_dp) * de_au
            n_val  = calc_cumulative_reaction_probability(ts, energy)

            weight = exp(- energy / k_b_t_au) * de_au
            integral_sum = integral_sum + n_val * weight
        end do

        ! Reactant partition functions (per unit volume for relative translation):
        ! Q_trans = (mu_{A,BC} * k_B * T / (2*pi*hbar^2))^(3/2)
        q_trans_reac = ((rmass%mu_a_bc_au * k_b_t_au) / (TWOPI * (HBAR_AU**2)))**1.5_dp
        ! Typical diatom rotational and vibrational partition functions at T
        q_rot_reac   = max(1.0_dp, k_b_t_au / (0.00027_dp)) ! Assume typical B_rot ~ 60 cm^-1
        q_vib_reac   = 1.0_dp / (1.0_dp - exp(- 0.020_dp / k_b_t_au)) ! Typical omega ~ 4400 cm^-1
        q_tot_reac   = q_trans_reac * q_rot_reac * q_vib_reac

        ! k(T) = integral / (2*pi*hbar * Q_tot)
        rate_au = integral_sum / (TWOPI * HBAR_AU * q_tot_reac)
        k_rate_cm3_s = rate_au * AU2CM3_S
    end function calc_canonical_rate_constant

    ! ==========================================================================
    ! calc_tst_wigner_rate:
    ! Evaluates conventional Transition State Theory rate with Wigner tunneling factor:
    !   kappa_W(T) = 1 + (1/24) * ( (hbar * omega^#) / (k_B * T) )^2
    ! ==========================================================================
    pure function calc_tst_wigner_rate(v_barrier_au, omega_im_au, temp_kelvin, a_prefactor) result(k_tst)
        real(dp), intent(in) :: v_barrier_au
        real(dp), intent(in) :: omega_im_au
        real(dp), intent(in) :: temp_kelvin
        real(dp), intent(in) :: a_prefactor    ! Pre-exponential Arrhenius factor (cm^3/s)
        real(dp)             :: k_tst

        real(dp) :: k_b_t_au, u_tunnel, kappa_w

        k_b_t_au = temp_kelvin * K2AU
        u_tunnel = (HBAR_AU * omega_im_au) / k_b_t_au

        ! Wigner tunneling correction
        kappa_w = 1.0_dp + (u_tunnel**2) / 24.0_dp

        ! Arrhenius rate: k = kappa * A * exp(- V_b / (k_B * T))
        k_tst = kappa_w * a_prefactor * exp(- v_barrier_au / k_b_t_au)
    end function calc_tst_wigner_rate

end module mod_hyperspherical_reactive
