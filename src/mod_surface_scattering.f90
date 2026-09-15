! ==============================================================================
! GeneralModule: mod_surface_scattering.f90
!
! Quantum Surface Scattering: Corrugated Surface Diffraction and
! Selective Adsorption Resonances (SAR)
!
! Theoretical Foundations:
!   1. 2D periodic corrugated surface diffraction (He/LiF, H/metal interfaces).
!      In-plane momentum conservation: K_G = K_i + G, G = (m*b1, n*b2).
!      Normal wavevector k_z,G^2 = 2*M*E - |K_G|^2 (open vs evanescent channels).
!   2. Hard Corrugated Surface (HCS) Eikonal approximation with Bessel decomposition:
!      S_(m,n) = (-i)^(|m|+|n|) * J_m(q_z * zeta_x / 2) * J_n(q_z * zeta_y / 2).
!   3. Selective Adsorption Resonances (SAR / Bound-state resonances):
!      Coupling of evanescent diffraction beams to bound states epsilon_v of the
!      attractive surface potential V_0(z) creates sharp Fano-type resonance dips.
!   4. Debye-Waller thermal attenuation factor exp(-2*W(T)) due to surface phonons.
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_surface_scattering
    use mod_constants, only: dp, PI, TWOPI, KB, AMU2AU, ANG2AU, AU2EV, EV2AU
    implicit none
    private

    public :: surface_lattice_t
    public :: surface_potential_t
    public :: diffraction_beam_t
    public :: init_surface_lattice
    public :: init_surface_potential_morse
    public :: calc_surface_diffraction_channels
    public :: calc_hcs_diffraction_probabilities
    public :: calc_selective_adsorption_resonance
    public :: calc_surface_debye_waller

    !> 2D periodic surface lattice configuration
    type :: surface_lattice_t
        real(dp) :: ax_bohr            !< Surface lattice constant along x [a.u.]
        real(dp) :: ay_bohr            !< Surface lattice constant along y [a.u.]
        real(dp) :: bx_bohr            !< Reciprocal lattice vector length b1 = 2*pi / ax [a.u.]
        real(dp) :: by_bohr            !< Reciprocal lattice vector length b2 = 2*pi / ay [a.u.]
        real(dp) :: zeta_x_bohr        !< Corrugation amplitude along x [a.u.]
        real(dp) :: zeta_y_bohr        !< Corrugation amplitude along y [a.u.]
        real(dp) :: m_substrate_amu    !< Substrate atomic mass [amu]
        real(dp) :: debye_temp_k       !< Surface Debye temperature [K]
    end type surface_lattice_t

    !> Lateral average attractive surface potential (Morse model)
    type :: surface_potential_t
        real(dp) :: well_depth_au      !< Surface well depth D [a.u.]
        real(dp) :: range_param_au     !< Exponential range parameter alpha [a.u.]
        real(dp) :: z_eq_bohr          !< Equilibrium distance z_e [a.u.]
        integer  :: n_bound_states     !< Number of bound states supported
    end type surface_potential_t

    !> Properties of a 2D diffraction channel G = (m, n)
    type :: diffraction_beam_t
        integer  :: m                  !< Reciprocal lattice order along x
        integer  :: n                  !< Reciprocal lattice order along y
        real(dp) :: k_gx_au            !< In-plane wavevector component x [a.u.]
        real(dp) :: k_gy_au            !< In-plane wavevector component y [a.u.]
        real(dp) :: kz_au              !< Normal wavevector component kz [a.u.]
        logical  :: is_open            !< True if channel is propagating (open), False if evanescent
        real(dp) :: probability        !< Normalized diffraction intensity P_(m,n)
        real(dp) :: theta_f_deg        !< Outgoing polar scattering angle [deg]
        real(dp) :: phi_f_deg          !< Outgoing azimuthal scattering angle [deg]
    end type diffraction_beam_t

contains

    !> Initialize 2D rectangular surface lattice (e.g. LiF(001), a = 2.84 Angstrom)
    pure subroutine init_surface_lattice(lat, ax_ang, ay_ang, zeta_x_ang, zeta_y_ang, &
                                         m_sub_amu, debye_temp_k)
        type(surface_lattice_t), intent(out) :: lat
        real(dp), intent(in) :: ax_ang, ay_ang
        real(dp), intent(in) :: zeta_x_ang, zeta_y_ang
        real(dp), intent(in) :: m_sub_amu
        real(dp), intent(in) :: debye_temp_k

        lat%ax_bohr = max(0.5_dp, ax_ang * ANG2AU)
        lat%ay_bohr = max(0.5_dp, ay_ang * ANG2AU)
        lat%bx_bohr = TWOPI / lat%ax_bohr
        lat%by_bohr = TWOPI / lat%ay_bohr
        lat%zeta_x_bohr = zeta_x_ang * ANG2AU
        lat%zeta_y_bohr = zeta_y_ang * ANG2AU
        lat%m_substrate_amu = max(1.0_dp, m_sub_amu)
        lat%debye_temp_k = max(1.0_dp, debye_temp_k)
    end subroutine init_surface_lattice

    !> Initialize lateral average surface Morse potential V_0(z) = D*(exp(-2*alpha*z) - 2*exp(-alpha*z))
    pure subroutine init_surface_potential_morse(pot, well_depth_mev, range_inv_ang, &
                                                mass_amu, n_bound)
        type(surface_potential_t), intent(out) :: pot
        real(dp), intent(in)  :: well_depth_mev  !< Well depth in meV (e.g. 7.5 meV for He/LiF)
        real(dp), intent(in)  :: range_inv_ang   !< Range parameter alpha in 1/Angstrom (e.g. 1.1 A^-1)
        real(dp), intent(in)  :: mass_amu        !< Incident projectile mass [amu]
        integer,  intent(out) :: n_bound         !< Number of bound states

        real(dp) :: d_au, alpha_au, mass_au, param_xi

        d_au = (well_depth_mev * 1.0e-3_dp) * EV2AU
        alpha_au = range_inv_ang / ANG2AU
        mass_au = mass_amu * AMU2AU

        pot%well_depth_au = max(1.0e-8_dp, d_au)
        pot%range_param_au = max(1.0e-4_dp, alpha_au)
        pot%z_eq_bohr = 0.0_dp

        ! Morse bound state parameter: xi = sqrt(2*M*D) / alpha
        param_xi = sqrt(2.0_dp * mass_au * pot%well_depth_au) / pot%range_param_au
        n_bound = max(1, int(param_xi - 0.5_dp) + 1)
        pot%n_bound_states = n_bound
    end subroutine init_surface_potential_morse

    !> Evaluate Morse potential bound state energies:
    !> epsilon_v = - D * [1 - alpha / sqrt(2*M*D) * (v + 1/2)]^2
    pure function calc_morse_bound_energy(pot, mass_amu, v_quant) result(e_bound_au)
        type(surface_potential_t), intent(in) :: pot
        real(dp), intent(in) :: mass_amu
        integer,  intent(in) :: v_quant
        real(dp) :: e_bound_au

        real(dp) :: mass_au, factor

        mass_au = mass_amu * AMU2AU
        factor = 1.0_dp - (pot%range_param_au / sqrt(2.0_dp * mass_au * pot%well_depth_au)) * &
                 (real(v_quant, dp) + 0.5_dp)
        if (factor > 0.0_dp) then
            e_bound_au = - pot%well_depth_au * (factor**2)
        else
            e_bound_au = 0.0_dp
        end if
    end function calc_morse_bound_energy

    !> Enumerate 2D diffraction channels G = (m, n) for given incident beam
    pure subroutine calc_surface_diffraction_channels(lat, mass_amu, energy_ev, &
                                                     theta_i_deg, phi_i_deg, max_order, &
                                                     n_channels, channels)
        type(surface_lattice_t), intent(in) :: lat
        real(dp), intent(in)  :: mass_amu
        real(dp), intent(in)  :: energy_ev
        real(dp), intent(in)  :: theta_i_deg, phi_i_deg
        integer,  intent(in)  :: max_order
        integer,  intent(out) :: n_channels
        type(diffraction_beam_t), intent(out) :: channels((2*max_order + 1)**2)

        real(dp) :: energy_au, mass_au, k_tot, th_i, ph_i
        real(dp) :: k_ix, k_iy, k_iz
        real(dp) :: gx, gy, k_gx, k_gy, kz2, k_norm, k_par
        integer  :: m, n, idx

        energy_au = energy_ev * EV2AU
        mass_au = mass_amu * AMU2AU
        k_tot = sqrt(2.0_dp * mass_au * max(1.0e-12_dp, energy_au))

        th_i = theta_i_deg * (PI / 180.0_dp)
        ph_i = phi_i_deg * (PI / 180.0_dp)

        ! In-plane incident wavevector
        k_ix = k_tot * sin(th_i) * cos(ph_i)
        k_iy = k_tot * sin(th_i) * sin(ph_i)
        k_iz = k_tot * cos(th_i)

        idx = 0
        do m = -max_order, max_order
            do n = -max_order, max_order
                idx = idx + 1
                gx = real(m, dp) * lat%bx_bohr
                gy = real(n, dp) * lat%by_bohr
                k_gx = k_ix + gx
                k_gy = k_iy + gy

                ! Normal momentum conservation: kz^2 = k_tot^2 - (k_gx^2 + k_gy^2)
                kz2 = (k_tot**2) - (k_gx**2 + k_gy**2)

                channels(idx)%m = m
                channels(idx)%n = n
                channels(idx)%k_gx_au = k_gx
                channels(idx)%k_gy_au = k_gy
                channels(idx)%probability = 0.0_dp

                if (kz2 >= 0.0_dp) then
                    channels(idx)%is_open = .true.
                    channels(idx)%kz_au = sqrt(kz2)
                    k_par = sqrt(k_gx**2 + k_gy**2)
                    channels(idx)%theta_f_deg = atan2(k_par, channels(idx)%kz_au) * (180.0_dp / PI)
                    channels(idx)%phi_f_deg   = atan2(k_gy, k_gx) * (180.0_dp / PI)
                else
                    channels(idx)%is_open = .false.
                    channels(idx)%kz_au = -sqrt(-kz2)  ! Imaginary decay
                    channels(idx)%theta_f_deg = 90.0_dp
                    channels(idx)%phi_f_deg   = atan2(k_gy, k_gx) * (180.0_dp / PI)
                end if
            end do
        end do
        n_channels = idx
    end subroutine calc_surface_diffraction_channels

    !> Compute Hard Corrugated Surface (HCS) Eikonal diffraction probabilities
    !> Uses Bessel decomposition: S_(m,n) ~ J_m(qz * zeta_x / 2) * J_n(qz * zeta_y / 2)
    pure subroutine calc_hcs_diffraction_probabilities(lat, mass_amu, energy_ev, &
                                                      theta_i_deg, phi_i_deg, max_order, &
                                                      n_channels, channels)
        type(surface_lattice_t), intent(in) :: lat
        real(dp), intent(in)  :: mass_amu
        real(dp), intent(in)  :: energy_ev
        real(dp), intent(in)  :: theta_i_deg, phi_i_deg
        integer,  intent(in)  :: max_order
        integer,  intent(out) :: n_channels
        type(diffraction_beam_t), intent(inout) :: channels((2*max_order + 1)**2)

        real(dp) :: energy_au, mass_au, k_tot, k_iz, q_z
        real(dp) :: arg_x, arg_y, j_m, j_n, unnorm_p, sum_prob
        integer  :: i

        energy_au = energy_ev * EV2AU
        mass_au = mass_amu * AMU2AU
        k_tot = sqrt(2.0_dp * mass_au * max(1.0e-12_dp, energy_au))
        k_iz = k_tot * cos(theta_i_deg * (PI / 180.0_dp))

        call calc_surface_diffraction_channels(lat, mass_amu, energy_ev, theta_i_deg, phi_i_deg, &
                                               max_order, n_channels, channels)

        sum_prob = 0.0_dp
        do i = 1, n_channels
            if (channels(i)%is_open) then
                ! Momentum transfer along z: qz = k_iz + kz,G
                q_z = k_iz + channels(i)%kz_au
                arg_x = 0.5_dp * q_z * lat%zeta_x_bohr
                arg_y = 0.5_dp * q_z * lat%zeta_y_bohr

                j_m = bessel_jn_approx(channels(i)%m, arg_x)
                j_n = bessel_jn_approx(channels(i)%n, arg_y)

                ! P_G = (kz,G / kiz) * |S_G|^2
                unnorm_p = (channels(i)%kz_au / max(1.0e-10_dp, k_iz)) * (j_m * j_n)**2
                channels(i)%probability = unnorm_p
                sum_prob = sum_prob + unnorm_p
            else
                channels(i)%probability = 0.0_dp
            end if
        end do

        ! Unitary renormalization over all open channels
        if (sum_prob > 1.0e-14_dp) then
            do i = 1, n_channels
                if (channels(i)%is_open) then
                    channels(i)%probability = channels(i)%probability / sum_prob
                end if
            end do
        end if
    end subroutine calc_hcs_diffraction_probabilities

    !> Simple robust pure Bessel function of integer order J_n(x) for small-to-moderate arguments
    pure function bessel_jn_approx(n_order, x_val) result(jn)
        integer,  intent(in) :: n_order
        real(dp), intent(in) :: x_val
        real(dp) :: jn

        integer :: n_abs, k
        real(dp) :: term, s, x_half

        n_abs = abs(n_order)
        x_half = 0.5_dp * x_val

        ! Taylor series expansion: J_n(x) = sum_{k=0}^M (-1)^k / (k! (n+k)!) * (x/2)^(n+2k)
        term = 1.0_dp
        do k = 1, n_abs
            term = term * x_half / real(k, dp)
        end do

        s = term
        do k = 1, 20
            term = - term * (x_half**2) / (real(k, dp) * real(n_abs + k, dp))
            s = s + term
            if (abs(term) < 1.0e-16_dp * abs(s)) exit
        end do

        if (n_order < 0 .and. mod(n_abs, 2) == 1) then
            jn = - s
        else
            jn = s
        end if
    end function bessel_jn_approx

    !> Evaluate Selective Adsorption Resonance (SAR):
    !> Checks if evanescent diffraction channel G matches Morse surface bound state epsilon_v,
    !> and computes Fano lineshape modulation of specular intensity.
    pure subroutine calc_selective_adsorption_resonance(lat, pot, mass_amu, energy_ev, &
                                                       theta_i_deg, phi_i_deg, m_res, n_res, &
                                                       v_bound, is_near_res, delta_e_mev, &
                                                       fano_specular_ratio)
        type(surface_lattice_t),   intent(in)  :: lat
        type(surface_potential_t), intent(in)  :: pot
        real(dp), intent(in)  :: mass_amu
        real(dp), intent(in)  :: energy_ev
        real(dp), intent(in)  :: theta_i_deg, phi_i_deg
        integer,  intent(in)  :: m_res, n_res     !< Evanescent channel (m, n)
        integer,  intent(in)  :: v_bound          !< Surface bound state index (v = 0, 1, ...)
        logical,  intent(out) :: is_near_res      !< True if within resonance window
        real(dp), intent(out) :: delta_e_mev      !< Energy mismatch Delta E = E_z,G - epsilon_v [meV]
        real(dp), intent(out) :: fano_specular_ratio !< Modulated I_specular / I_0

        real(dp) :: energy_au, mass_au, k_tot, th_i, ph_i
        real(dp) :: k_ix, k_iy, gx, gy, k_gx, k_gy, e_z_g, eps_v
        real(dp) :: gamma_res, q_fano, epsilon_red

        energy_au = energy_ev * EV2AU
        mass_au = mass_amu * AMU2AU
        k_tot = sqrt(2.0_dp * mass_au * max(1.0e-12_dp, energy_au))

        th_i = theta_i_deg * (PI / 180.0_dp)
        ph_i = phi_i_deg * (PI / 180.0_dp)
        k_ix = k_tot * sin(th_i) * cos(ph_i)
        k_iy = k_tot * sin(th_i) * sin(ph_i)

        gx = real(m_res, dp) * lat%bx_bohr
        gy = real(n_res, dp) * lat%by_bohr
        k_gx = k_ix + gx
        k_gy = k_iy + gy

        ! Normal energy of diffraction beam: E_z,G = E - |K_G|^2 / (2*M)
        e_z_g = energy_au - (k_gx**2 + k_gy**2) / (2.0_dp * mass_au)
        eps_v = calc_morse_bound_energy(pot, mass_amu, v_bound)

        delta_e_mev = (e_z_g - eps_v) * AU2EV * 1.0e3_dp

        ! Resonance line width Gamma ~ 0.2 meV, Fano q ~ -1.5 (typical for He/LiF)
        gamma_res = 0.25_dp
        q_fano    = - 1.20_dp

        epsilon_red = delta_e_mev / (0.5_dp * gamma_res)
        fano_specular_ratio = (q_fano + epsilon_red)**2 / (1.0_dp + epsilon_red**2)

        ! Normalized to baseline 1.0 off resonance
        fano_specular_ratio = fano_specular_ratio / (q_fano**2 + 1.0_dp)
        is_near_res = (abs(delta_e_mev) < 3.0_dp * gamma_res)
    end subroutine calc_selective_adsorption_resonance

    !> Calculate Debye-Waller attenuation factor exp(-2*W)
    !> 2*W(T) = 3 * hbar^2 * (Delta k_z)^2 * T / (M_sub * k_B * Theta_D^2)
    pure function calc_surface_debye_waller(lat, mass_amu, k_iz_au, kz_g_au, temp_k) result(dw_factor)
        type(surface_lattice_t), intent(in) :: lat
        real(dp), intent(in) :: mass_amu
        real(dp), intent(in) :: k_iz_au, kz_g_au
        real(dp), intent(in) :: temp_k
        real(dp) :: dw_factor

        real(dp) :: delta_kz_m, m_sub_kg, theta_d, t_surf
        real(dp) :: two_w

        ! Delta kz in SI (m^-1)
        delta_kz_m = (k_iz_au + abs(kz_g_au)) * (1.0_dp / 0.529177210903e-10_dp)
        m_sub_kg = lat%m_substrate_amu * 1.66053906660e-27_dp
        theta_d = lat%debye_temp_k
        t_surf = max(1.0_dp, temp_k)

        ! 2*W = 3 * hbar^2 * (Delta kz)^2 * T / (M * kB * Theta_D^2)
        two_w = (3.0_dp * (1.054571817e-34_dp**2) * (delta_kz_m**2) * t_surf) / &
                (m_sub_kg * KB * (theta_d**2))

        dw_factor = exp(-min(40.0_dp, two_w))
    end function calc_surface_debye_waller

end module mod_surface_scattering
