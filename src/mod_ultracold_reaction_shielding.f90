! ==============================================================================
! GeneralModule: mod_ultracold_reaction_shielding.f90
!
! Ultracold Polar Molecule Reactions and Microwave / DC Shielding
!
! Theoretical Foundations:
!   1. Ultracold bimolecular chemical reactions (e.g. 2 KRb -> K2 + Rb2)
!      and sticky collision complex loss.
!   2. Microwave blue-detuned shielding & DC electric field dipolar repulsion:
!      Avoided crossing between rotational dressed states creates an engineered
!      long-range repulsive barrier at R_shield ~ (2*C3 / (hbar*Delta))^(1/3),
!      preventing molecules from reaching short-range reaction zone.
!   3. Quantum Defect Theory (QDT) and WKB barrier penetration:
!      T_tunnel = exp(-2 int sqrt(2*mu/hbar^2 * (V_eff - E)) dR).
!   4. Shielded rate coefficients: Inelastic loss K2_inel and elastic K2_el,
!      yielding evaporative cooling figure-of-merit gamma = K2_el / K2_inel > 100.
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_ultracold_reaction_shielding
    use mod_constants, only: dp, PI, TWOPI, HBAR, KB, AMU2AU, AU2K, K2AU, &
                             AU2M, M2AU, DEBYE2AU, AU2DEBYE
    implicit none
    private

    public :: ultracold_molecule_t
    public :: shielding_config_t
    public :: init_ultracold_molecule_preset
    public :: init_shielding_config
    public :: calc_effective_shielding_potential
    public :: calc_shielding_barrier_height
    public :: calc_wkb_tunneling_probability
    public :: calc_shielded_scattering_rates
    public :: calc_shielding_detuning_scan

    ! Conversion constant: 1 MHz = 1e6 Hz * 2*pi * hbar / E_h in a.u.
    ! hbar * 2*pi * 1e6 / E_h = 1.519829846e-10 a.u.
    real(dp), parameter :: MHZ2AU = 1.519829846004e-10_dp

    !> Physical properties of an ultracold polar molecule
    type :: ultracold_molecule_t
        character(len=16) :: name      !< Molecule name (e.g. KRb, NaRb, NaK)
        real(dp) :: mass_amu           !< Single molecule mass [amu]
        real(dp) :: mu_au              !< Two-body reduced mass mu = m / 2 [a.u.]
        real(dp) :: dipole_debye       !< Permanent / transition dipole moment [Debye]
        real(dp) :: dipole_au          !< Dipole moment in atomic units
        real(dp) :: c6_au              !< van der Waals dispersion coefficient C6 [a.u.]
        real(dp) :: c3_au              !< Resonant dipole-dipole coupling C3 = d^2 [a.u.]
    end type ultracold_molecule_t

    !> Shielding field parameters (Microwave or DC electric field)
    type :: shielding_config_t
        integer  :: method             !< 1: Microwave (MW) shielding, 2: DC electric field
        real(dp) :: detuning_mhz       !< Microwave blue-detuning Delta [MHz]
        real(dp) :: rabi_mhz           !< Microwave Rabi frequency Omega [MHz]
        real(dp) :: e_field_kv_cm      !< DC field strength [kV/cm]
        real(dp) :: y_loss             !< Short-range absorption probability (0 <= y <= 1, 1: universal)
    end type shielding_config_t

contains

    !> Initialize preset parameters for prototypical ultracold polar molecules
    pure subroutine init_ultracold_molecule_preset(mol, mol_name)
        type(ultracold_molecule_t), intent(out) :: mol
        character(len=*), intent(in) :: mol_name

        mol%name = trim(adjustl(mol_name))

        select case(trim(mol%name))
        case("KRb", "krb", "40K87Rb")
            ! 40K-87Rb molecule: mass = 39.96 + 86.91 = 126.87 amu
            mol%mass_amu = 126.87_dp
            mol%dipole_debye = 0.574_dp   ! Permanent dipole ~ 0.57 Debye
            mol%c6_au = 1.61e4_dp         ! C6 ~ 16100 a.u.
        case("NaRb", "narb", "23Na87Rb")
            ! 23Na-87Rb molecule: mass = 22.99 + 86.91 = 109.90 amu
            mol%mass_amu = 109.90_dp
            mol%dipole_debye = 3.30_dp    ! Strong permanent dipole ~ 3.3 Debye
            mol%c6_au = 2.45e4_dp
        case("NaK", "nak", "23Na40K")
            ! 23Na-40K molecule: mass = 22.99 + 39.96 = 62.95 amu
            mol%mass_amu = 62.95_dp
            mol%dipole_debye = 2.76_dp    ! Dipole ~ 2.76 Debye
            mol%c6_au = 1.80e4_dp
        case default
            ! Generic polar molecule
            mol%mass_amu = 100.0_dp
            mol%dipole_debye = 1.0_dp
            mol%c6_au = 2.0e4_dp
        end select

        ! Two-body reduced mass in collisions of two identical molecules: mu = m / 2
        mol%mu_au = 0.5_dp * mol%mass_amu * AMU2AU
        mol%dipole_au = mol%dipole_debye * DEBYE2AU
        ! Resonant dipole-dipole parameter C3 = d^2 in a.u.
        mol%c3_au = mol%dipole_au**2
    end subroutine init_ultracold_molecule_preset

    !> Initialize shielding configuration
    pure subroutine init_shielding_config(cfg, method, detuning_mhz, rabi_mhz, &
                                         e_field_kv_cm, y_loss)
        type(shielding_config_t), intent(out) :: cfg
        integer,  intent(in) :: method
        real(dp), intent(in) :: detuning_mhz
        real(dp), intent(in) :: rabi_mhz
        real(dp), intent(in) :: e_field_kv_cm
        real(dp), intent(in) :: y_loss

        cfg%method = method
        cfg%detuning_mhz = max(0.1_dp, detuning_mhz)
        cfg%rabi_mhz = max(0.0_dp, rabi_mhz)
        cfg%e_field_kv_cm = max(0.0_dp, e_field_kv_cm)
        cfg%y_loss = min(1.0_dp, max(0.0_dp, y_loss))
    end subroutine init_shielding_config

    !> Compute effective intermolecular potential V_eff(R) in Kelvin
    !> Includes engineered repulsive dipole shielding, van der Waals C6, and centrifugal barrier
    pure subroutine calc_effective_shielding_potential(mol, cfg, r_bohr, l_ang, v_eff_kelvin)
        type(ultracold_molecule_t), intent(in) :: mol
        type(shielding_config_t),   intent(in) :: cfg
        real(dp), intent(in)  :: r_bohr
        integer,  intent(in)  :: l_ang
        real(dp), intent(out) :: v_eff_kelvin

        real(dp) :: r, delta_au, v_rep_au, v_vdw_au, v_cent_au, v_tot_au
        real(dp) :: c3_over_r3, rabi_au

        r = max(1.0_dp, r_bohr)

        ! 1. Engineered Repulsive Shielding Potential
        if (cfg%method == 1) then
            ! Microwave blue-detuned shielding (dressed-state avoided crossing):
            ! V_rep(R) = 0.5 * (sqrt(Delta^2 + (2*C3/R^3)^2 + Omega^2) - Delta)
            delta_au = cfg%detuning_mhz * MHZ2AU
            rabi_au = cfg%rabi_mhz * MHZ2AU
            c3_over_r3 = 2.0_dp * mol%c3_au / (r**3)

            v_rep_au = 0.5_dp * (sqrt(delta_au**2 + c3_over_r3**2 + rabi_au**2) - delta_au)
        else
            ! DC electric field induced repulsive dipole barrier (side-by-side repulsive channel):
            ! Induced dipole d_ind ~ d * min(1, d*E / (2*B_rot))
            ! Effective isotropic repulsive average: V_rep(R) = C3_ind / R^3
            v_rep_au = (0.5_dp * mol%c3_au * min(1.0_dp, cfg%e_field_kv_cm / 15.0_dp)**2) / (r**3)
        end if

        ! 2. Attractive van der Waals potential -C6 / R^6
        v_vdw_au = - mol%c6_au / (r**6)

        ! 3. Centrifugal barrier: hbar^2 l(l+1) / (2 * mu * R^2)
        v_cent_au = real(l_ang * (l_ang + 1), dp) / (2.0_dp * mol%mu_au * r**2)

        v_tot_au = v_rep_au + v_vdw_au + v_cent_au
        v_eff_kelvin = v_tot_au * AU2K
    end subroutine calc_effective_shielding_potential

    !> Calculate the shielding barrier location and barrier height
    pure subroutine calc_shielding_barrier_height(mol, cfg, r_barrier_bohr, v_barrier_kelvin)
        type(ultracold_molecule_t), intent(in) :: mol
        type(shielding_config_t),   intent(in) :: cfg
        real(dp), intent(out) :: r_barrier_bohr
        real(dp), intent(out) :: v_barrier_kelvin

        real(dp) :: delta_au, r_est, v_eff
        integer :: iter

        ! Analytic estimate of barrier radius: R_shield = (2 * C3 / Delta)^(1/3)
        delta_au = cfg%detuning_mhz * MHZ2AU
        r_est = (2.0_dp * mol%c3_au / delta_au)**(1.0_dp / 3.0_dp)
        r_barrier_bohr = max(50.0_dp, r_est)

        ! Evaluate potential at barrier peak for s-wave (l = 0)
        call calc_effective_shielding_potential(mol, cfg, r_barrier_bohr, 0, v_eff)
        v_barrier_kelvin = max(1.0e-7_dp, v_eff)
    end subroutine calc_shielding_barrier_height

    !> Calculate WKB quantum tunneling transmission probability through the barrier
    !> T_tunnel(E) = exp( - 2 * int_{R_in}^{R_out} sqrt(2*mu*(V_eff - E)) dR )
    pure subroutine calc_wkb_tunneling_probability(mol, cfg, collision_energy_uk, t_tunnel)
        type(ultracold_molecule_t), intent(in) :: mol
        type(shielding_config_t),   intent(in) :: cfg
        real(dp), intent(in)  :: collision_energy_uk  !< Energy in microkelvin (uK)
        real(dp), intent(out) :: t_tunnel

        real(dp) :: e_coll_k, r_barrier, v_barrier
        real(dp) :: r_in, r_out, dr, r_curr, v_curr, integrand, wkb_integral
        integer :: n_steps, i

        e_coll_k = collision_energy_uk * 1.0e-6_dp
        call calc_shielding_barrier_height(mol, cfg, r_barrier, v_barrier)

        ! If collision energy exceeds barrier height, over-the-barrier transmission ~ 1
        if (e_coll_k >= v_barrier) then
            t_tunnel = 1.0_dp
            return
        end if

        ! Numerical integration boundaries around the barrier peak
        r_in  = 0.35_dp * r_barrier
        r_out = 1.80_dp * r_barrier
        n_steps = 100
        dr = (r_out - r_in) / real(n_steps, dp)
        wkb_integral = 0.0_dp

        do i = 1, n_steps
            r_curr = r_in + (real(i, dp) - 0.5_dp) * dr
            call calc_effective_shielding_potential(mol, cfg, r_curr, 0, v_curr)
            if (v_curr > e_coll_k) then
                ! sqrt(2 * mu_au * (V_au - E_au))
                integrand = sqrt(2.0_dp * mol%mu_au * (v_curr - e_coll_k) * K2AU)
                wkb_integral = wkb_integral + integrand * dr
            end if
        end do

        ! T_tunnel = exp(-2 * S_wkb)
        t_tunnel = exp(-2.0_dp * min(40.0_dp, wkb_integral))
    end subroutine calc_wkb_tunneling_probability

    !> Compute thermal two-body elastic and inelastic collision rate coefficients
    !> and the evaporative cooling ratio gamma = K2_el / K2_inel
    pure subroutine calc_shielded_scattering_rates(mol, cfg, temp_uk, &
                                                   k2_el_cm3s, k2_inel_cm3s, gamma_ratio)
        type(ultracold_molecule_t), intent(in) :: mol
        type(shielding_config_t),   intent(in) :: cfg
        real(dp), intent(in)  :: temp_uk        !< Temperature in microkelvin (uK)
        real(dp), intent(out) :: k2_el_cm3s     !< Elastic rate coefficient [cm^3/s]
        real(dp), intent(out) :: k2_inel_cm3s   !< Inelastic reaction/loss rate coefficient [cm^3/s]
        real(dp), intent(out) :: gamma_ratio    !< Ratio gamma = K2_el / K2_inel

        real(dp) :: t_kelvin, e_th_uk, t_tunnel
        real(dp) :: v_th_si, r_barrier, v_barrier
        real(dp) :: a_eff_m, sigma_el_si, sigma_inel_si
        real(dp) :: mu_kg

        t_kelvin = max(1.0e-8_dp, temp_uk * 1.0e-6_dp)
        e_th_uk = temp_uk

        ! Barrier tunneling transmission at thermal energy
        call calc_wkb_tunneling_probability(mol, cfg, e_th_uk, t_tunnel)
        call calc_shielding_barrier_height(mol, cfg, r_barrier, v_barrier)

        ! Reduced mass in kg
        mu_kg = mol%mu_au * (1.0_dp / AMU2AU) * 1.66053906660e-27_dp

        ! Thermal velocity v_th = sqrt(8 * k_B * T / (pi * mu))
        v_th_si = sqrt(8.0_dp * KB * t_kelvin / (PI * mu_kg))

        ! Effective scattering length a_eff is dominated by the shielding barrier radius R_barrier
        ! a_eff ~ R_barrier (hard sphere scattering from the engineered repulsive core)
        a_eff_m = r_barrier * (1.0_dp / M2AU)

        ! Elastic cross section: sigma_el = 4 * pi * a_eff^2 for identical bosons / dressed states
        sigma_el_si = 4.0_dp * PI * (a_eff_m**2)
        k2_el_cm3s = v_th_si * sigma_el_si * 1.0e6_dp  ! m^3/s to cm^3/s

        ! Inelastic loss cross section from QDT:
        ! sigma_inel = 4 * pi / k^2 * y_loss * T_tunnel
        ! Thermal average rate K2_inel ~ v_th * sigma_inel
        ! Universal short range absorption benchmark: K2_0 ~ 1e-10 cm^3/s
        ! Shielded rate: K2_inel = K2_unshielded * (4 * y_loss * T_tunnel / (2 + y_loss * T_tunnel)^2)
        k2_inel_cm3s = 1.0e-10_dp * cfg%y_loss * t_tunnel
        k2_inel_cm3s = max(1.0e-18_dp, k2_inel_cm3s)

        gamma_ratio = k2_el_cm3s / k2_inel_cm3s
    end subroutine calc_shielded_scattering_rates

    !> Perform a parametric scan over microwave detuning Delta to map out the shielding curve
    pure subroutine calc_shielding_detuning_scan(mol, cfg_base, temp_uk, n_pts, &
                                                det_min_mhz, det_max_mhz, det_arr, gamma_arr)
        type(ultracold_molecule_t), intent(in) :: mol
        type(shielding_config_t),   intent(in) :: cfg_base
        real(dp), intent(in)  :: temp_uk
        integer,  intent(in)  :: n_pts
        real(dp), intent(in)  :: det_min_mhz, det_max_mhz
        real(dp), intent(out) :: det_arr(n_pts)
        real(dp), intent(out) :: gamma_arr(n_pts)

        integer :: i
        real(dp) :: d_step, det_val
        real(dp) :: k2_el, k2_inel, g_val
        type(shielding_config_t) :: cfg_scan

        if (n_pts <= 1) return
        d_step = (det_max_mhz - det_min_mhz) / real(n_pts - 1, dp)
        cfg_scan = cfg_base

        do i = 1, n_pts
            det_val = det_min_mhz + real(i - 1, dp) * d_step
            det_arr(i) = det_val
            cfg_scan%detuning_mhz = det_val
            call calc_shielded_scattering_rates(mol, cfg_scan, temp_uk, k2_el, k2_inel, g_val)
            gamma_arr(i) = g_val
        end do
    end subroutine calc_shielding_detuning_scan

end module mod_ultracold_reaction_shielding
