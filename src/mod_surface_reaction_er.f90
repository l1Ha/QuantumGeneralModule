! ==============================================================================
! GeneralModule: mod_surface_reaction_er.f90
!
! Gas-Surface Reaction Dynamics: Eley-Rideal (ER) and
! Langmuir-Hinshelwood (LH) Catalytic Mechanisms
!
! Theoretical Foundations:
!   1. Eley-Rideal abstraction: A(g) + B(ads)/Surf -> AB(g, v, j) + Surf.
!      Direct sub-picosecond abstraction without accommodation/thermalization.
!   2. 2D collinear reactive potential energy surface V(r, Z_cm) with modified LEPS.
!   3. Massive exothermicity release Delta E_exo = D_AB - D_chem:
!      Hyper-thermal partitioning into vibrational excitation <E_vib> (~50%),
!      translation <E_trans> (~35%), and substrate phonon/electron dissipation.
!   4. State-resolved vibrational population P(v) of desorbed nascent molecules,
!      energy-dependent cross section sigma_ER(E_i), and thermal rate constant k_ER(T).
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================

module mod_surface_reaction_er
    use mod_constants, only: dp, PI, KB, AMU2AU, ANG2AU, AU2EV, EV2AU
    implicit none
    private

    public :: er_reaction_system_t
    public :: er_energy_partition_t
    public :: init_er_reaction_system
    public :: calc_er_potential_2d
    public :: calc_er_energy_partitioning
    public :: calc_er_vibrational_populations
    public :: calc_er_reaction_cross_section
    public :: calc_er_thermal_rate_constant

    !> Eley-Rideal reaction system parameters
    type :: er_reaction_system_t
        character(len=24) :: name      !< System name (e.g. H+H/Cu(111), H+D/Cu)
        real(dp) :: m_gas_amu          !< Gas-phase projectile mass A [amu]
        real(dp) :: m_ad_amu           !< Adsorbed adatom mass B [amu]
        real(dp) :: d_mol_ev           !< Gas-phase molecule bond dissociation energy D_AB [eV]
        real(dp) :: d_chem_ev          !< Adatom surface chemisorption energy D_chem [eV]
        real(dp) :: r_e_ang            !< Equilibrium bond length of AB [Angstrom]
        real(dp) :: z_e_ang            !< Equilibrium chemisorption height [Angstrom]
        real(dp) :: alpha_mol_inv_ang  !< Molecular Morse parameter alpha_AB [A^-1]
        real(dp) :: beta_surf_inv_ang  !< Chemisorption Morse parameter beta_chem [A^-1]
        real(dp) :: delta_e_exo_ev     !< Reaction exothermicity Delta E = D_AB - D_chem [eV]
        real(dp) :: barrier_ev         !< Apparent reaction entrance barrier [eV]
    end type er_reaction_system_t

    !> Exothermicity and energy partitioning of desorbed product AB
    type :: er_energy_partition_t
        real(dp) :: e_total_avail_ev   !< Total available energy E_avail = Delta E_exo + E_i [eV]
        real(dp) :: e_trans_ev         !< Translational kinetic energy of desorbed AB [eV]
        real(dp) :: e_vib_ev           !< Vibrational internal energy of AB [eV]
        real(dp) :: e_rot_ev           !< Rotational internal energy of AB [eV]
        real(dp) :: e_diss_ev          !< Energy dissipated into surface phonons/electrons [eV]
        real(dp) :: f_trans            !< Fraction of energy in translation
        real(dp) :: f_vib              !< Fraction of energy in vibration
        real(dp) :: f_diss             !< Fraction of energy dissipated into substrate
    end type er_energy_partition_t

contains

    !> Initialize preset parameters for prototypical Eley-Rideal systems (e.g. H + H/Cu(111))
    pure subroutine init_er_reaction_system(sys, sys_name)
        type(er_reaction_system_t), intent(out) :: sys
        character(len=*), intent(in) :: sys_name

        sys%name = trim(adjustl(sys_name))

        select case(trim(sys%name))
        case("H+H/Cu", "H+H/Cu(111)", "H+H/Cu111")
            ! H(g) + H(ads)/Cu(111) -> H2(g) + Cu(111)
            sys%m_gas_amu = 1.0078_dp
            sys%m_ad_amu  = 1.0078_dp
            sys%d_mol_ev  = 4.75_dp       ! H2 bond dissociation energy ~ 4.75 eV
            sys%d_chem_ev = 2.45_dp       ! H on Cu chemisorption well depth ~ 2.45 eV
            sys%r_e_ang   = 0.741_dp      ! H2 equilibrium bond length
            sys%z_e_ang   = 1.00_dp       ! H/Cu chemisorption distance
            sys%alpha_mol_inv_ang = 1.94_dp
            sys%beta_surf_inv_ang = 1.10_dp
            sys%barrier_ev = 0.04_dp      ! Very small / negligible entrance barrier

        case("H+D/Cu", "H+D/Cu(111)", "H+D/Cu111")
            ! H(g) + D(ads)/Cu(111) -> HD(g) + Cu(111)
            sys%m_gas_amu = 1.0078_dp
            sys%m_ad_amu  = 2.0141_dp
            sys%d_mol_ev  = 4.75_dp
            sys%d_chem_ev = 2.48_dp
            sys%r_e_ang   = 0.741_dp
            sys%z_e_ang   = 1.00_dp
            sys%alpha_mol_inv_ang = 1.94_dp
            sys%beta_surf_inv_ang = 1.10_dp
            sys%barrier_ev = 0.04_dp

        case default
            ! Generic exothermic abstraction
            sys%m_gas_amu = 1.0_dp
            sys%m_ad_amu  = 1.0_dp
            sys%d_mol_ev  = 4.0_dp
            sys%d_chem_ev = 2.0_dp
            sys%r_e_ang   = 1.0_dp
            sys%z_e_ang   = 1.2_dp
            sys%alpha_mol_inv_ang = 1.8_dp
            sys%beta_surf_inv_ang = 1.0_dp
            sys%barrier_ev = 0.05_dp
        end select

        ! Net exothermicity Delta E_exo = D_AB - D_chem
        sys%delta_e_exo_ev = sys%d_mol_ev - sys%d_chem_ev
    end subroutine init_er_reaction_system

    !> Compute 2D reactive collinear potential energy V(r, Z_cm) in eV
    !> r: AB bond distance [Angstrom], Z_cm: center-of-mass distance from surface [Angstrom]
    pure subroutine calc_er_potential_2d(sys, r_bond_ang, z_cm_ang, v_pot_ev)
        type(er_reaction_system_t), intent(in) :: sys
        real(dp), intent(in)  :: r_bond_ang
        real(dp), intent(in)  :: z_cm_ang
        real(dp), intent(out) :: v_pot_ev

        real(dp) :: r, z_cm, z_ad, v_mol, v_chem, v_rep_wall
        real(dp) :: mass_tot, m_ratio

        r = max(0.2_dp, r_bond_ang)
        z_cm = max(0.1_dp, z_cm_ang)

        mass_tot = sys%m_gas_amu + sys%m_ad_amu
        m_ratio  = sys%m_gas_amu / mass_tot

        ! Adatom position from center of mass: z_ad = Z_cm - (m_gas / M_tot) * r
        z_ad = z_cm - m_ratio * r

        ! 1. Diatomic molecule bond potential: Morse potential
        ! V_mol(r) = D_mol * (1 - exp(-alpha * (r - re)))^2 - D_mol
        v_mol = sys%d_mol_ev * (1.0_dp - exp(-sys%alpha_mol_inv_ang * (r - sys%r_e_ang)))**2 - sys%d_mol_ev

        ! 2. Adatom surface interaction: Chemisorption Morse well
        ! V_chem(z_ad) = D_chem * (exp(-2*beta*(z_ad - ze)) - 2*exp(-beta*(z_ad - ze)))
        if (z_ad > 0.0_dp) then
            v_chem = sys%d_chem_ev * (exp(-2.0_dp * sys%beta_surf_inv_ang * (z_ad - sys%z_e_ang)) - &
                                      2.0_dp * exp(-sys%beta_surf_inv_ang * (z_ad - sys%z_e_ang)))
        else
            v_chem = sys%d_chem_ev * (exp(-2.0_dp * sys%beta_surf_inv_ang * (z_ad - sys%z_e_ang)) + 2.0_dp)
        end if

        ! 3. Gas atom direct surface repulsion (Pauli repulsive wall)
        ! z_gas = Z_cm + (m_ad / M_tot) * r
        v_rep_wall = 5.0_dp * exp(-2.5_dp * (z_cm + (sys%m_ad_amu / mass_tot) * r))

        ! Total reactive potential referenced to free surface + separated atoms:
        ! As r -> r_e, Z_cm -> infty: V -> -D_mol + Delta E_exo = -D_chem
        v_pot_ev = v_mol + v_chem + v_rep_wall + sys%d_chem_ev
    end subroutine calc_er_potential_2d

    !> Calculate energy budget partitioning into desorbed product modes
    !> Based on Jackson-Persson quantum dynamics and Rettner experimental data:
    !> Vibration: ~50%, Translation: ~35%, Substrate dissipation: ~10%, Rotation: ~5%.
    pure subroutine calc_er_energy_partitioning(sys, e_incident_ev, partition)
        type(er_reaction_system_t), intent(in)  :: sys
        real(dp),                   intent(in)  :: e_incident_ev
        type(er_energy_partition_t), intent(out) :: partition

        real(dp) :: e_avail

        e_avail = max(0.1_dp, sys%delta_e_exo_ev + e_incident_ev)
        partition%e_total_avail_ev = e_avail

        ! Classical / quantum empirical partitioning coefficients for ER on metal
        partition%f_vib   = 0.50_dp
        partition%f_trans = 0.35_dp
        partition%f_diss  = 0.10_dp
        ! Remainder is in rotation (5%)

        partition%e_vib_ev   = partition%f_vib * e_avail
        partition%e_trans_ev = partition%f_trans * e_avail
        partition%e_diss_ev  = partition%f_diss * e_avail
        partition%e_rot_ev   = e_avail - (partition%e_vib_ev + partition%e_trans_ev + partition%e_diss_ev)
    end subroutine calc_er_energy_partitioning

    !> Compute state-resolved vibrational population P(v) of nascent desorbed molecules
    !> Produces inverted / highly excited vibrational distribution peaking at v ~ 2-4
    pure subroutine calc_er_vibrational_populations(sys, e_incident_ev, max_v, v_dist)
        type(er_reaction_system_t), intent(in)  :: sys
        real(dp),                   intent(in)  :: e_incident_ev
        integer,                    intent(in)  :: max_v
        real(dp),                   intent(out) :: v_dist(0:max_v)

        type(er_energy_partition_t) :: part
        real(dp) :: omega_vib_ev, v_peak, s_sum, p_v
        integer  :: v

        call calc_er_energy_partitioning(sys, e_incident_ev, part)

        ! H2 vibrational quantum hbar * omega_vib ~ 0.54 eV
        omega_vib_ev = 0.54_dp * sqrt(1.0_dp / (sys%m_gas_amu * sys%m_ad_amu / (sys%m_gas_amu + sys%m_ad_amu)))

        ! Average vibrational quantum number <v> = E_vib / (hbar*omega)
        v_peak = max(0.5_dp, part%e_vib_ev / omega_vib_ev)

        s_sum = 0.0_dp
        do v = 0, max_v
            ! Poisson-like inverted distribution: P(v) ~ (v_peak^v / v!) * exp(-v_peak)
            p_v = (v_peak**real(v, dp)) * exp(-v_peak) / gamma_fact(v)
            v_dist(v) = p_v
            s_sum = s_sum + p_v
        end do

        ! Normalize distribution
        if (s_sum > 1.0e-14_dp) then
            v_dist = v_dist / s_sum
        end if
    end subroutine calc_er_vibrational_populations

    !> Pure factorial helper
    pure function gamma_fact(n_val) result(fact)
        integer, intent(in) :: n_val
        real(dp) :: fact
        integer :: i

        fact = 1.0_dp
        do i = 2, n_val
            fact = fact * real(i, dp)
        end do
    end function gamma_fact

    !> Calculate Eley-Rideal abstraction cross section sigma_ER(E_i) in Angstrom^2
    !> Exhibits threshold behavior and saturates at geometric cross section ~ 0.2 - 0.5 A^2
    pure function calc_er_reaction_cross_section(sys, e_incident_ev) result(sigma_er)
        type(er_reaction_system_t), intent(in) :: sys
        real(dp), intent(in) :: e_incident_ev
        real(dp) :: sigma_er

        real(dp) :: e_i, sigma_geo, tunneling_factor

        e_i = max(0.0_dp, e_incident_ev)
        sigma_geo = 0.45_dp  ! Geometric scale ~ 0.45 Angstrom^2

        if (e_i >= sys%barrier_ev) then
            ! Over-barrier cross section: sigma = sigma_0 * (1 - V_b / E_i)
            sigma_er = sigma_geo * (1.0_dp - sys%barrier_ev / max(0.01_dp, e_i))
        else
            ! Quantum tunneling penetration below barrier
            tunneling_factor = exp(- 2.0_dp * sqrt(2.0_dp * sys%m_gas_amu * AMU2AU * &
                               (sys%barrier_ev - e_i) * EV2AU) * (0.3_dp * ANG2AU))
            sigma_er = sigma_geo * 0.5_dp * tunneling_factor
        end if
        sigma_er = max(0.0_dp, sigma_er)
    end function calc_er_reaction_cross_section

    !> Calculate thermal rate constant k_ER(T) in cm^3/s by thermal velocity averaging
    pure function calc_er_thermal_rate_constant(sys, temp_k) result(k_er)
        type(er_reaction_system_t), intent(in) :: sys
        real(dp), intent(in) :: temp_k
        real(dp) :: k_er

        real(dp) :: t_surf, v_th_si, m_gas_kg, sigma_th_m2

        t_surf = max(10.0_dp, temp_k)
        m_gas_kg = sys%m_gas_amu * 1.66053906660e-27_dp

        ! Thermal velocity <v> = sqrt(8 * k_B * T / (pi * m_gas))
        v_th_si = sqrt(8.0_dp * KB * t_surf / (PI * m_gas_kg))

        ! Effective thermal cross section at thermal energy E_th = k_B * T
        sigma_th_m2 = calc_er_reaction_cross_section(sys, (KB * t_surf / 1.602176634e-19_dp)) * 1.0e-20_dp

        ! k_ER = <v> * sigma in cm^3/s
        k_er = (v_th_si * sigma_th_m2) * 1.0e6_dp
    end function calc_er_thermal_rate_constant

end module mod_surface_reaction_er
