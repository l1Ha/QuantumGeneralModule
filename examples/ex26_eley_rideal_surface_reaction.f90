! ==============================================================================
! GeneralModule Example 26: Gas-Surface Eley-Rideal Catalytic Reaction
!
! Scientific Engineering Application:
!  1. H(g) + H(ads)/Cu(111) -> H2(g, v, j) direct Eley-Rideal reaction
!  2. Hyperthermal exothermicity partitioning (vibration ~50%, translation ~35%)
!  3. Vibrational state population distribution and vibrational inversion
!  4. Reaction cross section sigma_ER(E_i) and thermal rate constant k_ER(T)
!
! Reference:
!   Rettner, J. Chem. Phys. 101, 1529 (1994); Phys. Rev. Lett. 69, 383 (1992).
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================
program ex26_eley_rideal_surface_reaction
    use mod_constants, only: dp
    use mod_surface_reaction_er
    implicit none

    type(er_reaction_system_t)  :: er_sys
    type(er_energy_partition_t) :: ep
    real(dp) :: p_vib(0:5)
    integer  :: u_pop, u_rate, v, step
    real(dp) :: e_inc, sigma_er, temp_k, k_rate

    print *, "================================================================"
    print *, " Example 26: Surface Eley-Rideal Catalytic Reaction Dynamics    "
    print *, "================================================================"

    ! 1. Initialize H(g) + H(ads)/Cu(111) system
    call init_er_reaction_system(er_sys, "H+H/Cu(111)")

    write(*, '(A, A)')        " Reaction System             : ", trim(er_sys%name)
    write(*, '(A, F8.2, A)')  " Molecular Dissociation D_AB : ", er_sys%d_mol_ev, " eV"
    write(*, '(A, F8.2, A)')  " Chemisorption Energy D_chem : ", er_sys%d_chem_ev, " eV"
    write(*, '(A, F8.2, A)')  " Reaction Exothermicity      : ", er_sys%delta_e_exo_ev, " eV"
    print *, "----------------------------------------------------------------"

    ! 2. Energy Partitioning and Vibrational Population Inversion at E_i = 0.15 eV
    e_inc = 0.15_dp
    call calc_er_energy_partitioning(er_sys, e_incident_ev=e_inc, partition=ep)
    call calc_er_vibrational_populations(er_sys, e_incident_ev=e_inc, max_v=5, v_dist=p_vib)

    write(*, '(A, F8.2, A)')  " Incident Kinetic Energy     : ", e_inc, " eV"
    write(*, '(A, F8.2, A)')  " Available Energy E_avail    : ", ep%e_total_avail_ev, " eV"
    write(*, '(A, F8.2, A)')  " Vibrational Excitation E_vib: ", ep%e_vib_ev, " eV"
    write(*, '(A, F8.2, A)')  " Translational Energy E_trans: ", ep%e_trans_ev, " eV"
    write(*, '(A, F8.2, A)')  " Surface Bath Dissipation    : ", ep%e_diss_ev, " eV"
    print *, "----------------------------------------------------------------"

    open(newunit=u_pop, file="ex26_er_reaction_energy_vibrational.dat", &
         status="replace", action="write")
    write(u_pop, '(A)') "# GeneralModule Example 26: ER Product Vibrational Populations"
    write(u_pop, '(A)') "# System: H(g) + H(ads)/Cu(111) -> H2(g, v) + Cu(111)"
    write(u_pop, '(A)') "# Col 1: Vibrational quantum number v"
    write(u_pop, '(A)') "# Col 2: Vibrational population P(v)"

    do v = 0, 5
        write(u_pop, '(I4, ES16.6)') v, p_vib(v)
        write(*, '(A, I2, A, F8.4)') " Population P(v = ", v, ") = ", p_vib(v)
    end do
    close(u_pop)
    print *, " Saved vibrational populations to ex26_er_reaction_energy_vibrational.dat"

    ! 3. Reaction Cross Section sigma_ER(E_i) and Thermal Rate Constant k_ER(T)
    open(newunit=u_rate, file="ex26_er_cross_section_rate.dat", &
         status="replace", action="write")
    write(u_rate, '(A)') "# GeneralModule Example 26: Cross Section & Rate Constant"
    write(u_rate, '(A)') "# Col 1: Incident energy E_i (eV)"
    write(u_rate, '(A)') "# Col 2: ER reaction cross section sigma_ER (Angstrom^2)"
    write(u_rate, '(A)') "# Col 3: Substrate temperature T (K)"
    write(u_rate, '(A)') "# Col 4: ER thermal rate constant k_ER (cm^3 / s)"

    do step = 1, 100
        e_inc  = real(step, dp) * 0.01_dp             ! 0.01 to 1.0 eV
        temp_k = 100.0_dp + real(step, dp) * 9.0_dp   ! 100 to 1000 K

        sigma_er = calc_er_reaction_cross_section(er_sys, e_incident_ev=e_inc)
        k_rate   = calc_er_thermal_rate_constant(er_sys, temp_k=temp_k)

        write(u_rate, '(4ES16.6)') e_inc, sigma_er, temp_k, k_rate
    end do
    close(u_rate)
    print *, " Saved cross sections and rate constants to ex26_er_cross_section_rate.dat"
    print *, " Example 26 completed successfully."

end program ex26_eley_rideal_surface_reaction
