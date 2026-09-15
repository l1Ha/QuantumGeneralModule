! ==============================================================================
! GeneralModule Example 23: Ultracold Polar Molecule Reaction Shielding
!
! Features:
!  1. 40K-87Rb ultracold polar molecule preset
!  2. Microwave blue-detuned engineered repulsive shielding potential barrier
!  3. Quantum Defect Theory (QDT) and WKB quantum tunneling suppression
!  4. Elastic rate K2_el, inelastic reaction rate K2_inel, and evaporative
!     cooling figure of merit gamma = K2_el / K2_inel
!
! Standard: Fortran 2008
! ==============================================================================
program ex23_ultracold_molecule_shielding
    use mod_constants, only: dp
    use mod_ultracold_reaction_shielding
    implicit none

    type(ultracold_molecule_t) :: krb
    type(shielding_config_t)   :: cfg
    integer  :: u_pot, u_rates, ir, id
    integer, parameter :: N_R = 150, N_DET = 25
    real(dp) :: r_bohr, v_eff_k, r_bar, v_bar_k
    real(dp) :: det_mhz, k2_el, k2_inel, gamma_val
    real(dp) :: temp_uk

    print *, "================================================================"
    print *, " Example 23: Ultracold Polar Molecule Reaction Shielding        "
    print *, "================================================================"

    ! 1. Initialize KRb molecule
    call init_ultracold_molecule_preset(krb, "KRb")

    write(*, '(A, F8.2, A)') " Polar Molecule Species           : 40K-87Rb"
    write(*, '(A, F8.2, A)') " Molecular Mass                   : ", krb%mass_amu, " amu"
    write(*, '(A, F8.3, A)') " Permanent Electric Dipole Moment : ", krb%dipole_debye, " Debye"
    write(*, '(A, ES10.2, A)')" van der Waals Dispersion C6      : ", krb%c6_au, " a.u."
    print *, "----------------------------------------------------------------"

    ! 2. Configure Microwave Shielding
    call init_shielding_config(cfg, method=1, detuning_mhz=15.0_dp, rabi_mhz=5.0_dp, &
                              e_field_kv_cm=0.0_dp, y_loss=1.0_dp)
    call calc_shielding_barrier_height(krb, cfg, r_bar, v_bar_k)

    write(*, '(A, F8.1, A)') " MW Blue Detuning Delta           : ", cfg%detuning_mhz, " MHz"
    write(*, '(A, F8.1, A)') " MW Rabi Frequency Omega          : ", cfg%rabi_mhz, " MHz"
    write(*, '(A, F8.1, A)') " Shielding Barrier Location R_bar : ", r_bar, " a0"
    write(*, '(A, F8.2, A)') " Barrier Height V_bar             : ", v_bar_k * 1.0e6_dp, " uK"
    print *, "----------------------------------------------------------------"

    ! 3. Output Effective Potential Curve V_eff(R)
    open(newunit=u_pot, file="ex23_shielding_potential.dat", status="replace", action="write")
    write(u_pot, '(A)') "# GeneralModule Example 23: Shielded Intermolecular Potential V_eff(R)"
    write(u_pot, '(A)') "# Col 1: Distance R (Bohr)"
    write(u_pot, '(A)') "# Col 2: Effective Potential V_eff (uK)"

    do ir = 1, N_R
        r_bohr = 50.0_dp + real(ir - 1, dp) * 15.0_dp
        call calc_effective_shielding_potential(krb, cfg, r_bohr, 0, v_eff_k)
        write(u_pot, '(2ES16.6)') r_bohr, v_eff_k * 1.0e6_dp
    end do
    close(u_pot)
    print *, " Potential curve saved to ex23_shielding_potential.dat"
    print *, "----------------------------------------------------------------"

    ! 4. Scan Microwave Detuning Delta at T = 0.5 uK
    temp_uk = 0.5_dp
    open(newunit=u_rates, file="ex23_shielding_rates.dat", status="replace", action="write")
    write(u_rates, '(A)') "# GeneralModule Example 23: Shielded Collision Rates vs Detuning"
    write(u_rates, '(A)') "# Col 1: MW Detuning Delta (MHz)"
    write(u_rates, '(A)') "# Col 2: Elastic Rate K2_el (cm^3/s)"
    write(u_rates, '(A)') "# Col 3: Inelastic Reaction Rate K2_inel (cm^3/s)"
    write(u_rates, '(A)') "# Col 4: Good-to-Bad Ratio gamma = K2_el / K2_inel"

    write(*, '(A)') " Delta (MHz) | K2_el (cm^3/s) | K2_inel (cm^3/s) |    gamma (Ratio)"
    write(*, '(A)') "-------------+----------------+------------------+------------------"

    do id = 1, N_DET
        det_mhz = 2.0_dp + real(id - 1, dp) * 1.2_dp
        cfg%detuning_mhz = det_mhz
        call calc_shielded_scattering_rates(krb, cfg, temp_uk, k2_el, k2_inel, gamma_val)

        write(u_rates, '(4ES16.6)') det_mhz, k2_el, k2_inel, gamma_val
        if (mod(id, 4) == 1) then
            write(*, '(F12.2, " | ", ES14.4, " | ", ES16.4, " | ", ES17.4)') &
                det_mhz, k2_el, k2_inel, gamma_val
        end if
    end do
    close(u_rates)

    print *, "----------------------------------------------------------------"
    print *, " SUCCESS: Shielded collision rates saved to ex23_shielding_rates.dat"
    print *, "================================================================"

end program ex23_ultracold_molecule_shielding
