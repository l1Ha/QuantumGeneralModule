! ==============================================================================
! GeneralModule Example 24: Rydberg Blockade Dynamics & Many-Body Scars
!
! Features:
!  1. 87Rb 70S Rydberg state C6 van der Waals interaction and blockade radius Rb
!  2. Two-atom quantum dynamics: unblockaded (R > Rb) vs blockaded (R < Rb)
!  3. Collective entangled W-state excitation and double-excitation suppression
!  4. 1D chain quantum many-body scar oscillations of Z2 staggered order
!
! Standard: Fortran 2008
! ==============================================================================
program ex24_rydberg_blockade_dynamics
    use mod_constants, only: dp
    use mod_rydberg_blockade
    implicit none

    type(rydberg_atom_t) :: rb70
    type(rydberg_array_config_t) :: chain_cfg
    real(dp) :: r_b, rabi_mhz, det_mhz
    integer  :: u_two, u_scar, i
    integer, parameter :: N_T = 150
    real(dp) :: t_arr(N_T), pg_blk(N_T), ps_blk(N_T), pd_blk(N_T)
    real(dp) :: pg_unblk(N_T), ps_unblk(N_T), pd_unblk(N_T)
    real(dp) :: t_scar(N_T), z2_scar(N_T)

    print *, "================================================================"
    print *, " Example 24: Rydberg Blockade Dynamics & Many-Body Scars        "
    print *, "================================================================"

    ! 1. Initialize 87Rb 70S Rydberg state
    call init_rydberg_atom(rb70, "87Rb", n_principal=70, l_orbital=0)

    rabi_mhz = 2.0_dp
    det_mhz  = 0.0_dp
    r_b = calc_rydberg_blockade_radius(rb70, rabi_mhz)

    write(*, '(A, F10.1, A)') " van der Waals C6 / h        : ", rb70%c6_mhz_um6, " MHz*um^6"
    write(*, '(A, F10.2, A)') " Rydberg Radiative Lifetime  : ", rb70%lifetime_us, " us"
    write(*, '(A, F10.2, A)') " Laser Rabi Frequency Omega  : ", rabi_mhz, " MHz"
    write(*, '(A, F10.2, A)') " Rydberg Blockade Radius R_b : ", r_b, " um"
    print *, "----------------------------------------------------------------"

    ! 2. Two-Atom Dynamics: Blockaded (R = 4.0 um < R_b) vs Unblockaded (R = 15.0 um > R_b)
    call calc_two_atom_dynamics(rb70, spacing_um=4.0_dp, rabi_mhz=rabi_mhz, detuning_mhz=det_mhz, &
                                t_max_us=1.5_dp, n_steps=N_T, t_arr=t_arr, &
                                p_g=pg_blk, p_single=ps_blk, p_double=pd_blk)

    call calc_two_atom_dynamics(rb70, spacing_um=15.0_dp, rabi_mhz=rabi_mhz, detuning_mhz=det_mhz, &
                                t_max_us=1.5_dp, n_steps=N_T, t_arr=t_arr, &
                                p_g=pg_unblk, p_single=ps_unblk, p_double=pd_unblk)

    open(newunit=u_two, file="ex24_two_atom_blockade.dat", status="replace", action="write")
    write(u_two, '(A)') "# GeneralModule Example 24: Two-Atom Rydberg Blockade Dynamics"
    write(u_two, '(A)') "# Col 1: Time t (us)"
    write(u_two, '(A)') "# Col 2: P_gg (Blockaded, R=4 um)"
    write(u_two, '(A)') "# Col 3: P_W  (Blockaded, R=4 um)"
    write(u_two, '(A)') "# Col 4: P_rr (Blockaded, R=4 um)"
    write(u_two, '(A)') "# Col 5: P_gg (Unblockaded, R=15 um)"
    write(u_two, '(A)') "# Col 6: P_W  (Unblockaded, R=15 um)"
    write(u_two, '(A)') "# Col 7: P_rr (Unblockaded, R=15 um)"

    write(*, '(A)') "  t (us)  | P_W (Blk, 4um) | P_rr (Blk, 4um) | P_W (Unblk,15um) | P_rr (Unblk,15um)"
    write(*, '(A)') "----------+----------------+-----------------+------------------+------------------"

    do i = 1, N_T
        write(u_two, '(7ES16.6)') t_arr(i), pg_blk(i), ps_blk(i), pd_blk(i), &
                                  pg_unblk(i), ps_unblk(i), pd_unblk(i)
        if (mod(i, 20) == 1) then
            write(*, '(F9.3, " | ", F14.4, " | ", ES15.4, " | ", F16.4, " | ", F17.4)') &
                t_arr(i), ps_blk(i), pd_blk(i), ps_unblk(i), pd_unblk(i)
        end if
    end do
    close(u_two)
    print *, " Two-atom dynamics saved to ex24_two_atom_blockade.dat"
    print *, "----------------------------------------------------------------"

    ! 3. Many-Body Quantum Scars on a 1D Blockaded Chain
    call init_rydberg_array(chain_cfg, n_atoms=10, spacing_um=5.0_dp, &
                            rabi_mhz=rabi_mhz, detuning_mhz=det_mhz, boundary_cond=2)
    call calc_rydberg_scar_dynamics(chain_cfg, rb70, t_max_us=2.5_dp, n_steps=N_T, &
                                    t_arr=t_scar, z2_order_arr=z2_scar)

    open(newunit=u_scar, file="ex24_rydberg_scar_dynamics.dat", status="replace", action="write")
    write(u_scar, '(A)') "# GeneralModule Example 24: Many-Body Quantum Scar Dynamics"
    write(u_scar, '(A)') "# Col 1: Time t (us)"
    write(u_scar, '(A)') "# Col 2: Staggered Z2 Order Parameter O_Z2(t)"

    do i = 1, N_T
        write(u_scar, '(2ES16.6)') t_scar(i), z2_scar(i)
    end do
    close(u_scar)

    print *, " SUCCESS: Scar dynamics saved to ex24_rydberg_scar_dynamics.dat"
    print *, "================================================================"

end program ex24_rydberg_blockade_dynamics
