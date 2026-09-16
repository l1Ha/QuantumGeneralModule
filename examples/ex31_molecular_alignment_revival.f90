!> \file ex31_molecular_alignment_revival.f90
!> \brief Example 31: Laser-Induced Molecular Alignment Revivals & Optical Centrifuge Superrotors
!> \author LiHao
program ex31_molecular_alignment_revival
    use mod_constants, only: dp
    use mod_molecular_alignment
    implicit none

    type(rotor_molecule_t) :: n2_mol, co2_mol
    integer, parameter :: N_PTS = 100
    real(dp), dimension(N_PTS) :: t_grid_ps, cos2_trace
    integer  :: stat, it, j_super
    real(dp) :: max_align, t_max_align
    logical  :: is_dissoc
    real(dp) :: e_rot_ev

    print '(A)', "================================================================"
    print '(A)', " Example 31: Molecular Alignment Revivals & Optical Centrifuge "
    print '(A)', "================================================================"

    call init_rotor_molecule(n2_mol, "N2", b_rot_cm1=1.9982_dp, delta_alpha_au=6.70_dp, &
                             dipole_debye=0.0_dp, stat=stat)
    call init_rotor_molecule(co2_mol, "CO2", b_rot_cm1=0.3902_dp, delta_alpha_au=14.2_dp, &
                             dipole_debye=0.0_dp, stat=stat)

    print '(A)', " 1. Rotor Molecule Properties:"
    print '(A, F7.3, A, F8.3, A)', "   N2:  B_e = ", n2_mol%b_rot_cm1, " cm^-1, T_rev = ", &
        n2_mol%t_rev_ps, " ps"
    print '(A, F7.3, A, F8.3, A)', "   CO2: B_e = ", co2_mol%b_rot_cm1, " cm^-1, T_rev = ", &
        co2_mol%t_rev_ps, " ps"

    ! 2. 模拟 N2 的非绝热短脉冲对齐 (100 fs 脉冲, 10 ps 演化)
    print '(A)', " 2. Simulating N2 Laser Alignment Transient (100 fs, 3e13 W/cm^2)..."
    call simulate_laser_induced_alignment(n2_mol, laser_i0_wcm2=3.0e13_dp, pulse_fwhm_fs=100.0_dp, &
                                         temp_k=0.0_dp, j_max=14, n_time=N_PTS, t_span_ps=10.0_dp, &
                                         t_grid_ps=t_grid_ps, cos2_trace=cos2_trace, stat=stat)

    max_align = -1.0_dp
    t_max_align = 0.0_dp
    do it = 1, N_PTS
        if (cos2_trace(it) > max_align) then
            max_align = cos2_trace(it)
            t_max_align = t_grid_ps(it)
        end if
    end do

    print '(A, F8.4, A, F6.2, A)', "   Peak Alignment <cos^2 theta> = ", max_align, &
        " at t = ", t_max_align, " ps"
    print '(A, F8.4)', "   Initial Isotropic Baseline  = ", cos2_trace(1)

    ! 3. 光学离心机超转子加速
    print '(A)', " 3. Optical Centrifuge Acceleration (beta = 0.5 THz/ps, T = 100 ps):"
    call calc_optical_centrifuge_kick(n2_mol, chirp_beta_thz2=0.5_dp, pulse_duration_ps=100.0_dp, &
                                     j_superrotor=j_super, stat=stat)
    print '(A, I3, A)', "   N2 accelerated to Superrotor state:   J = ", j_super, " hbar"

    call calc_superrotor_dissociation(n2_mol, j_rot=j_super, d_e_ev=9.76_dp, &
                                      is_dissociated=is_dissoc, rot_energy_ev=e_rot_ev)
    print '(A, F7.3, A, L1)', "   Rotational Kinetic Energy: E_rot = ", e_rot_ev, &
        " eV, Centrifugal Dissociation: ", is_dissoc

    print '(A)', "================================================================"
    print '(A)', " Example 31 completed successfully.                             "
    print '(A)', "================================================================"

end program ex31_molecular_alignment_revival
