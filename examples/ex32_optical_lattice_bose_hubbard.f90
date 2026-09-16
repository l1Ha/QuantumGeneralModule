!> \file ex32_optical_lattice_bose_hubbard.f90
!> \brief Example 32: Optical Lattice Bloch Bands, Bose-Hubbard Mapping & Bloch Oscillations
!> \author LiHao
program ex32_optical_lattice_bose_hubbard
    use mod_constants, only: dp
    use mod_optical_lattice_hubbard
    implicit none

    type(optical_lattice_t) :: latt
    type(bose_hubbard_param_t) :: bh
    integer  :: stat, is
    real(dp) :: s_val, t_bloch, omega_bloch, p_lz
    real(dp) :: force_grav
    real(dp), dimension(6) :: s_grid

    print '(A)', "================================================================"
    print '(A)', " Example 32: Ultracold Optical Lattice & Bose-Hubbard Dynamics  "
    print '(A)', "================================================================"

    print '(A)', " 1. 87Rb Optical Lattice Parameter Scan (lambda = 1064 nm):"
    print '(A)', "    Depth s (E_R) | J/h (Hz)   | U/h (Hz)   | U/J Ratio | Phase Regime"
    print '(A)', "   --------------------------------------------------------------"

    s_grid = [2.0_dp, 4.0_dp, 8.0_dp, 12.0_dp, 16.0_dp, 20.0_dp]

    do is = 1, 6
        s_val = s_grid(is)
        call init_optical_lattice(latt, mass_amu=86.909_dp, lambda_nm=1064.0_dp, s_depth=s_val, stat=stat)
        call calc_bose_hubbard_parameters(latt, a_s_bohr=100.0_dp, omega_perp_hz=1500.0_dp, bh=bh, stat=stat)

        if (bh%is_mott_candidate) then
            print '(4X, F9.1, 4X, F9.2, 4X, F9.2, 4X, F8.2, 4X, A)', &
                s_val, bh%j_hopping_hz, bh%u_onsite_hz, bh%u_over_j_ratio, "Mott Insulator"
        else
            print '(4X, F9.1, 4X, F9.2, 4X, F9.2, 4X, F8.2, 4X, A)', &
                s_val, bh%j_hopping_hz, bh%u_onsite_hz, bh%u_over_j_ratio, "Superfluid"
        end if
    end do

    ! 2. 重力驱动的原子布洛赫振荡
    force_grav = (86.909_dp * 1.66053906660e-27_dp) * 9.80665_dp
    call init_optical_lattice(latt, mass_amu=86.909_dp, lambda_nm=1064.0_dp, s_depth=10.0_dp, stat=stat)
    call calc_bloch_oscillation_dynamics(latt, force_grav, t_bloch_ms=t_bloch, &
                                        omega_bloch_hz=omega_bloch, p_lz_tunnel=p_lz, stat=stat)

    print '(A)', " 2. Vertical Gravitational Bloch Oscillations (s = 10 E_R):"
    print '(A, F7.3, A, F8.1, A)', "   Oscillation Period T_B = ", t_bloch, " ms, Frequency nu_B = ", &
        (1.0e3_dp / t_bloch), " Hz"
    print '(A, ES14.4)', "   Landau-Zener Interband Tunneling P_LZ = ", p_lz

    print '(A)', "================================================================"
    print '(A)', " Example 32 completed successfully.                             "
    print '(A)', "================================================================"

end program ex32_optical_lattice_bose_hubbard
