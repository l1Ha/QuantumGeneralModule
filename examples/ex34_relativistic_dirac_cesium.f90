!> \file ex34_relativistic_dirac_cesium.f90
!> \brief Example 34: Relativistic Dirac Atomic Structure, Fine Structure Splitting, and D1/D2 Lines
!> \author LiHao
program ex34_relativistic_dirac_cesium
    use mod_constants, only: dp
    use mod_relativistic_atomic
    implicit none

    type(dirac_state_t) :: s_6s12, p_6p12, p_6p32, d_5d32, d_5d52
    real(dp) :: delta_fs_ev, delta_fs_cm1
    real(dp) :: f_d1, f_d2, ratio_d2_d1
    real(dp), parameter :: Z_CS = 55.0_dp
    real(dp), parameter :: Z_ION = 1.0_dp
    real(dp), parameter :: ALPHA_CORE_CS = 19.0_dp
    real(dp), parameter :: R_CUT_CS = 2.0_dp

    print '(A)', "================================================================"
    print '(A)', " Example 34: Relativistic Dirac Atomic Structure (Cesium Cs I)  "
    print '(A)', "================================================================"

    ! 1. Solve relativistic bound states
    call solve_radial_dirac_eigenvalue(Z_CS, Z_ION, ALPHA_CORE_CS, R_CUT_CS, &
                                       n_princ=6, kappa=-1, state=s_6s12)
    call solve_radial_dirac_eigenvalue(Z_CS, Z_ION, ALPHA_CORE_CS, R_CUT_CS, &
                                       n_princ=6, kappa=1, state=p_6p12)
    call solve_radial_dirac_eigenvalue(Z_CS, Z_ION, ALPHA_CORE_CS, R_CUT_CS, &
                                       n_princ=6, kappa=-2, state=p_6p32)
    call solve_radial_dirac_eigenvalue(Z_CS, Z_ION, ALPHA_CORE_CS, R_CUT_CS, &
                                       n_princ=5, kappa=2, state=d_5d32)
    call solve_radial_dirac_eigenvalue(Z_CS, Z_ION, ALPHA_CORE_CS, R_CUT_CS, &
                                       n_princ=5, kappa=-3, state=d_5d52)

    print '(A)', " 1. Calculated Relativistic Dirac Valence Energies:"
    print '(A)', "    State  | kappa |  l  | 2*j | Energy (a.u.) | Energy (eV)  | Defect mu"
    print '(A)', "   ----------------------------------------------------------------------"
    print '(4X, A6, 4X, I2, 4X, I1, 4X, I1, 4X, F11.5, 4X, F9.4, 4X, F7.3)', &
        "6s_1/2", s_6s12%kappa, s_6s12%l_orb, s_6s12%two_j, s_6s12%energy_au, &
        s_6s12%energy_ev, s_6s12%quantum_defect
    print '(4X, A6, 4X, I2, 4X, I1, 4X, I1, 4X, F11.5, 4X, F9.4, 4X, F7.3)', &
        "6p_1/2", p_6p12%kappa, p_6p12%l_orb, p_6p12%two_j, p_6p12%energy_au, &
        p_6p12%energy_ev, p_6p12%quantum_defect
    print '(4X, A6, 4X, I2, 4X, I1, 4X, I1, 4X, F11.5, 4X, F9.4, 4X, F7.3)', &
        "6p_3/2", p_6p32%kappa, p_6p32%l_orb, p_6p32%two_j, p_6p32%energy_au, &
        p_6p32%energy_ev, p_6p32%quantum_defect
    print '(4X, A6, 4X, I2, 4X, I1, 4X, I1, 4X, F11.5, 4X, F9.4, 4X, F7.3)', &
        "5d_3/2", d_5d32%kappa, d_5d32%l_orb, d_5d32%two_j, d_5d32%energy_au, &
        d_5d32%energy_ev, d_5d32%quantum_defect
    print '(4X, A6, 4X, I2, 4X, I1, 4X, I1, 4X, F11.5, 4X, F9.4, 4X, F7.3)', &
        "5d_5/2", d_5d52%kappa, d_5d52%l_orb, d_5d52%two_j, d_5d52%energy_au, &
        d_5d52%energy_ev, d_5d52%quantum_defect

    ! 2. Fine-Structure Splitting
    call calc_dirac_fine_structure_splitting(p_6p12, p_6p32, delta_fs_ev, delta_fs_cm1)
    print '(A)', ""
    print '(A)', " 2. Fine Structure Splitting (6p Doublet):"
    print '(A, ES12.5, A, F9.3, A)', "    Delta E(6p_3/2 - 6p_1/2) = ", delta_fs_ev, &
        " eV (", delta_fs_cm1, " cm^-1)"

    ! 3. Relativistic E1 Transition Oscillator Strengths (Cesium D1 & D2 lines)
    ! Radial overlap for 6s -> 6p is typically ~ 4.5 a.u. in alkali atoms
    call calc_dirac_e1_matrix_element(s_6s12, p_6p12, r_overlap_au=3.80_dp, osc_strength=f_d1)
    call calc_dirac_e1_matrix_element(s_6s12, p_6p32, r_overlap_au=3.80_dp, osc_strength=f_d2)
    ratio_d2_d1 = f_d2 / max(1.0e-15_dp, f_d1)

    print '(A)', ""
    print '(A)', " 3. Relativistic Dipole Oscillator Strengths (Cesium D1 & D2 Lines):"
    print '(A, F8.4)', "    f(6s_1/2 -> 6p_1/2) [Cs D1] = ", f_d1
    print '(A, F8.4)', "    f(6s_1/2 -> 6p_3/2) [Cs D2] = ", f_d2
    print '(A, F8.4, A)', "    Multiplet branching ratio f(D2) / f(D1) = ", ratio_d2_d1, &
        " (expected ~ 1.0 - 2.0)"

    print '(A)', "================================================================"
    print '(A)', " Example 34 completed successfully.                             "
    print '(A)', "================================================================"

end program ex34_relativistic_dirac_cesium
