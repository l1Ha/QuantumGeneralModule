!> \file ex33_rph_variational_transition_state.f90
!> \brief Example 33: Reaction Path Hamiltonian (RPH), CVT, and Eckart Tunneling Dynamics
!> \author LiHao
program ex33_rph_variational_transition_state
    use mod_constants, only: dp
    use mod_reaction_path_hamiltonian
    implicit none

    type(rph_path_t) :: path
    integer  :: stat, it, ts_idx
    real(dp) :: t_curr, s_opt, rate_tst, rate_cvt, kappa_tunnel, rate_final
    real(dp), dimension(5) :: t_grid

    print '(A)', "================================================================"
    print '(A)', " Example 33: Reaction Path Hamiltonian & Variational TST (CVT)  "
    print '(A)', "================================================================"

    call init_rph_benchmark_reaction(path, reaction_type=1, stat=stat)
    ts_idx = (path%n_points + 1) / 2

    print '(A)', " 1. Reaction Path Properties:"
    print '(A, F6.3, A, F7.1, A)', "   Barrier Height V_0 = ", path%barrier_height_ev, &
        " eV, Imaginary Frequency = ", path%imag_freq_ts_cm1, " cm^-1"
    print '(A, I2, A)', "   Number of Generalized Orthogonal Normal Modes = ", &
        path%points(1)%n_modes, " modes"

    print '(A)', " 2. Temperature Dependence & Quantum Tunneling Enhancement:"
    print '(A)', "    T (K)   | k_TST (s^-1) | k_CVT (s^-1) | kappa(T) | k_total (s^-1)"
    print '(A)', "   ------------------------------------------------------------------"

    t_grid = [200.0_dp, 300.0_dp, 400.0_dp, 600.0_dp, 1000.0_dp]

    do it = 1, 5
        t_curr = t_grid(it)
        call calc_generalized_tst_rate(path, ts_idx, t_curr, rate_tst)
        call calc_cvt_rate_constant(path, t_curr, s_opt, rate_cvt, stat)
        call calc_eckart_tunneling_factor(path%barrier_height_ev, path%imag_freq_ts_cm1, &
                                          t_curr, kappa_tunnel)
        rate_final = kappa_tunnel * rate_cvt

        print '(4X, F6.1, 4X, ES11.4, 4X, ES11.4, 4X, F7.2, 4X, ES11.4)', &
            t_curr, rate_tst, rate_cvt, kappa_tunnel, rate_final
    end do

    print '(A)', "================================================================"
    print '(A)', " Example 33 completed successfully.                             "
    print '(A)', "================================================================"

end program ex33_rph_variational_transition_state
