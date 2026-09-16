!> \file test_reaction_path_hamiltonian.f90
!> \brief Unit tests for Reaction Path Hamiltonian (RPH), CVT, and Eckart tunneling
!> \author LiHao
program test_reaction_path_hamiltonian
    use mod_constants, only: dp
    use mod_reaction_path_hamiltonian
    implicit none

    type(rph_path_t) :: path
    real(dp) :: rate_ts, rate_cvt, s_opt
    real(dp) :: kappa_300, kappa_1000
    integer  :: n_pass, n_total, stat, ts_idx

    n_pass = 0
    n_total = 0

    print '(A)', "=================================================="
    print '(A)', "  GeneralModule Unit Tests: RPH & Variational TST "
    print '(A)', "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Benchmark Reaction Path & Saddle Point (H + H2)
    ! --------------------------------------------------------------------------
    call init_rph_benchmark_reaction(path, reaction_type=1, stat=stat)
    n_total = n_total + 1
    ts_idx = (path%n_points + 1) / 2  ! s = 0 point

    if (stat == 0 .and. abs(path%points(ts_idx)%v_pot_ev - 0.42_dp) < 1.0e-3_dp .and. &
        abs(path%points(ts_idx)%s_coord) < 1.0e-12_dp) then
        print '(A, F6.3, A, F7.1, A)', " [PASS] RPH path initialized: V_barrier = ", &
            path%barrier_height_ev, " eV, imag freq = ", path%imag_freq_ts_cm1, " cm^-1"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] RPH path initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Generalized TST Rate Temperature Dependence
    ! --------------------------------------------------------------------------
    call calc_generalized_tst_rate(path, s_idx=ts_idx, temp_k=300.0_dp, rate_gtst=rate_ts)
    n_total = n_total + 1
    ! At 300K, exp(-0.42 eV / kB*T) ~ 8.7e-8, rate ~ 1e5 - 1e7 s^-1
    if (rate_ts > 1.0e4_dp .and. rate_ts < 1.0e8_dp) then
        print '(A, ES12.4, A)', " [PASS] Conventional TST rate at 300K: k_TST = ", rate_ts, " s^-1"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] TST rate out of physical range"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Variational Principle (k_CVT <= k_TST at saddle point)
    ! --------------------------------------------------------------------------
    call calc_cvt_rate_constant(path, temp_k=300.0_dp, s_opt=s_opt, rate_cvt=rate_cvt, stat=stat)
    n_total = n_total + 1
    ! By variational definition: CVT rate <= conventional TST rate at s=0
    if (stat == 0 .and. rate_cvt <= rate_ts * 1.000001_dp .and. rate_cvt > 0.0_dp) then
        print '(A, ES12.4, A, F6.3, A)', " [PASS] Variational CVT bound verified: k_CVT = ", &
            rate_cvt, " s^-1 at s* = ", s_opt, " a0*amu^1/2"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Variational transition state inequality violated"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Eckart Tunneling Enhancement at 300 K
    ! --------------------------------------------------------------------------
    call calc_eckart_tunneling_factor(barrier_height_ev=0.42_dp, imag_freq_cm1=1500.0_dp, &
                                      temp_k=300.0_dp, kappa_tunnel=kappa_300)
    n_total = n_total + 1
    ! For H-atom transfer with 1500 cm^-1 barrier at 300K, kappa > 2.5
    if (kappa_300 > 2.0_dp .and. kappa_300 < 50.0_dp) then
        print '(A, F6.2)', " [PASS] Quantum tunneling factor at 300K: kappa(300K) = ", kappa_300
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Eckart tunneling factor out of physical range"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: High-Temperature Classical Limit (kappa -> 1 as T -> inf)
    ! --------------------------------------------------------------------------
    call calc_eckart_tunneling_factor(barrier_height_ev=0.42_dp, imag_freq_cm1=1500.0_dp, &
                                      temp_k=1500.0_dp, kappa_tunnel=kappa_1000)
    n_total = n_total + 1
    ! At 1500 K, u_TS << 1, kappa should be close to 1
    if (kappa_1000 >= 1.0_dp .and. kappa_1000 < kappa_300 .and. kappa_1000 < 1.30_dp) then
        print '(A, F6.3, A)', " [PASS] Classical high-temperature limit: kappa(1500K) = ", &
            kappa_1000, " (near classical 1.0)"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Classical high-temperature tunneling limit error"
    end if

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', " RPH & VTST Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', " SUCCESS: All Reaction Path Hamiltonian & VTST tests passed."
    else
        error stop " Test failures detected in test_reaction_path_hamiltonian."
    end if

end program test_reaction_path_hamiltonian
