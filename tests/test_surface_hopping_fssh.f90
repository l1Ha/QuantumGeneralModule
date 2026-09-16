!> \file test_surface_hopping_fssh.f90
!> \brief Unit tests for Tully Fewest Switches Surface Hopping (FSSH)
!> \author LiHao
program test_surface_hopping_fssh
    use mod_constants, only: dp
    use mod_surface_hopping_fssh
    implicit none

    type(tully_model_t) :: model_sac, model_dac
    type(fssh_trajectory_t) :: traj
    real(dp) :: e1, e2, f1, f2, d12
    real(dp) :: e_init, e_final
    real(dp) :: t1, t2, r1, r2, total_prob
    real(dp) :: dt
    integer  :: n_pass, n_total, stat, step

    n_pass = 0
    n_total = 0

    print '(A)', "=================================================="
    print '(A)', "  GeneralModule Unit Tests: Surface Hopping FSSH  "
    print '(A)', "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Model Initialization & Avoided Crossing Gap
    ! --------------------------------------------------------------------------
    call init_tully_model(model_sac, TULLY_SAC, stat)
    call init_tully_model(model_dac, TULLY_DAC, stat)

    call calc_adiabatic_surface_and_nacv(model_sac, 0.0_dp, e1, e2, f1, f2, d12)
    n_total = n_total + 1
    ! At x = 0, gap = 2 * C = 2 * 0.005 = 0.010 a.u.
    if (abs((e2 - e1) - 0.010_dp) < 1.0e-5_dp .and. abs(d12) > 0.1_dp) then
        print '(A, F8.4, A, F8.4, A)', " [PASS] Tully SAC avoided crossing gap = ", &
            (e2 - e1), " a.u., Peak NACV d12 = ", d12, " a.u.^-1"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Tully model avoided crossing parameters incorrect"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: FSSH Trajectory Initialization & Initial Norm
    ! --------------------------------------------------------------------------
    call init_fssh_trajectory(traj, mass=2000.0_dp, x0=-5.0_dp, p0=20.0_dp, initial_state=1)
    n_total = n_total + 1
    if (traj%active_state == 1 .and. abs(abs(traj%c(1))**2 - 1.0_dp) < 1.0e-12_dp) then
        print '(A)', " [PASS] FSSH trajectory initialized on lower surface with unitary state norm"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] FSSH trajectory initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Total Energy Conservation Under Continuous Evolution
    ! --------------------------------------------------------------------------
    call calc_adiabatic_surface_and_nacv(model_sac, traj%x, e1, e2, f1, f2, d12)
    e_init = (traj%p**2) / (2.0_dp * traj%mass) + e1
    dt = 5.0_dp  ! ~0.12 fs

    do step = 1, 100
        call propagate_fssh_step(model_sac, traj, dt, rand_val=0.999_dp)  ! Force no-hop
    end do

    call calc_adiabatic_surface_and_nacv(model_sac, traj%x, e1, e2, f1, f2, d12)
    e_final = (traj%p**2) / (2.0_dp * traj%mass) + e1

    n_total = n_total + 1
    ! Energy conservation within 0.1% over 100 steps
    if (abs((e_final - e_init) / e_init) < 1.0e-3_dp) then
        print '(A, ES12.4, A, ES12.4, A)', " [PASS] Nuclear energy conserved: E_init = ", &
            e_init, " a.u., E_final = ", e_final, " a.u."
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Energy drift during FSSH propagation"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Tully 1 (SAC) Ensemble Unitarity & Branching Ratio
    ! --------------------------------------------------------------------------
    ! Run 200 trajectories at P = 20 a.u. (E_kin = 0.1 a.u. > gap)
    call run_fssh_ensemble(model_sac, initial_p=20.0_dp, initial_state=1, n_trajs=200, &
                           dt=8.0_dp, max_steps=1200, t1_prob=t1, t2_prob=t2, &
                           r1_prob=r1, r2_prob=r2, rng_seed=42)

    total_prob = t1 + t2 + r1 + r2
    n_total = n_total + 1
    if (abs(total_prob - 1.0_dp) < 1.0e-12_dp .and. t1 > 0.3_dp .and. t2 > 0.05_dp) then
        print '(A, F6.3, A, F6.3, A, F6.3, A)', " [PASS] Tully SAC ensemble branching: T1 = ", &
            t1, ", T2 = ", t2, ", Unitarity Sum = ", total_prob
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] FSSH ensemble branching or unitarity violation"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Ehrenfest Mean-Field Wave Function Normalization
    ! --------------------------------------------------------------------------
    call init_fssh_trajectory(traj, mass=2000.0_dp, x0=-4.0_dp, p0=15.0_dp, initial_state=1)
    do step = 1, 150
        call propagate_ehrenfest_step(model_sac, traj, dt=5.0_dp)
    end do

    n_total = n_total + 1
    if (abs((abs(traj%c(1))**2 + abs(traj%c(2))**2) - 1.0_dp) < 1.0e-10_dp) then
        print '(A, F10.8)', " [PASS] Ehrenfest mean-field norm conserved: |c1|^2 + |c2|^2 = ", &
            abs(traj%c(1))**2 + abs(traj%c(2))**2
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Ehrenfest norm conservation failure"
    end if

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', " Surface Hopping Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', " SUCCESS: All surface hopping FSSH tests passed."
    else
        error stop " Test failures detected in test_surface_hopping_fssh."
    end if

end program test_surface_hopping_fssh
