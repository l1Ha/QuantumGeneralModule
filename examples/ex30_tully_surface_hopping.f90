!> \file ex30_tully_surface_hopping.f90
!> \brief Example 30: Tully Fewest Switches Surface Hopping (FSSH) Non-Adiabatic Dynamics
!> \author LiHao
program ex30_tully_surface_hopping
    use mod_constants, only: dp
    use mod_surface_hopping_fssh
    implicit none

    type(tully_model_t) :: model_sac, model_dac
    real(dp) :: p_in, t1, t2, r1, r2
    integer  :: stat, ip
    real(dp), dimension(4) :: p_grid

    print '(A)', "================================================================"
    print '(A)', " Example 30: Tully Fewest Switches Surface Hopping (FSSH)       "
    print '(A)', "================================================================"

    call init_tully_model(model_sac, TULLY_SAC, stat)
    call init_tully_model(model_dac, TULLY_DAC, stat)

    print '(A)', " 1. Tully Model 1 (Simple Avoided Crossing - SAC):"
    print '(A)', "    Momentum P (a.u.) | E_kin (eV) | Trans T1  | Trans T2  | Reflect R"
    print '(A)', "   --------------------------------------------------------------"

    p_grid = [10.0_dp, 15.0_dp, 20.0_dp, 25.0_dp]

    do ip = 1, 4
        p_in = p_grid(ip)
        call run_fssh_ensemble(model_sac, initial_p=p_in, initial_state=1, n_trajs=150, &
                               dt=6.0_dp, max_steps=1200, t1_prob=t1, t2_prob=t2, &
                               r1_prob=r1, r2_prob=r2, rng_seed=100 + ip)
        print '(4X, F12.1, 4X, F9.3, 4X, F8.3, 4X, F8.3, 4X, F8.3)', &
            p_in, (p_in**2 / (2.0_dp * 2000.0_dp)) * 27.211386_dp, t1, t2, (r1 + r2)
    end do

    print '(A)', " 2. Tully Model 2 (Dual Avoided Crossing - DAC Stueckelberg Interference):"
    p_in = 20.0_dp
    call run_fssh_ensemble(model_dac, initial_p=p_in, initial_state=1, n_trajs=150, &
                           dt=6.0_dp, max_steps=1400, t1_prob=t1, t2_prob=t2, &
                           r1_prob=r1, r2_prob=r2, rng_seed=202)
    print '(A, F6.1, A, F6.3, A, F6.3, A)', "    At P = ", p_in, " a.u.: T1 = ", t1, &
        ", T2 = ", t2, " (Stueckelberg oscillation region)"

    print '(A)', "================================================================"
    print '(A)', " Example 30 completed successfully.                             "
    print '(A)', "================================================================"

end program ex30_tully_surface_hopping
