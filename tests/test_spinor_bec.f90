! ==============================================================================
! GeneralModule: test_spinor_bec.f90
!
! Unit Test Suite: Ultracold Spinor Bose-Einstein Condensate Dynamics
!
! Standard: Fortran 2008
! ==============================================================================

program test_spinor_bec
    use mod_constants, only: dp
    use mod_spinor_bec
    implicit none

    integer :: n_pass, n_total
    integer :: stat, step
    type(spinor_param_t) :: rb87_param, na23_param
    type(spinor_state_t) :: state
    real(dp) :: qz_1g, qz_2g, init_norm, init_mz, dt_au

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Spinor BEC Dynamics    "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Spinor Presets & Ferromagnetic / Polar Classification
    ! --------------------------------------------------------------------------
    call init_spinor_preset("87Rb", b_field_gauss=1.0_dp, density_cm3=1.0e14_dp, &
                            param=rb87_param, stat=stat)

    n_total = n_total + 1
    if (stat == 0 .and. rb87_param%is_ferromagnet .and. rb87_param%c2_au < 0.0_dp) then
        print *, " [PASS] 87Rb correctly identified as Ferromagnetic (c_2 < 0): c_2 = ", &
                 rb87_param%c2_au, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] 87Rb ferromagnetic classification failed"
    end if

    call init_spinor_preset("23Na", b_field_gauss=1.0_dp, density_cm3=1.0e14_dp, &
                            param=na23_param, stat=stat)

    n_total = n_total + 1
    if (stat == 0 .and. (.not. na23_param%is_ferromagnet) .and. na23_param%c2_au > 0.0_dp) then
        print *, " [PASS] 23Na correctly identified as Antiferromagnetic/Polar (c_2 > 0): c_2 = ", &
                 na23_param%c2_au, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] 23Na polar classification failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Quadratic Zeeman Effect (QZE) B^2 Scaling
    ! --------------------------------------------------------------------------
    qz_1g = calc_quadratic_zeeman_shift(1.0_dp, 6.83468_dp)
    qz_2g = calc_quadratic_zeeman_shift(2.0_dp, 6.83468_dp)

    n_total = n_total + 1
    if (abs(qz_2g / qz_1g - 4.0_dp) < 1.0e-5_dp) then
        print *, " [PASS] Quadratic Zeeman shift scales exactly as B^2: q(2G)/q(1G) = ", &
                 qz_2g / qz_1g, " ~ 4.0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] QZE B^2 scaling violation: ", qz_2g / qz_1g
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Total Norm and Magnetization Conservation in RK4 Propagation
    ! --------------------------------------------------------------------------
    ! Initial state mostly in m=0 with tiny seed in m=+/-1:
    state%zeta_0  = cmplx(sqrt(0.98_dp), 0.0_dp, kind=dp)
    state%zeta_p1 = cmplx(sqrt(0.01_dp), 0.0_dp, kind=dp)
    state%zeta_m1 = cmplx(sqrt(0.01_dp), 0.0_dp, kind=dp)

    state%pop_0  = abs(state%zeta_0)**2
    state%pop_p1 = abs(state%zeta_p1)**2
    state%pop_m1 = abs(state%zeta_m1)**2
    state%total_norm = state%pop_0 + state%pop_p1 + state%pop_m1
    state%magnetization_mz = state%pop_p1 - state%pop_m1

    init_norm = state%total_norm
    init_mz   = state%magnetization_mz

    ! Propagate 1000 RK4 steps
    dt_au = 1.0e4_dp  ! ~ 0.24 ps per step
    do step = 1, 1000
        call propagate_spinor_sma_rk4(rb87_param, dt_au, state)
    end do

    n_total = n_total + 1
    if (abs(state%total_norm - init_norm) < 1.0e-12_dp) then
        print *, " [PASS] Total spinor norm strictly conserved: Norm = ", state%total_norm, &
                 " (Err < 1e-12)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Norm conservation violated: ", state%total_norm
    end if

    n_total = n_total + 1
    if (abs(state%magnetization_mz - init_mz) < 1.0e-12_dp) then
        print *, " [PASS] Magnetization m_z strictly conserved: m_z = ", state%magnetization_mz, &
                 " (Err < 1e-12)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Magnetization m_z conservation violated: ", state%magnetization_mz
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Coherent Spin Mixing Population Dynamics
    ! --------------------------------------------------------------------------
    ! Since c2 < 0, population should flow from m=0 into m=+/-1 states
    n_total = n_total + 1
    if (state%pop_p1 > 0.01_dp .and. state%pop_m1 > 0.01_dp) then
        print *, " [PASS] Spin mixing dynamics active: P_+1 grew from 0.01 to ", state%pop_p1
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Spin mixing inactive"
    end if

    print *, "--------------------------------------------------"
    print *, "Spinor BEC Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All spinor BEC dynamics tests passed."
    else
        stop 1
    end if

end program test_spinor_bec
