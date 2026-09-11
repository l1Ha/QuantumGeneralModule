!> \brief 波包时间推进与积分器单元测试
program test_propagators
    use general_module
    implicit none

    integer :: n_tests = 0
    integer :: n_passed = 0
    real(dp) :: tol = 1.0e-3_dp

    ! Bloch 测试变量
    real(dp) :: r_bloch(3), rabi, detun, dt, t, bloch_norm
    integer :: step, n_steps

    ! Split-Operator 测试变量
    integer, parameter :: nx = 128
    complex(dp) :: psi(nx)
    real(dp) :: v_pot(nx), x_grid(nx), dx, x_min, x_max, mass, norm_initial, norm_final, sigma
    integer :: i

    ! CAP 吸收边界测试变量
    type(absorbing_boundary_t) :: cap_obj
    complex(dp) :: v_cap
    real(dp) :: norm_after_cap

    print '(A)', "=================================================="
    print '(A)', "     GeneralModule Unit Tests: Propagators        "
    print '(A)', "=================================================="

    ! ----------------------------------------------------
    ! 1. 光学 Bloch 方程数值求解测试 (Rabi 振荡与布居反转)
    ! ----------------------------------------------------
    r_bloch = [0.0_dp, 0.0_dp, -1.0_dp] ! 基态 w = -1
    rabi = 1.0_dp
    detun = 0.0_dp
    dt = 0.001_dp
    n_steps = int(PI / (rabi * dt)) ! 演化至 pi 脉冲时间 t = pi/rabi

    do step = 1, n_steps
        call solve_bloch_two_level(r_bloch, rabi, detun, dt)
    end do

    ! pi 脉冲后，系统应处于激发态 w = +1.0
    call assert_close("Bloch Pi-pulse inversion w=+1", 1.0_dp, r_bloch(3), tol)

    ! 检验 Bloch 矢量模长守恒 |R| = 1
    bloch_norm = sqrt(sum(r_bloch**2))
    call assert_close("Bloch norm conservation |R|=1", 1.0_dp, bloch_norm, 1.0e-8_dp)

    ! ----------------------------------------------------
    ! 2. Split-Operator 一维波包演化范数守恒测试
    ! ----------------------------------------------------
    x_min = -15.0_dp
    x_max = 15.0_dp
    dx = (x_max - x_min) / real(nx, dp)
    mass = 1.0_dp
    sigma = 1.0_dp

    do i = 1, nx
        x_grid(i) = x_min + real(i - 1, dp) * dx
        ! 谐振子势能 V(x) = 0.5 * k * x^2
        v_pot(i) = 0.5_dp * (x_grid(i)**2)
        ! 高斯初态波包
        psi(i) = (1.0_dp / (PI * sigma**2)**0.25_dp) * exp(-0.5_dp * (x_grid(i) / sigma)**2)
    end do

    ! 初始范数
    norm_initial = sum(abs(psi)**2) * dx
    call assert_close("Initial wavepacket norm = 1", 1.0_dp, norm_initial, 1.0e-4_dp)

    ! 演化 50 步
    dt = 0.05_dp
    do step = 1, 50
        call propagate_split_operator_1d(psi, v_pot, dx, mass, dt)
    end do

    ! 结束范数检验
    norm_final = sum(abs(psi)**2) * dx
    call assert_close("Split-operator norm conservation", norm_initial, norm_final, 1.0e-5_dp)

    ! ----------------------------------------------------
    ! 3. 复吸收边界 (CAP) 测试
    ! ----------------------------------------------------
    call cap_init(r_start=10.0_dp, r_end=15.0_dp, strength=0.5_dp, cap_obj=cap_obj)

    ! 边界内无吸收
    v_cap = cap_evaluate(5.0_dp, cap_obj)
    call assert_close("CAP zero inside boundary", 0.0_dp, abs(v_cap), 1.0e-12_dp)

    ! 吸收区存在虚部势
    v_cap = cap_evaluate(12.5_dp, cap_obj)
    if (aimag(v_cap) < 0.0_dp) then
        call assert_close("CAP negative imaginary potential", 1.0_dp, 1.0_dp, 1.0e-12_dp)
    else
        call assert_close("CAP negative imaginary potential", 1.0_dp, 0.0_dp, 1.0e-12_dp)
    end if

    ! ----------------------------------------------------
    ! 4. 玻尔兹曼转动权重测试
    ! ----------------------------------------------------
    block
        real(dp) :: w_rot(0:20), z_part
        call boltzmann_rotational_weights(temp_k=300.0_dp, b_rot=1.0e-5_dp, j_max=20, weights=w_rot, z_rot=z_part)
        call assert_close("Boltzmann weights partition sum = 1", 1.0_dp, sum(w_rot), 1.0e-10_dp)

        ! T = 0 K 纯基态
        call boltzmann_rotational_weights(temp_k=0.0_dp, b_rot=1.0e-5_dp, j_max=20, weights=w_rot)
        call assert_close("Boltzmann ground state at 0K", 1.0_dp, w_rot(0), 1.0e-12_dp)
    end block

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Propagator Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some tests did not pass!"
        stop 1
    else
        print '(A)', "SUCCESS: All propagator tests passed."
    end if

contains

    subroutine assert_close(name, expected, actual, eps)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, eps
        real(dp) :: diff

        n_tests = n_tests + 1
        diff = abs(expected - actual)
        if (diff <= eps) then
            n_passed = n_passed + 1
            print '(A, A35, A)', " [PASS] ", name, ""
        else
            print '(A, A35, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_propagators
