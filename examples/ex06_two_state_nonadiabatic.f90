!> \brief 示例 6: 双态非绝热势能面避差穿越动力学与分支比模拟
!> \details 模拟核波包穿过经典 Tully 避差交叉点（Avoided Crossing），计算无绝热跃迁与各通道末态布居。
program ex06_two_state_nonadiabatic
    use general_module
    implicit none

    integer, parameter :: nx = 512
    integer, parameter :: n_steps = 1200
    real(dp) :: x_min, x_max, dx, mass, dt, t
    real(dp) :: x_0, p_0, sigma_x, velocity
    real(dp) :: v11(nx), v22(nx), v12(nx), x_grid(nx)
    complex(dp) :: psi1(nx), psi2(nx)
    real(dp) :: pop1, pop2, ratio, p_lz, delta_slope
    integer :: i, step, file_unit

    print '(A)', "=========================================================="
    print '(A)', "  Example 06: Two-State Non-Adiabatic Avoided Crossing    "
    print '(A)', "=========================================================="

    ! 1. 空间网格与物理模型参数 (Tully Simple Avoided Crossing 模型, a.u.)
    x_min = -25.0_dp
    x_max = 25.0_dp
    dx = (x_max - x_min) / real(nx, dp)
    mass = 2000.0_dp         ! 核质量 ~1.1 amu
    dt = 2.0_dp              ! 时间步长 (a.u.) ~ 0.05 fs

    do i = 1, nx
        x_grid(i) = x_min + real(i - 1, dp) * dx
        ! 绝热对角势与耦合势 Tully Model 1
        if (x_grid(i) > 0.0_dp) then
            v11(i) = 0.01_dp * (1.0_dp - exp(-1.6_dp * x_grid(i)))
        else
            v11(i) = -0.01_dp * (1.0_dp - exp(1.6_dp * x_grid(i)))
        end if
        v22(i) = -v11(i)
        v12(i) = 0.005_dp * exp(-x_grid(i)**2)
    end do

    ! 2. 准备初始波包 (位于态 1，x_0 = -10 a.u.，初动量 p_0 > 0)
    x_0 = -10.0_dp
    p_0 = 20.0_dp            ! 初始动量
    velocity = p_0 / mass
    sigma_x = 1.5_dp

    do i = 1, nx
        psi1(i) = (1.0_dp / (PI * sigma_x**2)**0.25_dp) * &
                  exp(-0.5_dp * ((x_grid(i) - x_0) / sigma_x)**2) * &
                  exp(EYE * p_0 * x_grid(i))
        psi2(i) = (0.0_dp, 0.0_dp)
    end do

    ! 3. Landau-Zener 理论跃迁几率估计
    ! 在 x=0 处: V_12 = 0.005, dV1/dx = 0.016, dV2/dx = -0.016 -> delta_slope = 0.032
    delta_slope = 0.032_dp
    p_lz = landau_zener_probability(0.005_dp, velocity, delta_slope)

    print '(A, F10.4)',      " Initial Wavepacket Position x_0: ", x_0
    print '(A, F10.4)',      " Initial Wavepacket Velocity v_0: ", velocity
    print '(A, F10.5)',      " Theoretical Landau-Zener P_LZ:  ", p_lz
    print '(A)', "----------------------------------------------------------"

    ! 4. 时间演化循环并记录两通道布居
    open(newunit=file_unit, file="nonadiabatic_dynamics.dat", status="replace", action="write")
    write(file_unit, '(A)') "# Time(fs)  Pop_State1  Pop_State2  Transition_Ratio"

    do step = 1, n_steps
        t = real(step, dp) * dt

        ! 2-通道分裂算符推进
        call propagate_split_operator_2channel(psi1, psi2, v11, v22, v12, dx, mass, dt)

        ! 统计通道几率
        if (mod(step, 10) == 0) then
            call calculate_channel_populations(psi1, psi2, dx, pop1, pop2, ratio)
            write(file_unit, '(4ES16.8)') t * AU2FS, pop1, pop2, ratio
        end if
    end do
    close(file_unit)

    call calculate_channel_populations(psi1, psi2, dx, pop1, pop2, ratio)
    print '(A, F10.4)', "Final State 1 Population: ", pop1
    print '(A, F10.4)', "Final State 2 Population: ", pop2
    print '(A, F10.4)', "Numerical Transition Ratio:", ratio
    print '(A)', "Dynamics output saved to: nonadiabatic_dynamics.dat"
    print '(A)', "=========================================================="

end program ex06_two_state_nonadiabatic
