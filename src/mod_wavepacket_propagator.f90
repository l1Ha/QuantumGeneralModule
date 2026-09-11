!> \brief 量子波包时间推进与常微分方程积分器核心模块
!> \details 包含二阶对称分裂算符（Split-Operator FFT）、四阶 Runge-Kutta（RK4）、
!>          Adams-Bashforth-Moulton (ABM) 4阶预估-校正法，以及光学 Bloch 方程数值求解器。
!> \author LiHao
module mod_wavepacket_propagator
    use mod_constants, only: dp, PI, TWOPI, EYE
    use mod_linear_algebra, only: fft_1d, fft_2d
    implicit none
    private

    public :: propagate_split_operator_1d
    public :: propagate_split_operator_2d
    public :: rk4_step
    public :: solve_bloch_two_level
    public :: abm4_step

contains

    !> \brief 一维二阶对称分裂算符步进 (Split-Operator FFT)
    !> \details psi(t+dt) = exp(-i*V*dt/2) * IFFT[ exp(-i*T(p)*dt) * FFT[ exp(-i*V*dt/2) * psi(t) ] ]
    !> \param[inout] psi 复数波函数向量 (n_x)
    !> \param[in] v_pot 局部坐标势能向量 V(x) (n_x)
    !> \param[in] dx 空间网格间距
    !> \param[in] mass 粒子约化质量 (a.u.)
    !> \param[in] dt 时间步长 (a.u.)
    subroutine propagate_split_operator_1d(psi, v_pot, dx, mass, dt)
        complex(dp), intent(inout) :: psi(:)
        real(dp), intent(in) :: v_pot(:)
        real(dp), intent(in) :: dx, mass, dt

        integer :: i, n
        real(dp) :: dp_k, p_val
        real(dp), allocatable :: p_kin(:)
        complex(dp), allocatable :: exp_v_half(:), exp_t(:)

        n = size(psi)
        allocate(p_kin(n))
        allocate(exp_v_half(n))
        allocate(exp_t(n))

        ! 1. 坐标表象动势能半步旋转算子: exp(-i * V(x) * dt / 2)
        do i = 1, n
            exp_v_half(i) = exp(-EYE * v_pot(i) * (0.5_dp * dt))
        end do

        ! 2. 动量空间网格与动能相位因子: p_k in [-pi/dx, pi/dx], exp(-i * p^2 / 2m * dt)
        dp_k = TWOPI / (real(n, dp) * dx)
        do i = 1, n
            if (i <= n / 2) then
                p_val = real(i - 1, dp) * dp_k
            else
                p_val = real(i - 1 - n, dp) * dp_k
            end if
            p_kin(i) = (p_val**2) / (2.0_dp * mass)
            exp_t(i) = exp(-EYE * p_kin(i) * dt)
        end do

        ! 第一步: 势能半步相乘
        psi = psi * exp_v_half

        ! 第二步: 变换到动量空间
        call fft_1d(psi, -1)

        ! 第三步: 动能整步相乘
        psi = psi * exp_t

        ! 第四步: 逆变换回坐标空间
        call fft_1d(psi, 1)

        ! 第五步: 势能后半步相乘
        psi = psi * exp_v_half

        deallocate(p_kin, exp_v_half, exp_t)
    end subroutine propagate_split_operator_1d

    !> \brief 通用标准 4 阶 Runge-Kutta (RK4) 推进单步
    !> \param[inout] y 状态向量 (dim)
    !> \param[in] t 当前时刻
    !> \param[in] dt 时间步长
    !> \param[in] f_deriv 导数函数接口 dy/dt = f(t, y)
    subroutine rk4_step(y, t, dt, f_deriv)
        real(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: t, dt
        interface
            pure function f_deriv(t_in, y_in) result(dydt)
                use mod_constants, only: dp
                real(dp), intent(in) :: t_in
                real(dp), intent(in) :: y_in(:)
                real(dp) :: dydt(size(y_in))
            end function f_deriv
        end interface

        integer :: n
        real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), y_tmp(:)

        n = size(y)
        allocate(k1(n), k2(n), k3(n), k4(n), y_tmp(n))

        k1 = f_deriv(t, y)
        y_tmp = y + 0.5_dp * dt * k1
        k2 = f_deriv(t + 0.5_dp * dt, y_tmp)
        y_tmp = y + 0.5_dp * dt * k2
        k3 = f_deriv(t + 0.5_dp * dt, y_tmp)
        y_tmp = y + dt * k3
        k4 = f_deriv(t + dt, y_tmp)

        y = y + (dt / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)

        deallocate(k1, k2, k3, k4, y_tmp)
    end subroutine rk4_step

    !> \brief 光学 Bloch 方程数值单步推进（求解二能级系统相干与布居反转）
    !> \details dR/dt = Omega x R，其中 R = (u, v, w)^T, Omega = (Omega_R, 0, Delta)^T
    !> \param[inout] r_bloch Bloch 矢量 [u, v, w]
    !> \param[in] rabi_freq 瞬时拉比频率 Omega_R(t)
    !> \param[in] detuning 瞬时有效失谐 Delta(t)
    !> \param[in] dt 时间步长
    subroutine solve_bloch_two_level(r_bloch, rabi_freq, detuning, dt)
        real(dp), intent(inout) :: r_bloch(3)
        real(dp), intent(in) :: rabi_freq, detuning, dt

        real(dp) :: k1(3), k2(3), k3(3), k4(3), r_tmp(3)

        ! k1
        k1 = bloch_torque(r_bloch, rabi_freq, detuning)
        ! k2
        r_tmp = r_bloch + 0.5_dp * dt * k1
        k2 = bloch_torque(r_tmp, rabi_freq, detuning)
        ! k3
        r_tmp = r_bloch + 0.5_dp * dt * k2
        k3 = bloch_torque(r_tmp, rabi_freq, detuning)
        ! k4
        r_tmp = r_bloch + dt * k3
        k4 = bloch_torque(r_tmp, rabi_freq, detuning)

        r_bloch = r_bloch + (dt / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end subroutine solve_bloch_two_level

    !> \brief Bloch 扭矩导数计算: dR/dt = Omega x R
    pure function bloch_torque(r, omega_r, delta) result(drdt)
        real(dp), intent(in) :: r(3)
        real(dp), intent(in) :: omega_r, delta
        real(dp) :: drdt(3)

        drdt(1) = -delta * r(2)
        drdt(2) = delta * r(1) - omega_r * r(3)
        drdt(3) = omega_r * r(2)
    end function bloch_torque

    !> \brief Adams-Bashforth-Moulton (ABM) 4阶预估-校正单步推进
    !> \param[inout] y 当前步解 y_n
    !> \param[in] fn 当前导数 f_n
    !> \param[in] fn_m1 前一步导数 f_{n-1}
    !> \param[in] fn_m2 前两步导数 f_{n-2}
    !> \param[in] fn_m3 前三步导数 f_{n-3}
    !> \param[in] dt 时间步长
    !> \param[in] f_eval 导数评估函数接口
    !> \param[in] t_next 下一步时刻 t_{n+1}
    subroutine abm4_step(y, fn, fn_m1, fn_m2, fn_m3, dt, t_next, f_eval)
        real(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: fn(:), fn_m1(:), fn_m2(:), fn_m3(:)
        real(dp), intent(in) :: dt, t_next
        interface
            pure function f_eval(t_in, y_in) result(dy)
                use mod_constants, only: dp
                real(dp), intent(in) :: t_in
                real(dp), intent(in) :: y_in(:)
                real(dp) :: dy(size(y_in))
            end function f_eval
        end interface

        integer :: n
        real(dp), allocatable :: y_pred(:), fn_p1(:)

        n = size(y)
        allocate(y_pred(n), fn_p1(n))

        ! Adams-Bashforth 4 阶预估: P = y_n + dt/24 * (55*fn - 59*fn-1 + 37*fn-2 - 9*fn-3)
        y_pred = y + (dt / 24.0_dp) * (55.0_dp * fn - 59.0_dp * fn_m1 + 37.0_dp * fn_m2 - 9.0_dp * fn_m3)

        ! 在预估点求导
        fn_p1 = f_eval(t_next, y_pred)

        ! Adams-Moulton 4 阶校正: y_{n+1} = y_n + dt/24 * (9*fn+1 + 19*fn - 5*fn-1 + fn-2)
        y = y + (dt / 24.0_dp) * (9.0_dp * fn_p1 + 19.0_dp * fn - 5.0_dp * fn_m1 + fn_m2)

        deallocate(y_pred, fn_p1)
    end subroutine abm4_step

    !> \brief 二维二阶对称分裂算符步进 (Split-Operator 2D FFT)
    !> \param[inout] psi 复数波函数二维数组 (nx, ny)
    !> \param[in] v_pot 局部坐标势能矩阵 V(x, y) (nx, ny)
    !> \param[in] dx, dy 坐标步长
    !> \param[in] mass_x, mass_y 对应自由度质量
    !> \param[in] dt 时间步长 (a.u.)
    subroutine propagate_split_operator_2d(psi, v_pot, dx, dy, mass_x, mass_y, dt)
        complex(dp), intent(inout) :: psi(:, :)
        real(dp), intent(in) :: v_pot(:, :)
        real(dp), intent(in) :: dx, dy, mass_x, mass_y, dt

        integer :: nx, ny, ix, iy
        real(dp) :: dpx, dpy, px, py, e_kin
        complex(dp), allocatable :: exp_v_half(:, :), exp_t(:, :)

        nx = size(psi, 1)
        ny = size(psi, 2)
        allocate(exp_v_half(nx, ny), exp_t(nx, ny))

        ! 坐标空间势能半步
        exp_v_half = exp(-EYE * v_pot * (0.5_dp * dt))

        ! 动量空间相位矩阵
        dpx = TWOPI / (real(nx, dp) * dx)
        dpy = TWOPI / (real(ny, dp) * dy)

        do iy = 1, ny
            if (iy <= ny / 2) then
                py = real(iy - 1, dp) * dpy
            else
                py = real(iy - 1 - ny, dp) * dpy
            end if
            do ix = 1, nx
                if (ix <= nx / 2) then
                    px = real(ix - 1, dp) * dpx
                else
                    px = real(ix - 1 - nx, dp) * dpx
                end if
                e_kin = (px**2) / (2.0_dp * mass_x) + (py**2) / (2.0_dp * mass_y)
                exp_t(ix, iy) = exp(-EYE * e_kin * dt)
            end do
        end do

        psi = psi * exp_v_half
        call fft_2d(psi, -1)
        psi = psi * exp_t
        call fft_2d(psi, 1)
        psi = psi * exp_v_half

        deallocate(exp_v_half, exp_t)
    end subroutine propagate_split_operator_2d

end module mod_wavepacket_propagator
