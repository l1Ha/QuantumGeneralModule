!> \file mod_surface_hopping_fssh.f90
!> \brief Tully 最少开关表面跳跃 (FSSH) 与非绝热混合量子-经典分子动力学模块
!> \details 实现 Tully's Fewest Switches Surface Hopping (FSSH) 算法，
!>          涵盖绝热/非绝热基底转换、非绝热导数耦合 (NACV)、核动量沿 NACV 方向能量重标度、
!>          经典 Tully 三大基准测试模型 (SAC, DAC, ECR)、以及平均场 Ehrenfest 动力学。
!> \author LiHao
module mod_surface_hopping_fssh
    use mod_constants, only: dp, PI, EYE
    implicit none
    private

    public :: fssh_trajectory_t
    public :: tully_model_t
    public :: TULLY_SAC, TULLY_DAC, TULLY_ECR
    public :: init_tully_model
    public :: calc_adiabatic_surface_and_nacv
    public :: init_fssh_trajectory
    public :: propagate_fssh_step
    public :: run_fssh_ensemble
    public :: propagate_ehrenfest_step

    integer, parameter :: TULLY_SAC = 1  !< 单避免交叉 (Simple Avoided Crossing)
    integer, parameter :: TULLY_DAC = 2  !< 双避免交叉与斯托克斯干涉 (Dual Avoided Crossing)
    integer, parameter :: TULLY_ECR = 3  !< 扩展耦合反射 (Extended Coupling with Reflection)

    !> Tully 解析势能模型参数结构体
    type :: tully_model_t
        integer  :: model_type
        real(dp) :: a, b, c, d, e0
    end type tully_model_t

    !> 单条 FSSH 混合量子-经典动力学轨迹
    type :: fssh_trajectory_t
        real(dp) :: mass                  !< 核质量 (a.u.)
        real(dp) :: x                     !< 经典核坐标 (a.u.)
        real(dp) :: p                     !< 经典核动量 (a.u.)
        integer  :: active_state          !< 当前活性绝热电子态 (1 或 2)
        complex(dp) :: c(2)               !< 电子量子振幅 c_1, c_2
        real(dp) :: time                  !< 演化时间 (a.u.)
        integer  :: hop_count             !< 累计跳跃次数
        logical  :: frustrated_hop        !< 最近一步是否发生禁阻跳跃
    end type fssh_trajectory_t

contains

    !> \brief 初始化 Tully 测试模型
    !> \param[out] model 模型参数结构体
    !> \param[in] model_type TULLY_SAC, TULLY_DAC, 或 TULLY_ECR
    !> \param[out] stat 状态码
    subroutine init_tully_model(model, model_type, stat)
        type(tully_model_t), intent(out) :: model
        integer, intent(in)              :: model_type
        integer, optional, intent(out)   :: stat

        if (present(stat)) stat = 0
        model%model_type = model_type

        select case (model_type)
        case (TULLY_SAC)
            ! Tully 1990 J. Chem. Phys. 93, 1061 Model 1: SAC
            model%a  = 0.01_dp
            model%b  = 1.6_dp
            model%c  = 0.005_dp
            model%d  = 1.0_dp
            model%e0 = 0.0_dp
        case (TULLY_DAC)
            ! Model 2: DAC (双重交叉，呈现量子干涉效应)
            model%a  = 0.10_dp
            model%b  = 0.28_dp
            model%c  = 0.015_dp
            model%d  = 0.06_dp
            model%e0 = 0.05_dp
        case (TULLY_ECR)
            ! Model 3: ECR (扩展耦合反射)
            model%a  = 0.0006_dp
            model%b  = 0.10_dp
            model%c  = 0.09_dp
            model%d  = 0.0_dp
            model%e0 = 0.0_dp
        case default
            if (present(stat)) stat = -1
        end select
    end subroutine init_tully_model

    !> \brief 计算二能级透生表象哈密顿量与空间导数
    pure subroutine calc_diabatic_hamiltonian(model, x, v11, v22, v12, dv11, dv22, dv12)
        type(tully_model_t), intent(in) :: model
        real(dp), intent(in)            :: x
        real(dp), intent(out)           :: v11, v22, v12
        real(dp), intent(out)           :: dv11, dv22, dv12

        select case (model%model_type)
        case (TULLY_SAC)
            if (x >= 0.0_dp) then
                v11 = model%a * (1.0_dp - exp(-model%b * x))
                dv11 = model%a * model%b * exp(-model%b * x)
            else
                v11 = -model%a * (1.0_dp - exp(model%b * x))
                dv11 = model%a * model%b * exp(model%b * x)
            end if
            v22 = -v11
            dv22 = -dv11
            v12 = model%c * exp(-model%d * (x**2))
            dv12 = -2.0_dp * model%d * x * v12

        case (TULLY_DAC)
            v11 = 0.0_dp
            dv11 = 0.0_dp
            v22 = -model%a * exp(-model%b * (x**2)) + model%e0
            dv22 = 2.0_dp * model%a * model%b * x * exp(-model%b * (x**2))
            v12 = model%c * exp(-model%d * (x**2))
            dv12 = -2.0_dp * model%d * x * v12

        case default
            v11 = model%a
            dv11 = 0.0_dp
            v22 = -model%a
            dv22 = 0.0_dp
            if (x < 0.0_dp) then
                v12 = model%b * exp(model%c * x)
                dv12 = model%b * model%c * exp(model%c * x)
            else
                v12 = model%b * (2.0_dp - exp(-model%c * x))
                dv12 = model%b * model%c * exp(-model%c * x)
            end if
        end select
    end subroutine calc_diabatic_hamiltonian

    !> \brief 对角化求绝热势能、绝热梯度力与非绝热导数耦合标量 d_12(x)
    pure subroutine calc_adiabatic_surface_and_nacv(model, x, e1, e2, f1, f2, d12)
        type(tully_model_t), intent(in) :: model
        real(dp), intent(in)            :: x
        real(dp), intent(out)           :: e1, e2   !< 绝热势能面上、下本征值
        real(dp), intent(out)           :: f1, f2   !< 绝热力 F = -dE/dx
        real(dp), intent(out)           :: d12      !< 非绝热导数耦合 <1|d/dx|2>

        real(dp) :: v11, v22, v12, dv11, dv22, dv12
        real(dp) :: tr, det_diff, gap

        call calc_diabatic_hamiltonian(model, x, v11, v22, v12, dv11, dv22, dv12)

        tr = 0.5_dp * (v11 + v22)
        det_diff = 0.5_dp * (v11 - v22)
        gap = sqrt(det_diff**2 + v12**2)

        ! 绝热本征能面 E_1 <= E_2
        e1 = tr - gap
        e2 = tr + gap

        ! 绝热力 F_k = - dE_k / dx
        f1 = - (0.5_dp * (dv11 + dv22) - (det_diff * 0.5_dp * (dv11 - dv22) + v12 * dv12) / max(1.0e-14_dp, gap))
        f2 = - (0.5_dp * (dv11 + dv22) + (det_diff * 0.5_dp * (dv11 - dv22) + v12 * dv12) / max(1.0e-14_dp, gap))

        ! 混合角 theta: tan(2*theta) = 2*v12 / (v11 - v22)
        ! 非绝热导数耦合: d_12 = d(theta)/dx = (v12 * (dv11 - dv22)/2 - det_diff * dv12) / (gap^2)
        d12 = (v12 * 0.5_dp * (dv11 - dv22) - det_diff * dv12) / max(1.0e-14_dp, gap**2)
    end subroutine calc_adiabatic_surface_and_nacv

    !> \brief 初始化单条 FSSH 轨迹
    subroutine init_fssh_trajectory(traj, mass, x0, p0, initial_state)
        type(fssh_trajectory_t), intent(out) :: traj
        real(dp), intent(in)                 :: mass, x0, p0
        integer, intent(in)                  :: initial_state

        traj%mass = mass
        traj%x = x0
        traj%p = p0
        traj%active_state = initial_state
        traj%time = 0.0_dp
        traj%hop_count = 0
        traj%frustrated_hop = .false.

        traj%c = (0.0_dp, 0.0_dp)
        if (initial_state == 1) then
            traj%c(1) = (1.0_dp, 0.0_dp)
        else
            traj%c(2) = (1.0_dp, 0.0_dp)
        end if
    end subroutine init_fssh_trajectory

    !> \brief 单步推进 FSSH 轨迹 (Velocity-Verlet 核推进 + 电子态演化 + 跳跃概率评估)
    subroutine propagate_fssh_step(model, traj, dt, rand_val, stat)
        type(tully_model_t), intent(in)     :: model
        type(fssh_trajectory_t), intent(inout) :: traj
        real(dp), intent(in)                :: dt
        real(dp), intent(in)                :: rand_val !< 外部注入均匀随机数 [0, 1)
        integer, optional, intent(out)      :: stat

        real(dp) :: e1, e2, f1, f2, d12
        real(dp) :: f_active, v_nuc, v_nacv
        real(dp) :: p_old, delta_e, kinetic_nacv
        real(dp) :: g_hop, p11, p22, rho_re
        complex(dp) :: c1_new, c2_new
        integer  :: target_state

        if (present(stat)) stat = 0
        traj%frustrated_hop = .false.

        ! 1. 计算当前点绝热信息与活性力
        call calc_adiabatic_surface_and_nacv(model, traj%x, e1, e2, f1, f2, d12)
        if (traj%active_state == 1) then
            f_active = f1
        else
            f_active = f2
        end if

        ! 2. 核坐标 Velocity-Verlet 前半步: x(t+dt) = x(t) + v*dt + 0.5*(F/M)*dt^2
        v_nuc = traj%p / traj%mass
        traj%x = traj%x + v_nuc * dt + 0.5_dp * (f_active / traj%mass) * (dt**2)
        traj%p = traj%p + 0.5_dp * f_active * dt

        ! 3. 计算新位置处的绝热力与非绝热耦合
        call calc_adiabatic_surface_and_nacv(model, traj%x, e1, e2, f1, f2, d12)
        if (traj%active_state == 1) then
            f_active = f1
        else
            f_active = f2
        end if
        traj%p = traj%p + 0.5_dp * f_active * dt
        v_nuc = traj%p / traj%mass

        ! 4. 电子态演化: i*hbar dc_k/dt = E_k c_k - i*hbar (v * d_kj) c_j
        ! 标量 1D 下 v * d_12
        v_nacv = v_nuc * d12

        ! 二阶显式 Cayley / 中点传播法
        c1_new = traj%c(1) * exp(-EYE * e1 * dt) - dt * v_nacv * traj%c(2)
        c2_new = traj%c(2) * exp(-EYE * e2 * dt) + dt * v_nacv * traj%c(1)
        ! 保持幺正归一
        traj%c(1) = c1_new / sqrt(abs(c1_new)**2 + abs(c2_new)**2)
        traj%c(2) = c2_new / sqrt(abs(c1_new)**2 + abs(c2_new)**2)

        ! 5. Tully 最少开关跳跃几率计算
        ! g_{j -> k} = max(0, -2 * dt * Re(c_j * c_k^* * (v * d_jk)) / |c_j|^2)
        if (traj%active_state == 1) then
            target_state = 2
            p11 = abs(traj%c(1))**2
            rho_re = real(traj%c(1) * conjg(traj%c(2)), dp)
            ! d_12 = - d_21
            g_hop = max(0.0_dp, -2.0_dp * dt * (rho_re * (-v_nacv)) / max(1.0e-12_dp, p11))
            delta_e = e2 - e1
        else
            target_state = 1
            p22 = abs(traj%c(2))**2
            rho_re = real(traj%c(2) * conjg(traj%c(1)), dp)
            g_hop = max(0.0_dp, -2.0_dp * dt * (rho_re * (v_nacv)) / max(1.0e-12_dp, p22))
            delta_e = e1 - e2
        end if

        ! 6. 随机跳跃判定与能量守恒动量重标度
        if (rand_val < g_hop) then
            ! 检验动能是否足以补偿绝热势能差: P^2/(2M) >= delta_e
            kinetic_nacv = (traj%p**2) / (2.0_dp * traj%mass)
            if (kinetic_nacv >= delta_e) then
                ! 成功跳跃: 沿运动方向调整动量 P_new = sign(P) * sqrt(2*M*(E_kin - delta_e))
                p_old = traj%p
                traj%p = sign(1.0_dp, p_old) * sqrt(2.0_dp * traj%mass * (kinetic_nacv - delta_e))
                traj%active_state = target_state
                traj%hop_count = traj%hop_count + 1
            else
                ! 禁阻跳跃 (Frustrated Hop): 动能不足，维持原态
                traj%frustrated_hop = .true.
            end if
        end if

        traj%time = traj%time + dt
    end subroutine propagate_fssh_step

    !> \brief 运行 FSSH 轨迹系综并统计透射与反射分支比
    subroutine run_fssh_ensemble(model, initial_p, initial_state, n_trajs, dt, max_steps, &
                                 t1_prob, t2_prob, r1_prob, r2_prob, rng_seed)
        type(tully_model_t), intent(in) :: model
        real(dp), intent(in)            :: initial_p
        integer, intent(in)             :: initial_state
        integer, intent(in)             :: n_trajs
        real(dp), intent(in)            :: dt
        integer, intent(in)             :: max_steps
        real(dp), intent(out)           :: t1_prob, t2_prob  !< 态 1, 2 的正向透射率
        real(dp), intent(out)           :: r1_prob, r2_prob  !< 态 1, 2 的反向反射率
        integer, intent(in)             :: rng_seed

        type(fssh_trajectory_t) :: traj
        integer  :: itraj, step
        real(dp) :: mass, x0, rand_val
        integer  :: seed_curr
        integer  :: n_t1, n_t2, n_r1, n_r2

        mass = 2000.0_dp  ! Tully 基准粒子质量 M = 2000 a.u.
        x0   = -8.0_dp    ! 初始发射位置 (远离碰撞区)

        n_t1 = 0; n_t2 = 0
        n_r1 = 0; n_r2 = 0
        seed_curr = rng_seed

        do itraj = 1, n_trajs
            call init_fssh_trajectory(traj, mass, x0, initial_p, initial_state)

            do step = 1, max_steps
                ! 简易线性同余随机数产生
                seed_curr = mod(seed_curr * 1664525 + 1013904223, 2147483647)
                rand_val = real(seed_curr, dp) / 2147483647.0_dp

                call propagate_fssh_step(model, traj, dt, rand_val)

                ! 碰撞退出判据
                if (traj%x > 8.0_dp .or. traj%x < -10.0_dp) exit
            end do

            ! 终态统计
            if (traj%x > 0.0_dp) then
                if (traj%active_state == 1) n_t1 = n_t1 + 1
                if (traj%active_state == 2) n_t2 = n_t2 + 1
            else
                if (traj%active_state == 1) n_r1 = n_r1 + 1
                if (traj%active_state == 2) n_r2 = n_r2 + 1
            end if
        end do

        t1_prob = real(n_t1, dp) / real(n_trajs, dp)
        t2_prob = real(n_t2, dp) / real(n_trajs, dp)
        r1_prob = real(n_r1, dp) / real(n_trajs, dp)
        r2_prob = real(n_r2, dp) / real(n_trajs, dp)
    end subroutine run_fssh_ensemble

    !> \brief 平均场 Ehrenfest 动力学单步推进对比求解器
    subroutine propagate_ehrenfest_step(model, traj, dt)
        type(tully_model_t), intent(in)     :: model
        type(fssh_trajectory_t), intent(inout) :: traj
        real(dp), intent(in)                :: dt

        real(dp) :: e1, e2, f1, f2, d12
        real(dp) :: f_ehrenfest, v_nuc, v_nacv
        complex(dp) :: c1_new, c2_new

        call calc_adiabatic_surface_and_nacv(model, traj%x, e1, e2, f1, f2, d12)

        ! 平均场力: F_Ehrenfest = |c1|^2 * F1 + |c2|^2 * F2 + 2 * Re(c1^* c2) * (E1 - E2) * d12
        f_ehrenfest = (abs(traj%c(1))**2) * f1 + (abs(traj%c(2))**2) * f2 + &
                      2.0_dp * real(conjg(traj%c(1)) * traj%c(2), dp) * (e1 - e2) * d12

        ! 核坐标 Velocity Verlet 推进
        v_nuc = traj%p / traj%mass
        traj%x = traj%x + v_nuc * dt + 0.5_dp * (f_ehrenfest / traj%mass) * (dt**2)
        traj%p = traj%p + 0.5_dp * f_ehrenfest * dt

        call calc_adiabatic_surface_and_nacv(model, traj%x, e1, e2, f1, f2, d12)
        f_ehrenfest = (abs(traj%c(1))**2) * f1 + (abs(traj%c(2))**2) * f2 + &
                      2.0_dp * real(conjg(traj%c(1)) * traj%c(2), dp) * (e1 - e2) * d12
        traj%p = traj%p + 0.5_dp * f_ehrenfest * dt
        v_nuc = traj%p / traj%mass

        ! 电子态幺正演化
        v_nacv = v_nuc * d12
        c1_new = traj%c(1) * exp(-EYE * e1 * dt) - dt * v_nacv * traj%c(2)
        c2_new = traj%c(2) * exp(-EYE * e2 * dt) + dt * v_nacv * traj%c(1)
        traj%c(1) = c1_new / sqrt(abs(c1_new)**2 + abs(c2_new)**2)
        traj%c(2) = c2_new / sqrt(abs(c1_new)**2 + abs(c2_new)**2)
        traj%time = traj%time + dt
    end subroutine propagate_ehrenfest_step

end module mod_surface_hopping_fssh
