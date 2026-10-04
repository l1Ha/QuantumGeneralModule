!> \brief 准经典轨迹动力学模块 (Quasi-Classical Trajectory, QCT)
!> \details 基于书中第17章与第18.8节规范实现原子-双原子反应散射经典轨迹模拟：
!>          - H + H2 经典 LEPS 解析势能面与解析梯度
!>          - 经典 Hamilton 方程与二阶交错 Velocity-Verlet 辛积分器
!>          - 反应产物 Jacobi 坐标分析与内部能量求解
!>          - 不透明度函数与 Monte Carlo 反应截面统计
!> \author LiHao
module mod_qct_dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_constants, only: PI, TWOPI, EV2AU
    implicit none
    private

    public :: qct_config_t
    public :: qct_trajectory_t
    public :: qct_result_t
    public :: qct_leps_param_t
    public :: init_qct_leps_param
    public :: eval_leps_energy_gradient
    public :: qct_init_trajectory
    public :: qct_propagate_step
    public :: qct_analyze_final_state
    public :: run_qct_ensemble
    public :: qct_cross_section

    !> \brief QCT 计算配置派生类型 (对应书中 17.8 节)
    type :: qct_config_t
        real(dp) :: e_coll      = 0.05512_dp  !< 碰撞动能 (Hartree, 约 1.5 eV)
        real(dp) :: b_max       = 3.0_dp      !< 最大碰撞参数截断 (bohr)
        real(dp) :: r_start     = 9.0_dp      !< 入射初始原子-双原子质心距离 (bohr)
        real(dp) :: r_end       = 11.0_dp     !< 反应产物分离截断判据 (bohr)
        real(dp) :: dt          = 3.0_dp      !< Velocity-Verlet 辛步长 (a.u.)
        integer  :: n_traj      = 200         !< 系综轨迹数
        integer  :: max_steps   = 8000        !< 轨迹步数预算上限
        integer  :: v_initial   = 0           !< 初始振动量子数
        integer  :: j_initial   = 0           !< 初始转动量子数
        integer  :: seed        = 20260101    !< 伪随机数种子
        real(dp) :: mass        = 1837.15_dp  !< 体系单原子质量 (m_e, H 原子)
    end type qct_config_t

    !> \brief 单条轨迹状态类型 (对应书中 17.8 节)
    type :: qct_trajectory_t
        real(dp) :: q(3, 3)                   !< 3个原子的三维笛卡尔坐标 (x,y,z; atom 1..3)
        real(dp) :: p(3, 3)                   !< 3个原子的三维笛卡尔动量
        real(dp) :: time        = 0.0_dp      !< 传播时间 (a.u.)
        real(dp) :: b_impact    = 0.0_dp      !< 该轨迹碰撞参数 b (bohr)
        real(dp) :: e_total     = 0.0_dp      !< 体系总能量 (Hartree)
        real(dp) :: e_internal  = 0.0_dp      !< 产物双原子内部能量 (Hartree)
        integer  :: v_final     = -1          !< 终态振动量子数
        integer  :: j_final     = -1          !< 终态转动量子数
        logical  :: reactive    = .false.     !< 是否发生反应 (交换配对)
        logical  :: converged   = .false.     !< 是否达到渐近分离区正常结束
    end type qct_trajectory_t

    !> \brief 系综统计结果类型 (对应书中 17.8 节)
    type :: qct_result_t
        real(dp) :: cross_section   = 0.0_dp  !< 反应截面 (bohr^2)
        real(dp) :: stat_error      = 0.0_dp  !< 二项分布标准统计误差 (bohr^2)
        integer  :: n_traj          = 0       !< 实际总轨迹数
        integer  :: n_reactive      = 0       !< 反应轨迹数
    end type qct_result_t

    !> \brief H3 LEPS 势能面参数类型
    type :: qct_leps_param_t
        real(dp) :: d_e   = 0.1744158_dp      !< Morse 阱深 (Hartree, 约 4.746 eV)
        real(dp) :: beta  = 1.0276617_dp      !< Morse 范围参数 (bohr^-1)
        real(dp) :: r_e   = 1.4014199_dp      !< Morse 平衡核间距 (bohr)
    end type qct_leps_param_t

contains

    !> \brief 初始化标准 H3 LEPS 势能面参数
    pure subroutine init_qct_leps_param(par)
        type(qct_leps_param_t), intent(out) :: par
        par%d_e  = 4.746_dp * EV2AU
        par%beta = 1.942_dp * 0.529177_dp
        par%r_e  = 0.7416_dp / 0.529177_dp
    end subroutine init_qct_leps_param

    !> \brief 计算 LEPS 势能与 3 个原子的解析受力 F = -grad(V)
    subroutine eval_leps_energy_gradient(q, par, v_pot, force)
        real(dp), intent(in)               :: q(3, 3)
        type(qct_leps_param_t), intent(in) :: par
        real(dp), intent(out)              :: v_pot
        real(dp), intent(out)              :: force(3, 3)

        real(dp) :: diff(3, 3), r(3), u_vec(3, 3)
        real(dp) :: x(3), em(3), q_int(3), j_int(3), dq(3), dj(3)
        real(dp) :: delta, sq, dv(3), coef(3), g
        integer  :: k, a_idx, b_idx
        integer, parameter :: p_a(3) = (/1, 2, 3/)
        integer, parameter :: p_b(3) = (/2, 3, 1/)

        ! 计算三对核间距向量与欧氏距离: pair 1:(1,2), pair 2:(2,3), pair 3:(3,1)
        do k = 1, 3
            a_idx = p_a(k)
            b_idx = p_b(k)
            diff(:, k) = q(:, b_idx) - q(:, a_idx)
            r(k) = max(sqrt(sum(diff(:, k)**2)), 1.0e-12_dp)
            u_vec(:, k) = diff(:, k) / r(k)
        end do

        ! 计算对势 Coulomb 积分 Q_k 与 Exchange 积分 J_k 及其导数
        do k = 1, 3
            x(k) = r(k) - par%r_e
            em(k) = exp(-par%beta * x(k))
            q_int(k) = 0.5_dp * par%d_e * (1.5_dp * em(k)**2 - em(k))
            j_int(k) = 0.25_dp * par%d_e * (em(k)**2 - 6.0_dp * em(k))
            dq(k) = 0.5_dp * par%d_e * par%beta * (em(k) - 3.0_dp * em(k)**2)
            dj(k) = 0.25_dp * par%d_e * par%beta * (6.0_dp * em(k) - 2.0_dp * em(k)**2)
        end do

        ! London 交换根: Delta = 1/2 [ (J1-J2)^2 + (J2-J3)^2 + (J3-J1)^2 ]
        delta = 0.5_dp * ((j_int(1) - j_int(2))**2 + &
                          (j_int(2) - j_int(3))**2 + &
                          (j_int(3) - j_int(1))**2)
        sq = max(sqrt(max(delta, 1.0e-300_dp)), 1.0e-15_dp)
        v_pot = sum(q_int) - sq

        ! 解析链式法则求导 dV/dr_k
        do k = 1, 3
            coef(k) = (3.0_dp * j_int(k) - sum(j_int)) / (2.0_dp * sq)
            dv(k) = dq(k) - coef(k) * dj(k)
        end do

        ! 映射回笛卡尔力 F = -grad(V)
        force = 0.0_dp
        do k = 1, 3
            a_idx = p_a(k)
            b_idx = p_b(k)
            g = dv(k)
            force(:, a_idx) = force(:, a_idx) + g * u_vec(:, k)
            force(:, b_idx) = force(:, b_idx) - g * u_vec(:, k)
        end do
    end subroutine eval_leps_energy_gradient

    !> \brief 单步 Velocity-Verlet 辛积分器
    subroutine qct_propagate_step(q, p, dt, mass, par, force)
        real(dp), intent(inout)            :: q(3, 3), p(3, 3)
        real(dp), intent(in)               :: dt, mass
        type(qct_leps_param_t), intent(in) :: par
        real(dp), intent(inout), optional  :: force(3, 3)

        real(dp) :: f_curr(3, 3), f_next(3, 3), v_dummy

        if (present(force)) then
            f_curr = force
        else
            call eval_leps_energy_gradient(q, par, v_dummy, f_curr)
        end if

        p = p + 0.5_dp * dt * f_curr
        q = q + dt * (p / mass)
        call eval_leps_energy_gradient(q, par, v_dummy, f_next)
        p = p + 0.5_dp * dt * f_next

        if (present(force)) force = f_next
    end subroutine qct_propagate_step

    !> \brief 随机数发生器生成服从标准正态分布的 3D 向量并归一化
    subroutine random_unit_vector(vec, seed)
        real(dp), intent(out) :: vec(3)
        integer, intent(inout) :: seed
        real(dp) :: u1, u2, r_sq
        integer :: i

        do i = 1, 3
            seed = mod(seed * 1664525 + 1013904223, 2147483647)
            u1 = max(real(seed, dp) / 2147483647.0_dp, 1.0e-12_dp)
            seed = mod(seed * 1664525 + 1013904223, 2147483647)
            u2 = real(seed, dp) / 2147483647.0_dp
            vec(i) = sqrt(-2.0_dp * log(u1)) * cos(TWOPI * u2)
        end do
        r_sq = sum(vec**2)
        if (r_sq > 1.0e-15_dp) vec = vec / sqrt(r_sq)
    end subroutine random_unit_vector

    !> \brief 初始化单条反应散射准经典轨迹 (不变环面抽样)
    subroutine qct_init_trajectory(cfg, traj, seed, par)
        type(qct_config_t), intent(in)      :: cfg
        type(qct_trajectory_t), intent(out) :: traj
        integer, intent(inout)              :: seed
        type(qct_leps_param_t), intent(in)  :: par

        real(dp) :: xi1, b_val, rx, speed, mu_big_r, mu_small_r, r_axis(3)
        real(dp) :: v_A(3), v_B(3), v_C(3), r_vec(3), v_dummy, f_init(3, 3)
        real(dp) :: r0, p_r0

        mu_big_r = (2.0_dp / 3.0_dp) * cfg%mass
        mu_small_r = 0.5_dp * cfg%mass

        ! 碰撞参数在圆盘内均匀抽样: b = b_max * sqrt(xi)
        seed = mod(seed * 1664525 + 1013904223, 2147483647)
        xi1 = real(seed, dp) / 2147483647.0_dp
        b_val = cfg%b_max * sqrt(xi1)
        traj%b_impact = b_val

        ! 初始相对位移与动量 (A 从 x > 0 向内沿 -x 方向发射)
        rx = sqrt(max(cfg%r_start**2 - b_val**2, 0.0_dp))
        r_vec = (/ rx, b_val, 0.0_dp /)
        speed = sqrt(2.0_dp * cfg%e_coll / mu_big_r)

        ! 分子取向随机各向同性抽样
        call random_unit_vector(r_axis, seed)

        ! 初始双原子内核间距取平衡点
        r0 = par%r_e
        p_r0 = 0.0_dp

        ! 转换为 3 个原子的质心系坐标与动量
        v_A = (/ -speed, 0.0_dp, 0.0_dp /) * (2.0_dp / 3.0_dp)
        v_B = (/  speed, 0.0_dp, 0.0_dp /) * (1.0_dp / 3.0_dp) - 0.5_dp * (p_r0 / mu_small_r) * r_axis
        v_C = (/  speed, 0.0_dp, 0.0_dp /) * (1.0_dp / 3.0_dp) + 0.5_dp * (p_r0 / mu_small_r) * r_axis

        traj%q(:, 1) = r_vec
        traj%q(:, 2) = -0.5_dp * r0 * r_axis
        traj%q(:, 3) =  0.5_dp * r0 * r_axis

        traj%p(:, 1) = cfg%mass * v_A
        traj%p(:, 2) = cfg%mass * v_B
        traj%p(:, 3) = cfg%mass * v_C

        call eval_leps_energy_gradient(traj%q, par, v_dummy, f_init)
        traj%e_total = sum(traj%p**2) / (2.0_dp * cfg%mass) + v_dummy
        traj%time = 0.0_dp
        traj%reactive = .false.
        traj%converged = .false.
    end subroutine qct_init_trajectory

    !> \brief 判定产物配对并指认终态物理量
    subroutine qct_analyze_final_state(traj, cfg, par)
        type(qct_trajectory_t), intent(inout) :: traj
        type(qct_config_t), intent(in)        :: cfg
        type(qct_leps_param_t), intent(in)    :: par

        real(dp) :: r12, r23, r31, r_min
        real(dp) :: r_vec(3), p_rel(3), mu_r, em, v_morse
        integer  :: p_pair(2)

        r12 = sqrt(sum((traj%q(:, 2) - traj%q(:, 1))**2))
        r23 = sqrt(sum((traj%q(:, 3) - traj%q(:, 2))**2))
        r31 = sqrt(sum((traj%q(:, 1) - traj%q(:, 3))**2))

        ! 找出间距最小的原子对，判定为结合的产物双原子
        r_min = min(r12, r23, r31)
        if (abs(r_min - r23) <= 1.0e-12_dp) then
            p_pair = (/ 2, 3 /)
            traj%reactive = .false.  ! 原初配对未发生改变 (非反应弹性/非弹性散射)
        else if (abs(r_min - r12) <= 1.0e-12_dp) then
            p_pair = (/ 1, 2 /)
            traj%reactive = .true.   ! 反应产物 AB + C
        else
            p_pair = (/ 3, 1 /)
            traj%reactive = .true.   ! 反应产物 AC + B
        end if

        mu_r = 0.5_dp * cfg%mass
        r_vec = traj%q(:, p_pair(2)) - traj%q(:, p_pair(1))
        p_rel = 0.5_dp * (traj%p(:, p_pair(2)) - traj%p(:, p_pair(1)))

        em = exp(-par%beta * (r_min - par%r_e))
        v_morse = par%d_e * (em**2 - 2.0_dp * em)
        traj%e_internal = sum(p_rel**2) / (2.0_dp * mu_r) + v_morse
    end subroutine qct_analyze_final_state

    !> \brief 执行全套 QCT 系综计算并汇总反应截面与统计误差
    subroutine run_qct_ensemble(cfg, res, par)
        type(qct_config_t), intent(in)     :: cfg
        type(qct_result_t), intent(out)    :: res
        type(qct_leps_param_t), intent(in), optional :: par

        type(qct_leps_param_t) :: p_use
        type(qct_trajectory_t) :: traj
        integer  :: itraj, step, seed_curr, n_rxn
        real(dp) :: f_curr(3, 3), r_sep, com(3)

        if (present(par)) then
            p_use = par
        else
            call init_qct_leps_param(p_use)
        end if

        res%n_traj = cfg%n_traj
        n_rxn = 0
        seed_curr = cfg%seed

        do itraj = 1, cfg%n_traj
            call qct_init_trajectory(cfg, traj, seed_curr, p_use)
            call eval_leps_energy_gradient(traj%q, p_use, r_sep, f_curr)

            do step = 1, cfg%max_steps
                call qct_propagate_step(traj%q, traj%p, cfg%dt, cfg%mass, p_use, f_curr)
                traj%time = traj%time + cfg%dt

                ! 判定产物分离截断
                com = 0.5_dp * (traj%q(:, 2) + traj%q(:, 3))
                r_sep = sqrt(sum((traj%q(:, 1) - com)**2))
                if (r_sep > cfg%r_end) then
                    traj%converged = .true.
                    exit
                end if
            end do

            call qct_analyze_final_state(traj, cfg, p_use)
            if (traj%reactive) n_rxn = n_rxn + 1
        end do

        res%n_reactive = n_rxn
        call qct_cross_section(n_rxn, cfg%n_traj, cfg%b_max, res%cross_section, res%stat_error)
    end subroutine run_qct_ensemble

    !> \brief 计算 Monte Carlo 反应截面与标准统计误差
    pure subroutine qct_cross_section(n_rxn, n_tot, b_max, sigma, err)
        integer, intent(in)   :: n_rxn, n_tot
        real(dp), intent(in)  :: b_max
        real(dp), intent(out) :: sigma, err
        real(dp) :: p_est

        if (n_tot <= 0) then
            sigma = 0.0_dp
            err   = 0.0_dp
            return
        end if

        p_est = real(n_rxn, dp) / real(n_tot, dp)
        sigma = PI * (b_max**2) * p_est
        err   = PI * (b_max**2) * sqrt(max(p_est * (1.0_dp - p_est), 0.0_dp) / real(n_tot, dp))
    end subroutine qct_cross_section

end module mod_qct_dynamics
