!> \brief 准经典轨迹动力学模块 (Quasi-Classical Trajectory, QCT)
!> \details 基于书中第17章与第18.8节规范实现原子-双原子反应散射经典轨迹模拟：
!>          - H + H2 经典 LEPS 解析势能面与解析梯度
!>          - 经典 Hamilton 方程与二阶交错 Velocity-Verlet 辛积分器
!>          - 反应产物 Jacobi 坐标分析与内部能量求解
!>          - 不透明度函数与 Monte Carlo 反应截面统计
!> \author LiHao
module mod_qct_dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use mod_constants, only: PI, TWOPI, EV2AU
    implicit none
    private

    public :: qct_config_t
    public :: qct_trajectory_t
    public :: qct_result_t
    public :: qct_leps_param_t
    public :: qct_sampler_t
    public :: init_qct_leps_param
    public :: eval_leps_energy_gradient
    public :: qct_sampler_build
    public :: qct_sample_trajectory
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

    !> \brief 不变环面初始采样器 (EBK 指认能量 + 一个振动周期的径向轨道表)
    !> \details 系综内全部轨迹共享同一 (v, j) 态, 轨道表只需构建一次; 各轨迹以
    !>          均匀时间 (等价均匀振动角) 在表上抽样相位, 避免"r 均匀抽样"的
    !>          经典偏置 (书中 17.2 节)。
    type :: qct_sampler_t
        real(dp) :: e_int = 0.0_dp            !< EBK 指认的内部振动能 (Hartree)
        real(dp), allocatable :: r_tab(:)     !< 周期内径向位置表 (bohr)
        real(dp), allocatable :: p_tab(:)     !< 对应径向动量表 (a.u.)
    end type qct_sampler_t

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
    !> \brief 确定性 64 位中间量 LCG 步进 (modulus 2^31-1, 无有符号溢出 UB)
    pure function lcg_next(seed) result(next_seed)
        integer, intent(in) :: seed
        integer :: next_seed
        integer(int64) :: s
        s = mod(int(seed, int64) * 1664525_int64 + 1013904223_int64, 2147483647_int64)
        next_seed = int(s)
    end function lcg_next

    subroutine random_unit_vector(vec, seed)
        real(dp), intent(out) :: vec(3)
        integer, intent(inout) :: seed
        real(dp) :: u1, u2, r_sq
        integer :: i

        do i = 1, 3
            seed = lcg_next(seed)
            u1 = max(real(seed, dp) / 2147483647.0_dp, 1.0e-12_dp)
            seed = lcg_next(seed)
            u2 = real(seed, dp) / 2147483647.0_dp
            vec(i) = sqrt(-2.0_dp * log(u1)) * cos(TWOPI * u2)
        end do
        r_sq = sum(vec**2)
        if (r_sq > 1.0e-15_dp) vec = vec / sqrt(r_sq)
    end subroutine random_unit_vector

    !> \brief 3 维叉积
    pure function cross3(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross3

    !> \brief Morse 势 (以离解极限为零参考)
    pure function qct_morse_potential(r, par) result(v)
        real(dp), intent(in) :: r
        type(qct_leps_param_t), intent(in) :: par
        real(dp) :: v, ex
        ex = exp(-par%beta * (r - par%r_e))
        v = par%d_e * (ex * ex - 2.0_dp * ex)
    end function qct_morse_potential

    !> \brief 有效径向势 V_eff(r) = V_Morse(r) + j2/(2 mu r^2), j2 = j(j+1)
    pure function qct_eff_potential(r, j2, par, mu_r) result(v)
        real(dp), intent(in) :: r, j2, mu_r
        type(qct_leps_param_t), intent(in) :: par
        real(dp) :: v
        v = qct_morse_potential(r, par) + j2 / (2.0_dp * mu_r * r * r)
    end function qct_eff_potential

    !> \brief 有效径向力 (解析式): F = -dV_eff/dr
    pure function qct_eff_force(r, j2, par, mu_r) result(f)
        real(dp), intent(in) :: r, j2, mu_r
        type(qct_leps_param_t), intent(in) :: par
        real(dp) :: f, ex
        ex = exp(-par%beta * (r - par%r_e))
        f = -2.0_dp * par%d_e * par%beta * (ex - ex * ex) + j2 / (mu_r * r**3)
    end function qct_eff_force

    !> \brief 求有效势极小点与给定能量 e_int 的内外径向转折点
    pure subroutine qct_turning_points(e_int, j2, par, mu_r, r_in, r_out, ok)
        real(dp), intent(in) :: e_int, j2, mu_r
        type(qct_leps_param_t), intent(in) :: par
        real(dp), intent(out) :: r_in, r_out
        logical, intent(out) :: ok
        real(dp) :: r_bot, r_lo, r_hi, r_mid, f_lo
        integer  :: k

        ok = .false.
        ! 极小点扫描 (粗扫 + 三分细化)
        r_bot = par%r_e
        do k = 1, 200
            r_lo = 0.70_dp * par%r_e + (1.30_dp * par%r_e - 0.70_dp * par%r_e) * (k - 1) / 199.0_dp
            if (qct_eff_potential(r_lo, j2, par, mu_r) < qct_eff_potential(r_bot, j2, par, mu_r)) r_bot = r_lo
        end do
        if (e_int <= qct_eff_potential(r_bot, j2, par, mu_r) + 1.0e-12_dp) return

        ! 内转折点: 在 [0.2 r_e, r_bot] 上二分 V_eff - e_int = 0
        r_lo = 0.2_dp * par%r_e
        r_hi = r_bot
        f_lo = qct_eff_potential(r_lo, j2, par, mu_r) - e_int
        if (f_lo < 0.0_dp) return
        do k = 1, 120
            r_mid = 0.5_dp * (r_lo + r_hi)
            if (qct_eff_potential(r_mid, j2, par, mu_r) - e_int > 0.0_dp) then
                r_lo = r_mid
            else
                r_hi = r_mid
            end if
            if (r_hi - r_lo < 1.0e-13_dp) exit
        end do
        r_in = 0.5_dp * (r_lo + r_hi)

        ! 外转折点: 从 1.5 r_e 起把右端点逐步外扩直到 V_eff > e_int
        r_hi = 1.5_dp * par%r_e
        do k = 1, 200
            if (qct_eff_potential(r_hi, j2, par, mu_r) > e_int) exit
            r_hi = 1.5_dp * r_hi
        end do
        if (qct_eff_potential(r_hi, j2, par, mu_r) <= e_int) return
        r_lo = r_bot
        do k = 1, 120
            r_mid = 0.5_dp * (r_lo + r_hi)
            if (qct_eff_potential(r_mid, j2, par, mu_r) - e_int < 0.0_dp) then
                r_lo = r_mid
            else
                r_hi = r_mid
            end if
            if (r_hi - r_lo < 1.0e-13_dp) exit
        end do
        r_out = 0.5_dp * (r_lo + r_hi)
        ok = .true.
    end subroutine qct_turning_points

    !> \brief 径向作用量 J(E) = 2 ∫ sqrt(2 mu (E - V_eff)) dr (EBK 量子化条件)
    !> \details 端点处 sqrt 型奇异以变量替换 r = r_tp ± s^2 消除, 替换后被积
    !>          函数在 s=0 邻域为 s 的光滑偶函数, 复合 Simpson (每半区 2000 区间)
    !>          精度优于 1e-12 (相对量级)。
    pure function qct_radial_action(e_int, j2, par, mu_r) result(j_act)
        real(dp), intent(in) :: e_int, j2, mu_r
        type(qct_leps_param_t), intent(in) :: par
        real(dp) :: j_act
        real(dp) :: r_in, r_out, r_mid, s_max, ds, s, gg
        real(dp) :: sum_even, sum_odd, integ
        integer  :: k, n_seg
        logical  :: ok

        call qct_turning_points(e_int, j2, par, mu_r, r_in, r_out, ok)
        if (.not. ok) then
            j_act = -1.0_dp
            return
        end if
        r_mid = 0.5_dp * (r_in + r_out)
        n_seg = 2000
        j_act = 0.0_dp

        ! 左半区: r = r_in + s^2, s: 0 -> sqrt(r_mid - r_in), dr = 2 s ds
        ! 端点 s=0 (转折点) 被积函数为零; 端点 s=s_max (r 中点) 非零必须计入
        s_max = sqrt(r_mid - r_in)
        ds = s_max / real(n_seg, dp)
        sum_even = 0.0_dp
        sum_odd = 0.0_dp
        do k = 1, n_seg - 1
            s = real(k, dp) * ds
            gg = qct_radius_integrand(r_in + s * s, e_int, j2, par, mu_r) * 2.0_dp * s
            if (mod(k, 2) == 0) then
                sum_even = sum_even + gg
            else
                sum_odd = sum_odd + gg
            end if
        end do
        gg = qct_radius_integrand(r_mid, e_int, j2, par, mu_r) * 2.0_dp * s_max
        integ = ds / 3.0_dp * (sum_odd * 4.0_dp + sum_even * 2.0_dp + gg)
        j_act = j_act + integ

        ! 右半区: r = r_out - s^2, s: 0 -> sqrt(r_out - r_mid)
        s_max = sqrt(r_out - r_mid)
        ds = s_max / real(n_seg, dp)
        sum_even = 0.0_dp
        sum_odd = 0.0_dp
        do k = 1, n_seg - 1
            s = real(k, dp) * ds
            gg = qct_radius_integrand(r_out - s * s, e_int, j2, par, mu_r) * 2.0_dp * s
            if (mod(k, 2) == 0) then
                sum_even = sum_even + gg
            else
                sum_odd = sum_odd + gg
            end if
        end do
        gg = qct_radius_integrand(r_mid, e_int, j2, par, mu_r) * 2.0_dp * s_max
        integ = ds / 3.0_dp * (sum_odd * 4.0_dp + sum_even * 2.0_dp + gg)
        j_act = j_act + integ

        j_act = 2.0_dp * j_act
    end function qct_radial_action

    !> \brief sqrt 型被积函数 sqrt(2 mu (E - V_eff(r))) (端点为 0, 由替换消除奇异)
    pure function qct_radius_integrand(r, e_int, j2, par, mu_r) result(g)
        real(dp), intent(in) :: r, e_int, j2, mu_r
        type(qct_leps_param_t), intent(in) :: par
        real(dp) :: g
        g = sqrt(max(2.0_dp * mu_r * (e_int - qct_eff_potential(r, j2, par, mu_r)), 0.0_dp))
    end function qct_radius_integrand

    !> \brief EBK/WKB 作用量量子化: 解 J(E) = 2 pi (v + 1/2) 得内部振动能
    !> \details 对纯 Morse 振子 (j=0) 该条件精确, 与 G(v) 一致到舍入; j>0 时
    !>          离心项使转动/振动耦合, 属书中 17.2 节说明的近似。
    subroutine qct_ebk_energy(v, j, par, mu_r, e_int, ok)
        integer, intent(in)  :: v, j
        type(qct_leps_param_t), intent(in) :: par
        real(dp), intent(in) :: mu_r
        real(dp), intent(out) :: e_int
        logical, intent(out) :: ok
        real(dp) :: j2, target, e_lo, e_hi, e_mid, r_in, r_out
        integer  :: k

        ok = .false.
        j2 = real(j * (j + 1), dp)
        target = TWOPI * (real(v, dp) + 0.5_dp)
        ! 下界: 有效势极小点 + 微小量
        e_lo = qct_eff_potential(par%r_e, j2, par, mu_r) + 1.0e-9_dp * par%d_e
        e_hi = -1.0e-7_dp * abs(par%d_e)
        if (qct_radial_action(e_lo, j2, par, mu_r) > target) return
        if (qct_radial_action(e_hi, j2, par, mu_r) < target) return
        do k = 1, 120
            e_mid = 0.5_dp * (e_lo + e_hi)
            if (qct_radial_action(e_mid, j2, par, mu_r) < target) then
                e_lo = e_mid
            else
                e_hi = e_mid
            end if
            if (e_hi - e_lo < 1.0e-14_dp * par%d_e) exit
        end do
        e_int = 0.5_dp * (e_lo + e_hi)
        call qct_turning_points(e_int, j2, par, mu_r, r_in, r_out, ok)
    end subroutine qct_ebk_energy

    !> \brief 构建一个振动周期的径向轨道表 (均匀时间 = 不变环面均匀振动角)
    subroutine qct_build_orbit_table(e_int, j2, par, mu_r, r_tab, p_tab, ok)
        real(dp), intent(in) :: e_int, j2, mu_r
        type(qct_leps_param_t), intent(in) :: par
        real(dp), allocatable, intent(out) :: r_tab(:), p_tab(:)
        logical, intent(out) :: ok
        integer, parameter :: n_tab = 4096
        real(dp) :: r_in, r_out, w0, t_guess, dt_o, t_period, dt_tab
        real(dp) :: r, p, a, r_old
        logical  :: left_turning
        integer  :: k, istep

        call qct_turning_points(e_int, j2, par, mu_r, r_in, r_out, ok)
        if (.not. ok) return

        w0 = par%beta * sqrt(2.0_dp * par%d_e / mu_r)
        t_guess = TWOPI / w0
        dt_o = t_guess / 20000.0_dp

        ! 周期搜索: 从外转折点以零径向动量出发, 积分至首次回归外转折点。
        ! 注意径向运动为约化质量动力学: dp/dt = F, dr/dt = p/mu_r。
        ! 回归判据采用"已越过内转折点、动量转正且逼近外转折点 1e-6 以内",
        ! 避免 r >= r_out 的严格比较在离散轨道上永不触发导致的搜索漂移。
        r = r_out * (1.0_dp - 1.0e-9_dp)
        p = 0.0_dp
        a = qct_eff_force(r, j2, par, mu_r)
        t_period = 0.0_dp
        left_turning = .false.
        do istep = 1, 1000000
            p = p + 0.5_dp * dt_o * a
            r_old = r
            r = r + dt_o * (p / mu_r)
            a = qct_eff_force(r, j2, par, mu_r)
            p = p + 0.5_dp * dt_o * a
            t_period = t_period + dt_o
            if (r_old < r_out - 1.0e-6_dp) left_turning = .true.
            if (left_turning .and. p > 0.0_dp .and. r > r_out - 1.0e-6_dp) exit
            if (t_period > 2.0_dp * t_guess) then
                ok = .false.
                return
            end if
        end do

        ! 表填充: 均匀时间步长的采样序列 (同一约化质量动力学)
        allocate(r_tab(n_tab), p_tab(n_tab))
        dt_tab = t_period / real(n_tab, dp)
        r = r_out * (1.0_dp - 1.0e-9_dp)
        p = 0.0_dp
        a = qct_eff_force(r, j2, par, mu_r)
        do k = 1, n_tab
            r_tab(k) = r
            p_tab(k) = p
            p = p + 0.5_dp * dt_tab * a
            r = r + dt_tab * (p / mu_r)
            a = qct_eff_force(r, j2, par, mu_r)
            p = p + 0.5_dp * dt_tab * a
        end do
        ok = .true.
    end subroutine qct_build_orbit_table

    !> \brief 初始化单条反应散射准经典轨迹 (不变环面抽样)
    !> \brief 构建不变环面采样器: EBK 指认 (v, j) 内部能量 + 一个振动周期轨道表
    subroutine qct_sampler_build(cfg, par, sampler, ok)
        type(qct_config_t), intent(in)      :: cfg
        type(qct_leps_param_t), intent(in)  :: par
        type(qct_sampler_t), intent(out)    :: sampler
        logical, intent(out)                :: ok
        real(dp) :: j2, mu_r

        mu_r = 0.5_dp * cfg%mass
        j2 = real(cfg%j_initial * (cfg%j_initial + 1), dp)
        call qct_ebk_energy(cfg%v_initial, cfg%j_initial, par, mu_r, sampler%e_int, ok)
        if (.not. ok) return
        call qct_build_orbit_table(sampler%e_int, j2, par, mu_r, &
                                   sampler%r_tab, sampler%p_tab, ok)
        if (.not. ok) then
            if (allocated(sampler%r_tab)) deallocate(sampler%r_tab)
            if (allocated(sampler%p_tab)) deallocate(sampler%p_tab)
        end if
    end subroutine qct_sampler_build

    !> \brief 从不变环面采样器抽取一条轨迹的初始条件 (均匀时间相位 + 各向同性取向)
    !> \details 相位在轨道表上均匀随机取点 (等价均匀振动角, 见书中 17.2 节);
    !>          转动方向 |j| = sqrt(j(j+1)) 随机各向同性, 轨道平面内相位随机;
    !>          碰撞参数在圆盘内均匀 (b = b_max sqrt(xi))。全部随机量来自局部
    !>          种子流, 给定 (seed, cfg, sampler) 完全确定, 与线程/调用顺序无关。
    subroutine qct_sample_trajectory(cfg, sampler, traj, seed, par)
        type(qct_config_t), intent(in)      :: cfg
        type(qct_sampler_t), intent(in)     :: sampler
        type(qct_trajectory_t), intent(out) :: traj
        integer, intent(inout)              :: seed
        type(qct_leps_param_t), intent(in)  :: par

        real(dp) :: xi1, b_val, rx, speed, mu_big_r, mu_small_r, r_axis(3)
        real(dp) :: v_A(3), v_B(3), v_C(3), r_vec(3), v_dummy, f_init(3, 3)
        real(dp) :: r0, p_r0, j_mag, phi, c_phi, s_phi
        real(dp) :: j_hat(3), j_vec(3), u0(3), w_vec(3), dot_w, nrm
        real(dp) :: p_r_vec(3), v_rel(3), r_dot(3)
        integer  :: idx, n_tab

        mu_big_r = (2.0_dp / 3.0_dp) * cfg%mass
        mu_small_r = 0.5_dp * cfg%mass

        ! 相位抽样: 轨道表上的均匀随机位置 (均匀时间)
        n_tab = size(sampler%r_tab)
        seed = lcg_next(seed)
        idx = 1 + int(real(seed, dp) / 2147483647.0_dp * real(n_tab - 1, dp) + 0.5_dp)
        idx = min(max(idx, 1), n_tab)
        r0 = sampler%r_tab(idx)
        p_r0 = sampler%p_tab(idx)

        ! 转动角动量方向 (各向同性) 与轨道平面内相位
        j_mag = sqrt(real(cfg%j_initial * (cfg%j_initial + 1), dp))
        call random_unit_vector(j_hat, seed)
        call random_unit_vector(w_vec, seed)
        dot_w = sum(w_vec * j_hat)
        u0 = w_vec - dot_w * j_hat
        nrm = sqrt(sum(u0 * u0))
        if (nrm < 1.0e-8_dp) then
            ! 万分之一概率级退化: 换用与 j_hat 最不平行的坐标轴构造正交基
            u0 = cross3(j_hat, (/ 1.0_dp, 0.0_dp, 0.0_dp /))
            nrm = sqrt(sum(u0 * u0))
            if (nrm < 1.0e-8_dp) u0 = cross3(j_hat, (/ 0.0_dp, 1.0_dp, 0.0_dp /))
            nrm = sqrt(sum(u0 * u0))
        end if
        u0 = u0 / nrm
        seed = lcg_next(seed)
        phi = TWOPI * real(seed, dp) / 2147483647.0_dp
        c_phi = cos(phi)
        s_phi = sin(phi)
        r_axis = u0 * c_phi + cross3(j_hat, u0) * s_phi
        r_axis = r_axis / sqrt(sum(r_axis * r_axis))
        j_vec = j_mag * j_hat

        ! 碰撞参数在圆盘内均匀抽样: b = b_max * sqrt(xi)
        seed = lcg_next(seed)
        xi1 = real(seed, dp) / 2147483647.0_dp
        b_val = cfg%b_max * sqrt(xi1)
        traj%b_impact = b_val

        ! 初始相对位移与动量 (A 从 x > 0 向内沿 -x 方向发射)
        rx = sqrt(max(cfg%r_start**2 - b_val**2, 0.0_dp))
        r_vec = (/ rx, b_val, 0.0_dp /)
        speed = sqrt(2.0_dp * cfg%e_coll / mu_big_r)
        r_dot = (/ -speed, 0.0_dp, 0.0_dp /)

        ! 径向动量矢量 = 振动项 + 转动项 j x r_axis / r0
        p_r_vec = p_r0 * r_axis + cross3(j_vec, r_axis) / r0
        v_rel = p_r_vec / mu_small_r

        ! 转换为 3 个原子的质心系坐标与动量 (等质量, 总动量为零)
        v_A = (2.0_dp / 3.0_dp) * r_dot
        v_B = -(1.0_dp / 3.0_dp) * r_dot - 0.5_dp * v_rel
        v_C = -(1.0_dp / 3.0_dp) * r_dot + 0.5_dp * v_rel

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
    end subroutine qct_sample_trajectory

    !> \brief 初始化单条反应散射准经典轨迹 (EBK 不变环面抽样)
    !> \details 便捷入口: 现场构建采样器后抽样一条轨迹。系综计算应改用
    !>          qct_sampler_build + qct_sample_trajectory 以避免重复建表。
    !>          若 EBK 指认失败 (非法 v/j), 退化为势阱底静止初始化以保证不中断。
    subroutine qct_init_trajectory(cfg, traj, seed, par)
        type(qct_config_t), intent(in)      :: cfg
        type(qct_trajectory_t), intent(out) :: traj
        integer, intent(inout)              :: seed
        type(qct_leps_param_t), intent(in)  :: par

        type(qct_sampler_t) :: sampler
        logical :: ok
        real(dp) :: v_dummy, f_init(3, 3), r_axis(3), speed_fb, mu_big_fb

        call qct_sampler_build(cfg, par, sampler, ok)
        if (ok) then
            call qct_sample_trajectory(cfg, sampler, traj, seed, par)
        else
            ! 回退: 势阱底静止 (与旧行为一致, 仅用于非法量子数输入)
            call random_unit_vector(r_axis, seed)
            mu_big_fb = (2.0_dp / 3.0_dp) * cfg%mass
            speed_fb = sqrt(2.0_dp * cfg%e_coll / mu_big_fb)
            traj%b_impact = 0.0_dp
            traj%q(:, 1) = (/ cfg%r_start, 0.0_dp, 0.0_dp /)
            traj%q(:, 2) = -0.5_dp * par%r_e * r_axis
            traj%q(:, 3) =  0.5_dp * par%r_e * r_axis
            traj%p(:, 1) = cfg%mass * (/ -speed_fb, 0.0_dp, 0.0_dp /) * (2.0_dp / 3.0_dp)
            traj%p(:, 2) = cfg%mass * (/ speed_fb, 0.0_dp, 0.0_dp /) * (1.0_dp / 3.0_dp)
            traj%p(:, 3) = traj%p(:, 2)
            call eval_leps_energy_gradient(traj%q, par, v_dummy, f_init)
            traj%e_total = sum(traj%p**2) / (2.0_dp * cfg%mass) + v_dummy
            traj%time = 0.0_dp
            traj%reactive = .false.
            traj%converged = .false.
        end if
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
    !> \details 轨迹间完全独立，是天然可并行的 Monte-Carlo 系综。为同时保证
    !>          (a) 任意线程数下的逐位可复现性与 (b) 与串行基线完全一致的随机
    !>          抽样序列，先串行预扫描把每条轨迹的 LCG 种子存入数组（每轨迹仅
    !>          4 字节），再以 OpenMP 并行循环各自独立初始化并传播。无 OpenMP
    !>          时哨兵指令被忽略，退化为原串行路径且结果不变。
    subroutine run_qct_ensemble(cfg, res, par)
        type(qct_config_t), intent(in)     :: cfg
        type(qct_result_t), intent(out)    :: res
        type(qct_leps_param_t), intent(in), optional :: par

        type(qct_leps_param_t) :: p_use
        type(qct_trajectory_t) :: traj
        type(qct_sampler_t) :: sampler
        integer, allocatable :: seeds(:)
        integer  :: itraj, step, seed_curr, seed_loc, n_rxn
        logical  :: ok_sampler
        real(dp) :: f_curr(3, 3), r_sep, com(3)

        if (present(par)) then
            p_use = par
        else
            call init_qct_leps_param(p_use)
        end if

        res%n_traj = cfg%n_traj
        n_rxn = 0

        ! 不变环面采样器 (系综内共享): EBK 能量与轨道表只构建一次
        call qct_sampler_build(cfg, p_use, sampler, ok_sampler)

        ! 每条轨迹独立的确定性种子流: seed_i = LCG^i(seed0), O(n) 预扫描,
        ! 与线程数及执行顺序无关, 保证任意并行配置下逐位可复现
        allocate(seeds(cfg%n_traj))
        seed_curr = cfg%seed
        do itraj = 1, cfg%n_traj
            seed_curr = lcg_next(seed_curr)
            seeds(itraj) = seed_curr
        end do

        !$omp parallel do default(shared) &
        !$omp& private(itraj, traj, step, seed_loc, f_curr, r_sep, com) &
        !$omp& reduction(+:n_rxn) schedule(dynamic, 64)
        do itraj = 1, cfg%n_traj
            seed_loc = seeds(itraj)
            if (ok_sampler) then
                call qct_sample_trajectory(cfg, sampler, traj, seed_loc, p_use)
            else
                call qct_init_trajectory(cfg, traj, seed_loc, p_use)
            end if
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
        !$omp end parallel do

        deallocate(seeds)
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
