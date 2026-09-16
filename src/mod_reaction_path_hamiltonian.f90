!> \file mod_reaction_path_hamiltonian.f90
!> \brief 多原子反应路径哈密顿量 (RPH) 与变分过渡态理论 (CVT/VTST) 模块
!> \details 涵盖内禀反应坐标 (IRC) 路径参数化、反应路径哈密顿量、
!>          广义正则变分过渡态理论 (CVT)、熵瓶颈优化寻找、
!>          以及 Eckart 势垒半经典量子隧穿修正系数 kappa(T)。
!> \author LiHao
module mod_reaction_path_hamiltonian
    use mod_constants, only: dp, PI, TWOPI, KB
    implicit none
    private

    public :: rph_point_t
    public :: rph_path_t
    public :: init_rph_benchmark_reaction
    public :: calc_generalized_tst_rate
    public :: calc_cvt_rate_constant
    public :: calc_eckart_tunneling_factor

    real(dp), parameter :: PLANCK_H_SI = 6.62607015e-34_dp
    real(dp), parameter :: KB_SI       = 1.380649e-23_dp
    real(dp), parameter :: EV2J        = 1.602176634e-19_dp
    real(dp), parameter :: CM12J       = 1.98644586e-23_dp

    !> 沿内禀反应坐标 (IRC) 的单点结构与振动参数
    type :: rph_point_t
        real(dp) :: s_coord           !< 质量加权反应坐标弧长 s (amu^1/2 * Bohr)
        real(dp) :: v_pot_ev          !< 路径经典势能标量 (eV)
        integer  :: n_modes           !< 正交垂直振动自由度数目
        real(dp), allocatable :: frequencies_cm1(:) !< 各垂直正则振动模频率 (cm^-1)
        real(dp) :: curvature         !< 路径曲率标量 |kappa(s)| (amu^-1/2 * Bohr^-1)
    end type rph_point_t

    !> 全反应路径集合
    type :: rph_path_t
        integer :: n_points
        type(rph_point_t), allocatable :: points(:)
        real(dp) :: v_reactants_ev    !< 反应物渐近势能 (eV)
        real(dp), allocatable :: freq_reactants_cm1(:) !< 反应物振动频率 (cm^-1)
        real(dp) :: barrier_height_ev !< 鞍点经典能垒高度 (eV)
        real(dp) :: imag_freq_ts_cm1  !< 过渡态虚频大小 |omega_TS| (cm^-1)
    end type rph_path_t

contains

    !> \brief 初始化经典基准反应路径 (如 H + H2 -> H2 + H 或 CH4 + H 氢迁移模型)
    subroutine init_rph_benchmark_reaction(path, reaction_type, stat)
        type(rph_path_t), intent(out)  :: path
        integer, intent(in)            :: reaction_type
        integer, optional, intent(out) :: stat

        integer, parameter :: N_PTS = 31
        integer  :: i, n_m
        real(dp) :: s_val, v_val, v_barrier, s_width

        if (present(stat)) stat = 0
        n_m = 2  ! 2 个垂直广义振动模
        path%n_points = N_PTS
        allocate(path%points(N_PTS))
        allocate(path%freq_reactants_cm1(n_m))

        select case (reaction_type)
        case (2)
            ! CH4 + H -> CH3 + H2 氢提取模型势
            v_barrier = 0.65_dp
            s_width   = 0.60_dp
            path%imag_freq_ts_cm1 = 1800.0_dp
            path%freq_reactants_cm1(1) = 2900.0_dp
            path%freq_reactants_cm1(2) = 1350.0_dp
        case default
            ! 基准双曲正割共线双原子氢转移势垒 H + H2
            v_barrier = 0.42_dp     ! 能垒 ~0.42 eV (~9.7 kcal/mol)
            s_width   = 0.80_dp     ! 势垒宽度尺度 (amu^1/2 * Bohr)
            path%imag_freq_ts_cm1 = 1500.0_dp  ! TS 反应坐标虚频 ~1500i cm^-1
            path%freq_reactants_cm1(1) = 2200.0_dp
            path%freq_reactants_cm1(2) = 1100.0_dp
        end select

        path%v_reactants_ev = 0.0_dp
        path%barrier_height_ev = v_barrier

        do i = 1, N_PTS
            s_val = -1.5_dp + real(i - 1, dp) * (3.0_dp / real(N_PTS - 1, dp))
            path%points(i)%s_coord = s_val

            ! 势能剖面: V(s) = V_barrier / cosh^2(s / s_width)
            v_val = v_barrier / (cosh(s_val / s_width)**2)
            path%points(i)%v_pot_ev = v_val

            path%points(i)%n_modes = n_m
            allocate(path%points(i)%frequencies_cm1(n_m))

            ! 垂直振动频率随反应路径的连续软化/硬化演变
            path%points(i)%frequencies_cm1(1) = 2200.0_dp - 400.0_dp * exp(-(s_val / 0.7_dp)**2)
            path%points(i)%frequencies_cm1(2) = 1100.0_dp + 300.0_dp * exp(-(s_val / 0.7_dp)**2)

            ! 反应路径曲率在靠近过渡态拐弯处呈现峰值
            path%points(i)%curvature = 0.35_dp * abs(s_val) * exp(-(s_val / 0.6_dp)**2)
        end do
    end subroutine init_rph_benchmark_reaction

    !> \brief 计算给定温度 T 下位于反应坐标 s 处的广义过渡态速率 k^GTST(T, s)
    !> \details k^GTST(T, s) = (k_B * T / h) * (Q_vib^TS(T, s) / Q_vib^R(T)) * exp(-Delta V(s) / (k_B * T))
    pure subroutine calc_generalized_tst_rate(path, s_idx, temp_k, rate_gtst)
        type(rph_path_t), intent(in) :: path
        integer, intent(in)          :: s_idx
        real(dp), intent(in)         :: temp_k
        real(dp), intent(out)        :: rate_gtst

        integer  :: m
        real(dp) :: beta_si, factor_pre, q_vib_ts, q_vib_react
        real(dp) :: delta_v_j, exp_boltz, u_k

        if (temp_k <= 1.0e-3_dp .or. s_idx < 1 .or. s_idx > path%n_points) then
            rate_gtst = 0.0_dp
            return
        end if

        beta_si = 1.0_dp / (KB_SI * temp_k)
        factor_pre = (KB_SI * temp_k) / PLANCK_H_SI  ! k_B * T / h (s^-1)

        ! 垂直振动配分函数乘积: Q_vib = prod 1 / (1 - exp(-h*nu / k_B*T))
        q_vib_ts = 1.0_dp
        do m = 1, path%points(s_idx)%n_modes
            u_k = (path%points(s_idx)%frequencies_cm1(m) * CM12J) * beta_si
            u_k = max(1.0e-6_dp, min(50.0_dp, u_k))
            q_vib_ts = q_vib_ts * (1.0_dp / (1.0_dp - exp(-u_k)))
        end do

        q_vib_react = 1.0_dp
        do m = 1, size(path%freq_reactants_cm1)
            u_k = (path%freq_reactants_cm1(m) * CM12J) * beta_si
            u_k = max(1.0e-6_dp, min(50.0_dp, u_k))
            q_vib_react = q_vib_react * (1.0_dp / (1.0_dp - exp(-u_k)))
        end do

        delta_v_j = (path%points(s_idx)%v_pot_ev - path%v_reactants_ev) * EV2J
        exp_boltz = exp(-min(100.0_dp, delta_v_j * beta_si))

        rate_gtst = factor_pre * (q_vib_ts / max(1.0e-30_dp, q_vib_react)) * exp_boltz
    end subroutine calc_generalized_tst_rate

    !> \brief 正则变分过渡态理论 (CVT) 寻找最小化通量瓶颈位置 s*(T) 与速率常数
    !> \param[in] path 反应路径
    !> \param[in] temp_k 温度 (K)
    !> \param[out] s_opt 变分最优几何瓶颈位置 s*(T)
    !> \param[out] rate_cvt CVT 速率常数 (s^-1)
    !> \param[out] stat 状态码
    subroutine calc_cvt_rate_constant(path, temp_k, s_opt, rate_cvt, stat)
        type(rph_path_t), intent(in)   :: path
        real(dp), intent(in)           :: temp_k
        real(dp), intent(out)          :: s_opt
        real(dp), intent(out)          :: rate_cvt
        integer, optional, intent(out) :: stat

        integer  :: i, min_idx
        real(dp) :: rate_curr, min_rate

        if (present(stat)) stat = 0
        if (temp_k <= 0.0_dp) then
            if (present(stat)) stat = -1
            s_opt = 0.0_dp
            rate_cvt = 0.0_dp
            return
        end if

        min_rate = 1.0e35_dp
        min_idx = 1

        do i = 1, path%n_points
            call calc_generalized_tst_rate(path, i, temp_k, rate_curr)
            if (rate_curr < min_rate) then
                min_rate = rate_curr
                min_idx = i
            end if
        end do

        s_opt = path%points(min_idx)%s_coord
        rate_cvt = min_rate
    end subroutine calc_cvt_rate_constant

    !> \brief 计算对称/非对称 Eckart 势垒的半经典量子隧穿系数 kappa(T)
    !> \details 采用 Wigner / Skodje-Truhlar 截断解析修正:
    !>          在高温趋近于 1.0 + (1/24) * (h*|nu_TS| / k_B*T)^2
    !>          在低温通过 Eckart 解析透射率积分计算显著增强因子
    !> \param[in] barrier_height_ev 经典势垒高度 (eV)
    !> \param[in] imag_freq_cm1 过渡态虚频绝对值 |omega_TS| (cm^-1)
    !> \param[in] temp_k 温度 (K)
    !> \param[out] kappa_tunnel 量子隧穿增强因子 (>= 1.0)
    pure subroutine calc_eckart_tunneling_factor(barrier_height_ev, imag_freq_cm1, temp_k, kappa_tunnel)
        real(dp), intent(in)  :: barrier_height_ev
        real(dp), intent(in)  :: imag_freq_cm1
        real(dp), intent(in)  :: temp_k
        real(dp), intent(out) :: kappa_tunnel

        real(dp) :: u_ts, alpha_param, beta_param
        real(dp) :: v_barrier_j, beta_si

        if (temp_k <= 1.0_dp .or. imag_freq_cm1 <= 0.0_dp .or. barrier_height_ev <= 0.0_dp) then
            kappa_tunnel = 1.0_dp
            return
        end if

        v_barrier_j = barrier_height_ev * EV2J
        beta_si = 1.0_dp / (KB_SI * temp_k)

        ! 无量纲势垒虚频参数 u_TS = h * |nu_TS| / (k_B * T)
        u_ts = (imag_freq_cm1 * CM12J) * beta_si
        alpha_param = v_barrier_j * beta_si

        ! Wigner 近似 (弱隧穿极限)
        if (u_ts < 1.0_dp) then
            kappa_tunnel = 1.0_dp + (u_ts**2) / 24.0_dp
        else
            ! Skodje & Truhlar Eckart 强隧穿解析参数化 (J. Phys. Chem. 85, 624):
            beta_param = TWOPI * alpha_param / u_ts
            if (beta_param > u_ts) then
                kappa_tunnel = (0.5_dp * u_ts / sin(0.5_dp * min(3.10_dp, u_ts))) - &
                               (u_ts / (beta_param - u_ts)) * exp(-(beta_param - u_ts))
            else
                kappa_tunnel = (beta_param / (u_ts - beta_param)) * exp(u_ts - beta_param)
            end if
            kappa_tunnel = max(1.0_dp, min(1.0e6_dp, kappa_tunnel))
        end if
    end subroutine calc_eckart_tunneling_factor

end module mod_reaction_path_hamiltonian
