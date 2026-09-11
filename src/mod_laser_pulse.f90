!> \brief 超快强激光脉冲时域合成与光场相互作用模块
!> \details 支持高斯、正弦平方、斜坡平顶、线性/高阶啁啾、双色场（omega + 2*omega）以及太赫兹（THz）脉冲串，
!>          提供动态 AC Stark 移动量计算与时间序列矢量化生成。
!> \author LiHao
module mod_laser_pulse
    use mod_constants, only: dp, PI, TWOPI, FS2AU, VM2AU
    implicit none
    private

    public :: PULSE_GAUSSIAN, PULSE_SIN2, PULSE_FLATTOP, PULSE_CHIRP, PULSE_TWOCOLOR, PULSE_THZ_TRAIN
    public :: pulse_config_t
    public :: pulse_envelope
    public :: pulse_electric_field
    public :: pulse_electric_field_2d
    public :: pulse_vector_potential
    public :: pulse_stark_shift
    public :: pulse_generate_timeseries

    ! 脉冲形状常量枚举
    integer, parameter :: PULSE_GAUSSIAN  = 1
    integer, parameter :: PULSE_SIN2      = 2
    integer, parameter :: PULSE_FLATTOP   = 3
    integer, parameter :: PULSE_CHIRP     = 4
    integer, parameter :: PULSE_TWOCOLOR  = 5
    integer, parameter :: PULSE_THZ_TRAIN = 6

    !> \brief 统一激光脉冲配置结构体
    type :: pulse_config_t
        integer  :: shape_type = PULSE_GAUSSIAN !< 脉冲包络类型
        real(dp) :: field_peak = 0.0_dp         !< 峰值场强 (a.u.)
        real(dp) :: freq_central = 0.0_dp       !< 载波中心频率 omega_0 (a.u.)
        real(dp) :: duration = 100.0_dp * FS2AU !< 脉冲持续时间/半高全宽 FWHM (a.u.)
        real(dp) :: t_center = 0.0_dp           !< 脉冲峰值中心时刻 t_0 (a.u.)
        real(dp) :: chirp_rate = 0.0_dp         !< 线性啁啾系数 beta (a.u.)
        real(dp) :: cep_phase = 0.0_dp          !< 载波包络相位 CEP (rad)
        real(dp) :: t_flattop = 0.0_dp          !< 平顶脉冲平台段持续时间 (a.u.)
        real(dp) :: ellipticity = 0.0_dp        !< 椭偏率 epsilon (-1到1，0为线偏振)
        ! 双色场参数
        real(dp) :: two_color_ratio = 0.0_dp    !< 二次谐波场强比 E_2w / E_w
        real(dp) :: two_color_phase = 0.0_dp    !< 双色相对相位 (rad)
        ! THz 脉冲串参数
        integer  :: train_count = 1             !< 脉冲序列中脉冲个数
        real(dp) :: train_delay = 0.0_dp        !< 脉冲串两两之间时间间隔 (a.u.)
    end type pulse_config_t

contains

    !> \brief 计算时刻 t 对应的脉冲无量纲瞬时包络 f(t) ∈ [0, 1]
    pure function pulse_envelope(t, cfg) result(env)
        real(dp), intent(in) :: t
        type(pulse_config_t), intent(in) :: cfg
        real(dp) :: env
        real(dp) :: dt, tau, t_ramp, k_dt
        integer :: k

        env = 0.0_dp
        dt = t - cfg%t_center
        tau = max(1.0e-12_dp, cfg%duration)

        select case(cfg%shape_type)
        case(PULSE_GAUSSIAN, PULSE_CHIRP, PULSE_TWOCOLOR)
            ! 高斯包络: exp(-2*ln(2) * (t/tau)^2)
            env = exp(-2.772588722239781_dp * (dt / tau)**2)

        case(PULSE_SIN2)
            ! sin^2 包络在 [-tau, tau] 之间
            if (abs(dt) <= tau) then
                env = sin(PI * (dt + tau) / (2.0_dp * tau))**2
            else
                env = 0.0_dp
            end if

        case(PULSE_FLATTOP)
            ! 梯形平顶脉冲: 上升沿 (tau) + 平顶 (t_flattop) + 下降沿 (tau)
            t_ramp = tau
            if (abs(dt) <= 0.5_dp * cfg%t_flattop) then
                env = 1.0_dp
            else if (abs(dt) <= 0.5_dp * cfg%t_flattop + t_ramp) then
                env = 0.5_dp * (1.0_dp + cos(PI * (abs(dt) - 0.5_dp * cfg%t_flattop) / t_ramp))
            else
                env = 0.0_dp
            end if

        case(PULSE_THZ_TRAIN)
            ! 多单周期脉冲串包络叠加
            do k = 0, cfg%train_count - 1
                k_dt = t - (cfg%t_center + real(k, dp) * cfg%train_delay)
                env = env + exp(-2.772588722239781_dp * (k_dt / tau)**2)
            end do
            env = min(1.0_dp, env)

        case default
            env = exp(-2.772588722239781_dp * (dt / tau)**2)
        end select
    end function pulse_envelope

    !> \brief 计算时刻 t 对应的瞬时激光电场标量值 E(t) (a.u.)
    pure function pulse_electric_field(t, cfg) result(efield)
        real(dp), intent(in) :: t
        type(pulse_config_t), intent(in) :: cfg
        real(dp) :: efield
        real(dp) :: dt, env, phase, e1, e2, k_dt
        integer :: k

        dt = t - cfg%t_center
        env = pulse_envelope(t, cfg)

        select case(cfg%shape_type)
        case(PULSE_CHIRP)
            ! 线性啁啾瞬时相位: omega_0 * t + 0.5 * beta * t^2 + phi
            phase = cfg%freq_central * dt + 0.5_dp * cfg%chirp_rate * dt**2 + cfg%cep_phase
            efield = cfg%field_peak * env * cos(phase)

        case(PULSE_TWOCOLOR)
            ! 双色场 E(t) = E_1 * cos(w*t + phi1) + E_2 * cos(2w*t + phi2)
            e1 = cfg%field_peak * env * cos(cfg%freq_central * dt + cfg%cep_phase)
            e2 = cfg%field_peak * cfg%two_color_ratio * env * &
                 cos(2.0_dp * cfg%freq_central * dt + cfg%two_color_phase)
            efield = e1 + e2

        case(PULSE_THZ_TRAIN)
            ! 太赫兹单周期串: E(t) = sum_k E_k * env_k * sin(w * t_k)
            efield = 0.0_dp
            do k = 0, cfg%train_count - 1
                k_dt = t - (cfg%t_center + real(k, dp) * cfg%train_delay)
                efield = efield + cfg%field_peak * exp(-2.772588722239781_dp * (k_dt / cfg%duration)**2) * &
                         sin(cfg%freq_central * k_dt + cfg%cep_phase)
            end do

        case default
            ! 标准载波振荡 E(t) = E_0 * f(t) * cos(omega*t + phi)
            efield = cfg%field_peak * env * cos(cfg%freq_central * dt + cfg%cep_phase)
        end select
    end function pulse_electric_field

    !> \brief 计算瞬时动态 AC Stark 能级位移量 (a.u.)
    !> \details Delta_S(t, theta) = -1/4 * [alpha_parallel * cos^2(theta) + alpha_perp * sin^2(theta)] * E^2(t)
    pure function pulse_stark_shift(t, cfg, alpha_parallel, alpha_perp, theta) result(shift)
        real(dp), intent(in) :: t
        type(pulse_config_t), intent(in) :: cfg
        real(dp), intent(in) :: alpha_parallel, alpha_perp, theta
        real(dp) :: shift
        real(dp) :: ef, costh, sinth, alpha_eff

        ef = pulse_electric_field(t, cfg)
        costh = cos(theta)
        sinth = sin(theta)
        alpha_eff = alpha_parallel * costh**2 + alpha_perp * sinth**2
        shift = -0.25_dp * alpha_eff * (ef**2)
    end function pulse_stark_shift

    !> \brief 矢量化生成完整时域电场与包络采样数组
    subroutine pulse_generate_timeseries(t_arr, cfg, e_arr, env_arr)
        real(dp), intent(in) :: t_arr(:)
        type(pulse_config_t), intent(in) :: cfg
        real(dp), intent(out) :: e_arr(:)
        real(dp), intent(out), optional :: env_arr(:)
        integer :: i, n

        n = size(t_arr)
        do i = 1, n
            e_arr(i) = pulse_electric_field(t_arr(i), cfg)
        end do

        if (present(env_arr)) then
            do i = 1, n
                env_arr(i) = pulse_envelope(t_arr(i), cfg)
            end do
        end if
    end subroutine pulse_generate_timeseries

    !> \brief 计算瞬时电场二维矢量分量 (Ex, Ey)
    pure subroutine pulse_electric_field_2d(t, cfg, ex, ey)
        real(dp), intent(in) :: t
        type(pulse_config_t), intent(in) :: cfg
        real(dp), intent(out) :: ex, ey
        real(dp) :: dt, env, phase, norm_ellip

        dt = t - cfg%t_center
        env = pulse_envelope(t, cfg)
        norm_ellip = 1.0_dp / sqrt(1.0_dp + cfg%ellipticity**2)

        phase = cfg%freq_central * dt + 0.5_dp * cfg%chirp_rate * dt**2 + cfg%cep_phase
        ex = cfg%field_peak * env * norm_ellip * cos(phase)
        ey = cfg%field_peak * env * norm_ellip * cfg%ellipticity * sin(phase)
    end subroutine pulse_electric_field_2d

    !> \brief 计算脉冲瞬时激光矢量势标量值 A(t) = - int_{-\infty}^t E(t') dt'
    !> \details 缓变包络近似下: A(t) = - E_0 / omega * f(t) * sin(omega*t + phi)
    pure function pulse_vector_potential(t, cfg) result(avec)
        real(dp), intent(in) :: t
        type(pulse_config_t), intent(in) :: cfg
        real(dp) :: avec
        real(dp) :: dt, env, phase, w

        w = max(1.0e-6_dp, cfg%freq_central)
        dt = t - cfg%t_center
        env = pulse_envelope(t, cfg)
        phase = w * dt + 0.5_dp * cfg%chirp_rate * dt**2 + cfg%cep_phase

        avec = -(cfg%field_peak / w) * env * sin(phase)
    end function pulse_vector_potential

end module mod_laser_pulse
