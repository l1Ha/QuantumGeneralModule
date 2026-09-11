!> \brief 高次谐波（HHG）与强场阿秒光电子能谱分析模块
!> \details 提供 Ehrenfest 偶极加速度、含时偶极矩、Gabor 时频小波变换、
!>          HHG 谐波功率谱及强场近似（SFA / Lewenstein 模型）偶极发射计算。
!> \author LiHao
module mod_hhg_spectra
    use mod_constants, only: dp, PI, TWOPI, EYE
    use mod_linear_algebra, only: fft_1d
    implicit none
    private

    public :: calculate_dipole_length
    public :: calculate_dipole_acceleration
    public :: hhg_power_spectrum
    public :: gabor_transform_point
    public :: lewenstein_sfa_dipole

contains

    !> \brief 计算长度表象瞬时偶极矩: d(t) = <psi(t)| x |psi(t)>
    pure function calculate_dipole_length(x_grid, dx, psi) result(d_val)
        real(dp), intent(in) :: x_grid(:)
        real(dp), intent(in) :: dx
        complex(dp), intent(in) :: psi(:)
        real(dp) :: d_val
        integer :: i, n

        d_val = 0.0_dp
        n = size(psi)
        do i = 1, n
            d_val = d_val + (abs(psi(i))**2) * x_grid(i) * dx
        end do
    end function calculate_dipole_length

    !> \brief 根据 Ehrenfest 定理计算瞬时偶极加速度:
    !>        a(t) = - <psi(t)| dV/dx |psi(t)> - E(t) * <psi(t)|psi(t)>
    pure function calculate_dipole_acceleration(dv_dx, dx, psi, e_field) result(a_val)
        real(dp), intent(in) :: dv_dx(:)
        real(dp), intent(in) :: dx
        complex(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: e_field
        real(dp) :: a_val
        real(dp) :: force_pot, norm_sq
        integer :: i, n

        force_pot = 0.0_dp
        norm_sq = 0.0_dp
        n = size(psi)

        do i = 1, n
            norm_sq = norm_sq + (abs(psi(i))**2) * dx
            force_pot = force_pot + (abs(psi(i))**2) * dv_dx(i) * dx
        end do

        a_val = -force_pot - e_field * norm_sq
    end function calculate_dipole_acceleration

    !> \brief 计算高次谐波功率发射谱: S(omega) = |FFT[ a(t) * W(t) ]|^2
    !> \param[in] a_t 偶极加速度时域序列 (n_steps)
    !> \param[in] dt 时间步长 (a.u.)
    !> \param[out] omega_arr 谐波频率网格 (n_steps / 2)
    !> \param[out] spectrum 谐波功率谱值 (n_steps / 2)
    subroutine hhg_power_spectrum(a_t, dt, omega_arr, spectrum)
        real(dp), intent(in) :: a_t(:)
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: omega_arr(:)
        real(dp), intent(out) :: spectrum(:)

        integer :: n, n_out, i
        complex(dp), allocatable :: work(:)
        real(dp) :: window, d_omega

        n = size(a_t)
        n_out = min(size(omega_arr), size(spectrum), n / 2)

        allocate(work(n))

        ! 加窗 (Hanning Window) 抑制首尾突变产生的频域泄露
        do i = 1, n
            window = 0.5_dp * (1.0_dp - cos(TWOPI * real(i - 1, dp) / real(n - 1, dp)))
            work(i) = cmplx(a_t(i) * window, 0.0_dp, kind=dp)
        end do

        ! 执行快速傅里叶变换
        call fft_1d(work, -1)

        d_omega = TWOPI / (real(n, dp) * dt)
        do i = 1, n_out
            omega_arr(i) = real(i - 1, dp) * d_omega
            spectrum(i) = (abs(work(i) * dt)**2)
        end do

        deallocate(work)
    end subroutine hhg_power_spectrum

    !> \brief 计算特定时刻 t_0 与频率 omega 的 Gabor 时频变换幅值（发射时间轮廓）
    !> \details G(omega, t_0) = int a(t) * exp(-(t - t_0)^2 / (2 * sigma^2)) * exp(i * omega * t) dt
    pure function gabor_transform_point(t_arr, a_t, dt, t_0, omega, sigma_gabor) result(amp)
        real(dp), intent(in) :: t_arr(:)
        real(dp), intent(in) :: a_t(:)
        real(dp), intent(in) :: dt, t_0, omega, sigma_gabor
        real(dp) :: amp
        complex(dp) :: integral_val
        real(dp) :: window, phase
        integer :: i, n

        integral_val = (0.0_dp, 0.0_dp)
        n = size(t_arr)

        do i = 1, n
            window = exp(-0.5_dp * ((t_arr(i) - t_0) / sigma_gabor)**2)
            phase = omega * t_arr(i)
            integral_val = integral_val + a_t(i) * window * exp(EYE * phase) * dt
        end do

        amp = abs(integral_val)
    end function gabor_transform_point

    !> \brief Lewenstein 强场近似 (SFA) 鞍点一维氢型原子含时偶极矩
    !> \param[in] t 当前时刻 (a.u.)
    !> \param[in] t_arr 历史时刻网格 (1:it)
    !> \param[in] a_field 矢量势历史数据 A(t') (1:it)
    !> \param[in] e_field 电场历史数据 E(t') (1:it)
    !> \param[in] ip_au 电离势 (a.u.)
    !> \param[in] dt 时间积分微元
    pure function lewenstein_sfa_dipole(t, t_arr, a_field, e_field, ip_au, dt) result(dip)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: t_arr(:)
        real(dp), intent(in) :: a_field(:)
        real(dp), intent(in) :: e_field(:)
        real(dp), intent(in) :: ip_au, dt
        real(dp) :: dip

        integer :: it, it_p, n_steps
        real(dp) :: tp, tau, p_st, s_act, action_int
        real(dp) :: d_rec, d_ion, d_tau
        complex(dp) :: sum_dip, integrand

        sum_dip = (0.0_dp, 0.0_dp)
        n_steps = size(t_arr)
        if (n_steps <= 1) then
            dip = 0.0_dp
            return
        end if

        ! 回溯积分至多 1.5 个激光光学周期 (约 150 a.u.)
        do it_p = max(1, n_steps - 300), n_steps - 1
            tp = t_arr(it_p)
            tau = t - tp
            if (tau <= 1.0e-3_dp) cycle

            ! 鞍点正则动量: p_s(t, t') = 1/tau * int_t'^t A(t'') dt''
            ! 梯形求和近似
            p_st = 0.5_dp * (a_field(n_steps) + a_field(it_p))

            ! 拟经典作用量: S(t, t') = int_t'^t [ (p_s - A(t''))^2 / 2 + Ip ] dt''
            action_int = 0.5_dp * ((p_st - a_field(it_p))**2 + (p_st - a_field(n_steps))**2)
            s_act = (action_int + ip_au) * tau

            ! 偶极矩阵元 (氢原子 1s 态解析连续态偶极跃迁)
            ! d(v) = -i * (2^3.5 * (2*Ip)^1.25 / pi) * v / (v^2 + 2*Ip)^3
            d_ion = (p_st - a_field(it_p)) / (((p_st - a_field(it_p))**2 + 2.0_dp * ip_au)**3 + 1.0e-12_dp)
            d_rec = (p_st - a_field(n_steps)) / (((p_st - a_field(n_steps))**2 + 2.0_dp * ip_au)**3 + 1.0e-12_dp)

            ! 积分项: d^*(p - A(t)) * d(p - A(t')) * E(t') * exp(-i * S) / (tau + i*eps)^1.5
            d_tau = (tau**2 + 0.1_dp)**0.75_dp
            integrand = cmplx(d_rec * d_ion * e_field(it_p) / d_tau, 0.0_dp, kind=dp) * exp(-EYE * s_act) * dt
            sum_dip = sum_dip + integrand
        end do

        dip = -2.0_dp * aimag(sum_dip)
    end function lewenstein_sfa_dipole

end module mod_hhg_spectra
