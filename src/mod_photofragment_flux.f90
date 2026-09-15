!> \brief 分子光解离动力学、碎片动能释放谱（KER）与自相关吸收谱模块
!> \details 提供波包自相关函数计算、光吸收截面谱（Heller 傅里叶积分）、
!>          渐近解离面能量分辨概率流通量法（Balint-Kurti / Atabek 通量公式）
!>          以及多通道解离分支比统计分析。
!> \author LiHao
module mod_photofragment_flux
    use mod_constants, only: dp, PI, TWOPI, EYE, HBAR, FS2AU, AU2EV
    implicit none
    private

    public :: calculate_autocorrelation
    public :: calculate_absorption_spectrum
    public :: energy_resolved_flux_amplitude
    public :: calculate_ker_spectrum
    public :: calculate_branching_ratios

contains

    !> \brief 计算初态与含时波包内积自相关函数 C(t) = <psi(0) | psi(t)>
    pure function calculate_autocorrelation(psi_0, psi_t, dx) result(c_val)
        complex(dp), intent(in) :: psi_0(:), psi_t(:)
        real(dp), intent(in) :: dx
        complex(dp) :: c_val

        c_val = sum(conjg(psi_0(:)) * psi_t(:)) * dx
    end function calculate_autocorrelation

    !> \brief 基于自相关函数傅里叶变换计算连续吸收谱 sigma(omega)
    !> \details sigma(omega) = C_0 * omega * Re \int_0^T C(t) * W(t) * exp(i*(E_0 + omega)*t) dt
    subroutine calculate_absorption_spectrum(t_arr, c_arr, e0_au, omega_arr, spectrum, dt)
        real(dp), intent(in) :: t_arr(:)          !< 时间序列 (a.u.)
        complex(dp), intent(in) :: c_arr(:)       !< 自相关函数 C(t) (n_t)
        real(dp), intent(in) :: e0_au             !< 初态基态能量 E_0 (a.u.)
        real(dp), intent(in) :: omega_arr(:)      !< 目标光子频率数组 (a.u.)
        real(dp), intent(out) :: spectrum(:)      !< 相对吸收截面谱 (n_w)
        real(dp), intent(in), optional :: dt

        integer :: n_t, n_w, it, iw
        real(dp) :: dt_val, t_val, t_max, w_hann, total_energy
        complex(dp) :: integ, integrand

        n_t = size(c_arr)
        n_w = size(omega_arr)

        if (present(dt)) then
            dt_val = dt
        else if (n_t > 1) then
            dt_val = t_arr(2) - t_arr(1)
        else
            dt_val = 1.0_dp
        end if

        t_max = t_arr(n_t)
        if (t_max <= 0.0_dp) t_max = real(n_t, dp) * dt_val

        do iw = 1, n_w
            integ = (0.0_dp, 0.0_dp)
            total_energy = e0_au + omega_arr(iw)

            do it = 1, n_t
                t_val = t_arr(it)
                ! Hanning 平滑时间窗消除有限时域截断振荡
                w_hann = 0.5_dp * (1.0_dp + cos(PI * t_val / t_max))
                integrand = c_arr(it) * w_hann * exp(EYE * total_energy * t_val)
                if (it == 1 .or. it == n_t) then
                    integ = integ + 0.5_dp * integrand * dt_val
                else
                    integ = integ + integrand * dt_val
                end if
            end do

            spectrum(iw) = max(0.0_dp, omega_arr(iw) * real(integ, dp))
        end do
    end subroutine calculate_absorption_spectrum

    !> \brief 在渐近解离面 R = R_d 处计算时域波包采样的时间-能量傅里叶振幅 A(E)
    !> \details A(E) = (1/sqrt(2*pi)) * \int_0^T psi(R_d, t) * W(t) * exp(i*E*t) dt
    subroutine energy_resolved_flux_amplitude(t_arr, psi_rd_t, e_grid, amp_e, dt)
        real(dp), intent(in) :: t_arr(:)        !< 时间序列 (a.u.)
        complex(dp), intent(in) :: psi_rd_t(:)  !< 在渐近面 R=R_d 记录的含时波包 (n_t)
        real(dp), intent(in) :: e_grid(:)       !< 动能/总能网格 (n_e)
        complex(dp), intent(out) :: amp_e(:)    !< 能量分辨连续谱振幅 A(E) (n_e)
        real(dp), intent(in), optional :: dt

        integer :: n_t, n_e, it, ie
        real(dp) :: dt_val, t_val, t_max, w_hann, factor
        complex(dp) :: integ, term

        n_t = size(psi_rd_t)
        n_e = size(e_grid)

        if (present(dt)) then
            dt_val = dt
        else if (n_t > 1) then
            dt_val = t_arr(2) - t_arr(1)
        else
            dt_val = 1.0_dp
        end if

        t_max = max(1.0_dp, t_arr(n_t))
        factor = 1.0_dp / sqrt(TWOPI)

        do ie = 1, n_e
            integ = (0.0_dp, 0.0_dp)
            do it = 1, n_t
                t_val = t_arr(it)
                w_hann = 0.5_dp * (1.0_dp + cos(PI * t_val / t_max))
                term = psi_rd_t(it) * w_hann * exp(EYE * e_grid(ie) * t_val)
                if (it == 1 .or. it == n_t) then
                    integ = integ + 0.5_dp * term * dt_val
                else
                    integ = integ + term * dt_val
                end if
            end do
            amp_e(ie) = factor * integ
        end do
    end subroutine energy_resolved_flux_amplitude

    !> \brief 计算光碎片动能释放谱 (Kinetic Energy Release, KER) P(E_k)
    !> \details 基于能量分辨通量振幅 A(E): P(E_k) = (k/m) * |A(E_k)|^2 = sqrt(2*E_k/m) * |A(E_k)|^2
    subroutine calculate_ker_spectrum(amp_e, e_kin_grid, mass, ker_spectrum)
        complex(dp), intent(in) :: amp_e(:)     !< 能量连续态傅里叶振幅 (n_e)
        real(dp), intent(in) :: e_kin_grid(:)   !< 碎片相对动能网格 E_k (a.u.)
        real(dp), intent(in) :: mass            !< 体系约化质量 (a.u.)
        real(dp), intent(out) :: ker_spectrum(:)!< 动能释放能谱 P(E_k) (n_e)

        integer :: n_e, ie
        real(dp) :: vel

        n_e = size(amp_e)
        do ie = 1, n_e
            if (e_kin_grid(ie) > 0.0_dp .and. mass > 0.0_dp) then
                vel = sqrt(2.0_dp * e_kin_grid(ie) / mass)
                ker_spectrum(ie) = vel * (abs(amp_e(ie))**2)
            else
                ker_spectrum(ie) = 0.0_dp
            end if
        end do
    end subroutine calculate_ker_spectrum

    !> \brief 计算多通道解离碎片总产额与分支比 (Branching Ratios)
    subroutine calculate_branching_ratios(ker_channels, e_grid, yields, ratios)
        real(dp), intent(in) :: ker_channels(:, :)  !< 各通道能谱 (n_e, n_channels)
        real(dp), intent(in) :: e_grid(:)           !< 动能格点 (n_e)
        real(dp), intent(out) :: yields(:)          !< 各通道积分产额 (n_channels)
        real(dp), intent(out) :: ratios(:)          !< 各通道分支百分比 [0, 1] (n_channels)

        integer :: n_e, n_ch, ic
        real(dp) :: de, total_yield

        n_e = size(e_grid)
        n_ch = size(ker_channels, 2)

        if (n_e > 1) then
            de = e_grid(2) - e_grid(1)
        else
            de = 1.0_dp
        end if

        yields = 0.0_dp
        do ic = 1, n_ch
            yields(ic) = sum(ker_channels(:, ic)) * de
        end do

        total_yield = sum(yields)
        if (total_yield > 1.0e-30_dp) then
            ratios = yields / total_yield
        else
            ratios = 0.0_dp
        end if
    end subroutine calculate_branching_ratios

end module mod_photofragment_flux
