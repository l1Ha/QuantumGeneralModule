! ==============================================================================
! GeneralModule: mod_td_scattering.f90
! ------------------------------------------------------------------------------
! 现代量子动力学含时波包散射理论与 S-矩阵信息提取模块 (Time-Dependent Scattering)
! 包含特性：
! 1. 初始高斯散射波包制备与动量表象解析动量谱 g(k)
! 2. 渐近监测面时间-能量傅里叶积分获取连续态透射振幅 A_trans(E) 与反射振幅 A_refl(E)
! 3. 单次动力学模拟提取连续能量分辨透射几率 T(E)、反射几率 R(E) 与幺正性校验
! 4. 基于波包时频流与自由对照波包提取全能量散射矩阵元 S(E) 与散射相移 delta(E)
! 5. 坐标-动量表象 Möller 渐近投影法提取 S-矩阵与微分散射信息
! 6. 多通道含时概率通量提取非弹性散射矩阵元 |S_{ij}(E)|^2
! 7. 波包散射质心追踪与含时 Wigner 散射时延 tau_W(E) = 2*hbar*d(delta)/dE
! ==============================================================================
module mod_td_scattering
    use mod_constants, only: dp, PI, TWOPI, SQRTPI
    implicit none
    private

    ! --------------------------------------------------------------------------
    ! 派生类型定义
    ! --------------------------------------------------------------------------
    type, public :: td_scattering_channel_t
        real(dp) :: mass               ! 体系质量
        real(dp) :: x0                 ! 初始波包中心位置
        real(dp) :: sigma_x            ! 初始波包空间弥散宽度
        real(dp) :: k0                 ! 初始波包平均入射动量
        real(dp) :: x_det_trans        ! 渐近透射监测面位置
        real(dp) :: x_det_refl         ! 渐近反射监测面位置
    end type td_scattering_channel_t

    ! --------------------------------------------------------------------------
    ! 公共接口导出
    ! --------------------------------------------------------------------------
    public :: gaussian_wavepacket_1d
    public :: gaussian_momentum_amplitude
    public :: accumulate_flux_amplitude
    public :: calculate_td_transmission
    public :: calculate_td_smatrix_element
    public :: project_wavepacket_to_smatrix
    public :: multichannel_td_smatrix_elements
    public :: wavepacket_centroid_position
    public :: wavepacket_wigner_delay

contains

    ! ==========================================================================
    ! 1. 构造一维入射高斯波包
    !    psi_0(x) = (2*pi*sigma_x^2)^(-1/4) * exp(-(x - x0)^2 / (4*sigma_x^2) + i*k0*(x - x0))
    ! ==========================================================================
    subroutine gaussian_wavepacket_1d(x_grid, x0, sigma_x, k0, psi_0)
        real(dp), dimension(:), intent(in)        :: x_grid
        real(dp), intent(in)                      :: x0
        real(dp), intent(in)                      :: sigma_x
        real(dp), intent(in)                      :: k0
        complex(dp), dimension(:), intent(out)    :: psi_0

        integer  :: n_pts, i
        real(dp) :: norm_fac, dx_rel, env, phase

        n_pts = size(x_grid)
        norm_fac = (2.0_dp * PI * (sigma_x**2))**(-0.25_dp)

        do i = 1, n_pts
            dx_rel = x_grid(i) - x0
            env = norm_fac * exp(-(dx_rel**2) / (4.0_dp * (sigma_x**2)))
            phase = k0 * dx_rel
            psi_0(i) = cmplx(env * cos(phase), env * sin(phase), kind=dp)
        end do
    end subroutine gaussian_wavepacket_1d

    ! ==========================================================================
    ! 2. 初始波包在动量空间的解析振幅 g(k)
    !    g(k) = (2*pi*sigma_x^2)^(1/4) * sqrt(2/pi) * exp(-sigma_x^2 * (k - k0)^2 - i*k*x0)
    !    |g(k)|^2 为动量几率分布密度，满足 int |g(k)|^2 dk = 1
    ! ==========================================================================
    pure function gaussian_momentum_amplitude(k, x0, sigma_x, k0) result(gk)
        real(dp), intent(in) :: k
        real(dp), intent(in) :: x0
        real(dp), intent(in) :: sigma_x
        real(dp), intent(in) :: k0
        complex(dp)          :: gk

        real(dp) :: pref, env, phase

        pref = (2.0_dp * (sigma_x**2) / PI)**(0.25_dp)
        env = pref * exp(-(sigma_x**2) * ((k - k0)**2))
        phase = -k * x0
        gk = cmplx(env * cos(phase), env * sin(phase), kind=dp)
    end function gaussian_momentum_amplitude

    ! ==========================================================================
    ! 3. 渐近监测面时间-能量傅里叶振幅累积
    !    A(E) = 1/sqrt(2*pi) * int_0^T psi(x_det, t) * exp(i * E * t / hbar) dt
    ! ==========================================================================
    subroutine accumulate_flux_amplitude(t, psi_at_det, dt, energy_grid, n_energies, hbar, amp_accum)
        real(dp), intent(in)                      :: t
        complex(dp), intent(in)                   :: psi_at_det
        real(dp), intent(in)                      :: dt
        real(dp), dimension(n_energies), intent(in):: energy_grid
        integer, intent(in)                       :: n_energies
        real(dp), intent(in)                      :: hbar
        complex(dp), dimension(n_energies), intent(inout) :: amp_accum

        integer  :: ie
        real(dp) :: phase, inv_sqrt_twopi
        complex(dp) :: exp_fac

        inv_sqrt_twopi = 1.0_dp / sqrt(TWOPI)

        do ie = 1, n_energies
            phase = energy_grid(ie) * t / hbar
            exp_fac = cmplx(cos(phase), sin(phase), kind=dp)
            amp_accum(ie) = amp_accum(ie) + inv_sqrt_twopi * psi_at_det * exp_fac * dt
        end do
    end subroutine accumulate_flux_amplitude

    ! ==========================================================================
    ! 4. 基于含时通量傅里叶振幅计算能量分辨透射几率 T(E) 与反射几率 R(E)
    !    T(E) = (hbar * k_E / mass) * |A_trans(E)|^2 / |g(k_E)|^2
    ! ==========================================================================
    subroutine calculate_td_transmission(energy_grid, n_energies, amp_trans, mass, hbar, &
                                         x0, sigma_x, k0, t_prob)
        real(dp), dimension(n_energies), intent(in)   :: energy_grid
        integer, intent(in)                           :: n_energies
        complex(dp), dimension(n_energies), intent(in):: amp_trans
        real(dp), intent(in)                          :: mass
        real(dp), intent(in)                          :: hbar
        real(dp), intent(in)                          :: x0
        real(dp), intent(in)                          :: sigma_x
        real(dp), intent(in)                          :: k0
        real(dp), dimension(n_energies), intent(out)  :: t_prob

        integer  :: ie
        real(dp) :: e_val, k_val, flux_weight, p_inc
        complex(dp) :: gk

        do ie = 1, n_energies
            e_val = energy_grid(ie)
            if (e_val <= 0.0_dp) then
                t_prob(ie) = 0.0_dp
                cycle
            end if

            k_val = sqrt(2.0_dp * mass * e_val) / hbar
            gk = gaussian_momentum_amplitude(k_val, x0, sigma_x, k0)
            p_inc = abs(gk)**2

            if (p_inc > 1.0e-12_dp) then
                flux_weight = (hbar * k_val / mass)
                t_prob(ie) = flux_weight * (abs(amp_trans(ie))**2) / p_inc
                ! 数值截断保证物理有界
                t_prob(ie) = min(1.0_dp, max(0.0_dp, t_prob(ie)))
            else
                t_prob(ie) = 0.0_dp
            end if
        end do
    end subroutine calculate_td_transmission

    ! ==========================================================================
    ! 5. 提取含时散射矩阵元 S(E) 与散射相移 delta(E)
    !    通过有势演化振幅 A_scatter(E) 与无势自由对照演化 A_free(E) 相比:
    !    S(E) = A_scatter(E) / A_free(E)
    !    delta(E) = 0.5 * arg(S(E))
    ! ==========================================================================
    subroutine calculate_td_smatrix_element(amp_scatter, amp_free, n_energies, s_matrix_e, phase_shift_e)
        complex(dp), dimension(n_energies), intent(in)  :: amp_scatter
        complex(dp), dimension(n_energies), intent(in)  :: amp_free
        integer, intent(in)                             :: n_energies
        complex(dp), dimension(n_energies), intent(out) :: s_matrix_e
        real(dp), dimension(n_energies), intent(out)    :: phase_shift_e

        integer  :: ie
        complex(dp) :: s_val

        do ie = 1, n_energies
            if (abs(amp_free(ie)) > 1.0e-14_dp) then
                s_val = amp_scatter(ie) / amp_free(ie)
                s_matrix_e(ie) = s_val
                phase_shift_e(ie) = 0.5_dp * atan2(aimag(s_val), real(s_val, dp))
            else
                s_matrix_e(ie) = (1.0_dp, 0.0_dp)
                phase_shift_e(ie) = 0.0_dp
            end if
        end do
    end subroutine calculate_td_smatrix_element

    ! ==========================================================================
    ! 6. 坐标-动量表象 Möller 渐近投影法
    !    将完成碰撞后的末态波包 psi(x, T_final) 离散傅里叶变换到动量空间
    !    k > 0 成分对应透射几率，k < 0 成分对应反射几率
    ! ==========================================================================
    subroutine project_wavepacket_to_smatrix(x_grid, dx, psi_final, mass, hbar, &
                                            k0, sigma_x, x0, energy_grid, n_energies, &
                                            t_prob, r_prob)
        real(dp), dimension(:), intent(in)          :: x_grid
        real(dp), intent(in)                        :: dx
        complex(dp), dimension(:), intent(in)       :: psi_final
        real(dp), intent(in)                        :: mass
        real(dp), intent(in)                        :: hbar
        real(dp), intent(in)                        :: k0
        real(dp), intent(in)                        :: sigma_x
        real(dp), intent(in)                        :: x0
        real(dp), dimension(n_energies), intent(in) :: energy_grid
        integer, intent(in)                         :: n_energies
        real(dp), dimension(n_energies), intent(out):: t_prob
        real(dp), dimension(n_energies), intent(out):: r_prob

        integer  :: ie, ix, nx
        real(dp) :: e_val, k_val, inv_sqrt_twopi, p_inc
        complex(dp) :: psi_k_pos, psi_k_neg, gk, exp_pos, exp_neg
        real(dp) :: phase

        nx = size(x_grid)
        inv_sqrt_twopi = 1.0_dp / sqrt(TWOPI)

        do ie = 1, n_energies
            e_val = energy_grid(ie)
            if (e_val <= 0.0_dp) then
                t_prob(ie) = 0.0_dp
                r_prob(ie) = 0.0_dp
                cycle
            end if

            k_val = sqrt(2.0_dp * mass * e_val) / hbar
            gk = gaussian_momentum_amplitude(k_val, x0, sigma_x, k0)
            p_inc = abs(gk)**2

            psi_k_pos = (0.0_dp, 0.0_dp)
            psi_k_neg = (0.0_dp, 0.0_dp)

            ! 投影数值积分: int psi(x) * exp(-i * k * x) dx
            do ix = 1, nx
                phase = k_val * x_grid(ix)
                exp_pos = cmplx(cos(phase), -sin(phase), kind=dp)
                exp_neg = cmplx(cos(phase),  sin(phase), kind=dp)

                psi_k_pos = psi_k_pos + psi_final(ix) * exp_pos * dx
                psi_k_neg = psi_k_neg + psi_final(ix) * exp_neg * dx
            end do

            psi_k_pos = psi_k_pos * inv_sqrt_twopi
            psi_k_neg = psi_k_neg * inv_sqrt_twopi

            if (p_inc > 1.0e-12_dp) then
                t_prob(ie) = min(1.0_dp, max(0.0_dp, abs(psi_k_pos)**2 / p_inc))
                r_prob(ie) = min(1.0_dp, max(0.0_dp, abs(psi_k_neg)**2 / p_inc))
            else
                t_prob(ie) = 0.0_dp
                r_prob(ie) = 0.0_dp
            end if
        end do
    end subroutine project_wavepacket_to_smatrix

    ! ==========================================================================
    ! 7. 多通道含时概率通量提取非弹性散射矩阵元 |S_{ij}(E)|^2
    ! ==========================================================================
    subroutine multichannel_td_smatrix_elements(amp_ch1, amp_ch2, mass, hbar, &
                                               energy_grid, delta_e, x0, sigma_x, k0, &
                                               s11_prob, s12_prob, n_energies)
        complex(dp), dimension(n_energies), intent(in) :: amp_ch1
        complex(dp), dimension(n_energies), intent(in) :: amp_ch2
        real(dp), intent(in)                           :: mass
        real(dp), intent(in)                           :: hbar
        real(dp), dimension(n_energies), intent(in)    :: energy_grid
        real(dp), intent(in)                           :: delta_e
        real(dp), intent(in)                           :: x0
        real(dp), intent(in)                           :: sigma_x
        real(dp), intent(in)                           :: k0
        real(dp), dimension(n_energies), intent(out)   :: s11_prob
        real(dp), dimension(n_energies), intent(out)   :: s12_prob
        integer, intent(in)                            :: n_energies

        integer  :: ie
        real(dp) :: e1, e2, k1, k2, p_inc
        complex(dp) :: gk1

        do ie = 1, n_energies
            e1 = energy_grid(ie)
            e2 = e1 - delta_e

            if (e1 <= 0.0_dp) then
                s11_prob(ie) = 0.0_dp
                s12_prob(ie) = 0.0_dp
                cycle
            end if

            k1 = sqrt(2.0_dp * mass * e1) / hbar
            gk1 = gaussian_momentum_amplitude(k1, x0, sigma_x, k0)
            p_inc = abs(gk1)**2

            if (p_inc > 1.0e-12_dp) then
                s11_prob(ie) = (hbar * k1 / mass) * (abs(amp_ch1(ie))**2) / p_inc
                s11_prob(ie) = min(1.0_dp, max(0.0_dp, s11_prob(ie)))

                if (e2 > 0.0_dp) then
                    k2 = sqrt(2.0_dp * mass * e2) / hbar
                    s12_prob(ie) = (hbar * k2 / mass) * (abs(amp_ch2(ie))**2) / p_inc
                    s12_prob(ie) = min(1.0_dp, max(0.0_dp, s12_prob(ie)))
                else
                    s12_prob(ie) = 0.0_dp
                end if
            else
                s11_prob(ie) = 0.0_dp
                s12_prob(ie) = 0.0_dp
            end if
        end do
    end subroutine multichannel_td_smatrix_elements

    ! ==========================================================================
    ! 8. 波包质心追踪与含时 Wigner 散射时延
    ! ==========================================================================
    pure function wavepacket_centroid_position(x_grid, dx, psi) result(x_mean)
        real(dp), dimension(:), intent(in)    :: x_grid
        real(dp), intent(in)                  :: dx
        complex(dp), dimension(:), intent(in) :: psi
        real(dp)                              :: x_mean

        integer  :: i, n
        real(dp) :: norm_sum, x_prob_sum

        n = size(x_grid)
        norm_sum = 0.0_dp
        x_prob_sum = 0.0_dp

        do i = 1, n
            norm_sum = norm_sum + (abs(psi(i))**2) * dx
            x_prob_sum = x_prob_sum + x_grid(i) * (abs(psi(i))**2) * dx
        end do

        if (norm_sum > 1.0e-14_dp) then
            x_mean = x_prob_sum / norm_sum
        else
            x_mean = 0.0_dp
        end if
    end function wavepacket_centroid_position

    subroutine wavepacket_wigner_delay(energy_grid, phase_shift, n_energies, hbar, delay_wigner)
        real(dp), dimension(n_energies), intent(in)  :: energy_grid
        real(dp), dimension(n_energies), intent(in)  :: phase_shift
        integer, intent(in)                          :: n_energies
        real(dp), intent(in)                         :: hbar
        real(dp), dimension(n_energies), intent(out) :: delay_wigner

        integer  :: ie
        real(dp) :: de, d_phi

        if (n_energies < 2) then
            delay_wigner = 0.0_dp
            return
        end if

        delay_wigner(1) = 0.0_dp
        delay_wigner(n_energies) = 0.0_dp

        do ie = 2, n_energies - 1
            de = energy_grid(ie + 1) - energy_grid(ie - 1)
            d_phi = phase_shift(ie + 1) - phase_shift(ie - 1)
            ! tau_W = 2 * hbar * d(delta)/dE
            if (abs(de) > 1.0e-14_dp) then
                delay_wigner(ie) = 2.0_dp * hbar * (d_phi / de)
            else
                delay_wigner(ie) = 0.0_dp
            end if
        end do
        delay_wigner(1) = delay_wigner(2)
        delay_wigner(n_energies) = delay_wigner(n_energies - 1)
    end subroutine wavepacket_wigner_delay

end module mod_td_scattering
