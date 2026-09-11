! ==============================================================================
! GeneralModule: mod_ti_scattering.f90
! ------------------------------------------------------------------------------
! 现代量子动力学非含时散射理论与超冷碰撞模块 (Time-Independent Scattering)
! 包含特性：
! 1. Riccati-Bessel/Neumann 特殊函数与任意角动量 l 的解析导数
! 2. 零能 Numerov 与 Johnson Log-Derivative 散射长度 (a_s) 求解
! 3. 有限正能量分波相移 delta_l(E)、反应矩阵 K_l、散射矩阵 S_l 与跃迁矩阵 T_l
! 4. 分波弹性截面、总截面与光学定理一致性自洽校验
! 5. 角度分辨微分散射截面 d(sigma)/d(Omega) (Legendre 多项式级数展开)
! 6. 超低能区有效力程展开 (ERE: a_s, r_0) 最小二乘拟合
! 7. 范德华长程色散平均散射长度 a_bar 与 Gribakin-Flambaum 半经典解析公式
! 8. 形状共振 (Shape Resonance) Wigner 时延分析与 Breit-Wigner 参数提取
! 9. 双通道非绝热耦合密耦定态散射矩阵 (Close-Coupling 2x2 S-Matrix)
! ==============================================================================
module mod_ti_scattering
    use mod_constants, only: dp, PI, TWOPI, HALFPI
    use mod_special_functions, only: legendre_poly
    implicit none
    private

    ! --------------------------------------------------------------------------
    ! 派生类型定义
    ! --------------------------------------------------------------------------
    type, public :: scattering_state_t
        real(dp)    :: energy          ! 碰撞质心动能 E (a.u.)
        real(dp)    :: k_wave          ! 渐近相对动量 k = sqrt(2*mu*E)/hbar
        integer     :: l               ! 轨道角动量量子数
        real(dp)    :: phase_shift     ! 分波相移 delta_l (rad)
        real(dp)    :: k_matrix        ! 反应矩阵元 K_l = tan(delta_l)
        complex(dp) :: s_matrix        ! 散射矩阵元 S_l = exp(2*i*delta_l)
        complex(dp) :: t_matrix        ! 跃迁矩阵元 T_l = S_l - 1
        real(dp)    :: cross_section   ! 分波弹性截面 sigma_l
    end type scattering_state_t

    type, public :: ere_result_t
        real(dp) :: a_s                ! s-波散射长度
        real(dp) :: r_0                ! 有效力程 (effective range)
        real(dp) :: sigma_zero         ! 零温极限弹性散射截面 4*pi*a_s^2
    end type ere_result_t

    type, public :: resonance_info_t
        real(dp) :: e_res              ! 共振能量 E_R (a.u.)
        real(dp) :: gamma_width        ! 共振线宽 Gamma (a.u.)
        real(dp) :: tau_max            ! 最大 Wigner 时延 (a.u.)
        real(dp) :: lifetime           ! 准束缚态寿命 hbar / Gamma
        real(dp) :: peak_cross_section ! 峰值弹性散射截面
    end type resonance_info_t

    ! 全同粒子统计与量子对称性参数
    integer, parameter, public :: PARTICLE_DISTINGUISHABLE             = 0
    integer, parameter, public :: PARTICLE_IDENTICAL_BOSON             = 1
    integer, parameter, public :: PARTICLE_IDENTICAL_FERMION_POLARIZED = 2
    integer, parameter, public :: PARTICLE_IDENTICAL_FERMION_UNPOLAR   = 3

    ! --------------------------------------------------------------------------
    ! 公共接口导出
    ! --------------------------------------------------------------------------
    public :: riccati_bessel_neumann
    public :: calc_scattering_length_numerov
    public :: calc_scattering_length_logder
    public :: calc_phase_shift_single_l
    public :: calc_partial_wave_cross_sections
    public :: optical_theorem_cross_section
    public :: calc_differential_cross_section
    public :: calc_differential_cross_section_identical
    public :: calc_transport_cross_sections
    public :: calc_cross_section_spectrum
    public :: calc_generalized_cross_sections
    public :: calc_differential_legendre_expansion
    public :: fit_effective_range_expansion
    public :: gribakin_flambaum_length
    public :: van_der_waals_mean_length
    public :: analyze_shape_resonance
    public :: calc_coupled_channel_smatrix_2x2

contains

    ! ==========================================================================
    ! 1. Riccati-Bessel 与 Riccati-Neumann 函数及其导数
    !    j_hat_l(x) = x * j_l(x), n_hat_l(x) = x * n_l(x)
    ! ==========================================================================
    pure subroutine riccati_bessel_neumann(l, x, jl, nl, d_jl, d_nl)
        integer, intent(in)   :: l
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: jl, nl, d_jl, d_nl

        real(dp) :: j_prev, j_curr, j_next
        real(dp) :: n_prev, n_curr, n_next
        real(dp) :: sin_x, cos_x
        integer  :: k

        if (x <= 1.0e-14_dp) then
            if (l == 0) then
                jl = 0.0_dp; nl = -1.0_dp; d_jl = 1.0_dp; d_nl = 0.0_dp
            else
                jl = 0.0_dp; nl = -1.0e30_dp; d_jl = 0.0_dp; d_nl = 1.0e30_dp
            end if
            return
        end if

        sin_x = sin(x)
        cos_x = cos(x)

        if (l == 0) then
            jl   = sin_x
            nl   = -cos_x
            d_jl = cos_x
            d_nl = sin_x
            return
        else if (l == 1) then
            jl   = sin_x / x - cos_x
            nl   = -cos_x / x - sin_x
            d_jl = cos_x / x - sin_x / (x * x) + sin_x
            d_nl = sin_x / x + cos_x / (x * x) - cos_x
            return
        end if

        ! 对于 l >= 2, 利用递推关系: f_{k+1}(x) = (2k+1)/x * f_k(x) - f_{k-1}(x)
        j_prev = sin_x
        j_curr = sin_x / x - cos_x

        n_prev = -cos_x
        n_curr = -cos_x / x - sin_x

        do k = 1, l - 1
            j_next = real(2 * k + 1, dp) / x * j_curr - j_prev
            n_next = real(2 * k + 1, dp) / x * n_curr - n_prev
            j_prev = j_curr
            j_curr = j_next
            n_prev = n_curr
            n_curr = n_next
        end do

        jl = j_curr
        nl = n_curr

        ! 导数关系: d/dx f_l(x) = f_{l-1}(x) - l/x * f_l(x)
        d_jl = j_prev - real(l, dp) / x * jl
        d_nl = n_prev - real(l, dp) / x * nl
    end subroutine riccati_bessel_neumann

    ! ==========================================================================
    ! 2. 零能 Numerov 算法计算 s-波散射长度 a_s
    !    求解: -hbar^2/(2mu) u''(r) + V(r) u(r) = 0
    !    在渐近区: u(r) -> C * (r - a_s) => a_s = r - u(r)/u'(r)
    ! ==========================================================================
    subroutine calc_scattering_length_numerov(r_grid, v_pot, mass, a_s, u_zero, stat)
        real(dp), dimension(:), intent(in)  :: r_grid
        real(dp), dimension(:), intent(in)  :: v_pot
        real(dp), intent(in)                :: mass
        real(dp), intent(out)               :: a_s
        real(dp), dimension(:), optional, intent(out) :: u_zero
        integer, optional, intent(out)      :: stat

        integer  :: n_pts, i
        real(dp) :: dr, dr2_12, q_curr, q_prev, q_next
        real(dp) :: c_curr, c_prev, c_next, d_u
        real(dp), allocatable :: u_wf(:)

        if (present(stat)) stat = 0
        n_pts = size(r_grid)
        if (n_pts < 5) then
            if (present(stat)) stat = -1
            a_s = 0.0_dp
            return
        end if

        dr = r_grid(2) - r_grid(1)
        dr2_12 = (dr * dr) / 12.0_dp

        allocate(u_wf(n_pts))
        u_wf = 0.0_dp

        ! 禁区原点初始边界条件
        u_wf(1) = 0.0_dp
        u_wf(2) = 1.0e-7_dp  ! 任意极小线性非零初值

        ! Numerov 推进: (1 + h^2/12 Q_{i+1}) u_{i+1} = 2(1 - 5h^2/12 Q_i) u_i - (1 + h^2/12 Q_{i-1}) u_{i-1}
        ! 此处 Q(r) = -2*mu/hbar^2 * V(r)
        do i = 2, n_pts - 1
            q_prev = -2.0_dp * mass * v_pot(i - 1)
            q_curr = -2.0_dp * mass * v_pot(i)
            q_next = -2.0_dp * mass * v_pot(i + 1)

            c_prev = 1.0_dp + dr2_12 * q_prev
            c_curr = 2.0_dp * (1.0_dp - 5.0_dp * dr2_12 * q_curr)
            c_next = 1.0_dp + dr2_12 * q_next

            u_wf(i + 1) = (c_curr * u_wf(i) - c_prev * u_wf(i - 1)) / c_next

            ! 避免超长积分溢出
            if (abs(u_wf(i + 1)) > 1.0e20_dp) then
                u_wf(1:i + 1) = u_wf(1:i + 1) * 1.0e-15_dp
            end if
        end do

        ! 在渐近外边界计算对数导数截距:
        ! u'(r_N) = (3 u_N - 4 u_{N-1} + u_{N-2}) / (2 dr)
        d_u = (3.0_dp * u_wf(n_pts) - 4.0_dp * u_wf(n_pts - 1) + u_wf(n_pts - 2)) / (2.0_dp * dr)

        if (abs(d_u) > 1.0e-14_dp) then
            a_s = r_grid(n_pts) - u_wf(n_pts) / d_u
        else
            a_s = 1.0e30_dp  ! 刚好处于零能束缚态共振极点
        end if

        if (present(u_zero)) then
            if (size(u_zero) == n_pts) then
                u_zero = u_wf
            end if
        end if

        deallocate(u_wf)
    end subroutine calc_scattering_length_numerov

    ! ==========================================================================
    ! 3. Johnson 对数导数比值法 (Log-Derivative Ratio Propagator)
    !    直接推进比值 R_i = u_i / u_{i-1}，天然抗拒强深阱经典禁区的指数溢出
    !    渐近截距公式: a_s = r_N - dr / (1 - 1/R_N)
    ! ==========================================================================
    subroutine calc_scattering_length_logder(r_grid, v_pot, mass, a_s, stat)
        real(dp), dimension(:), intent(in) :: r_grid
        real(dp), dimension(:), intent(in) :: v_pot
        real(dp), intent(in)               :: mass
        real(dp), intent(out)              :: a_s
        integer, optional, intent(out)     :: stat

        integer  :: n_pts, i
        real(dp) :: dr, dr2_12
        real(dp) :: w_prev, w_curr, w_next, r_ratio

        if (present(stat)) stat = 0
        n_pts = size(r_grid)
        if (n_pts < 5) then
            if (present(stat)) stat = -1
            a_s = 0.0_dp
            return
        end if

        dr = r_grid(2) - r_grid(1)
        dr2_12 = (dr * dr) / 12.0_dp

        ! 起始边界比值
        r_ratio = 2.0_dp  ! 线性起始 u(2)/u(1)

        do i = 2, n_pts - 1
            w_prev = dr2_12 * (-2.0_dp * mass * v_pot(i - 1))
            w_curr = dr2_12 * (-2.0_dp * mass * v_pot(i))
            w_next = dr2_12 * (-2.0_dp * mass * v_pot(i + 1))

            ! Johnson 比值推进公式:
            ! R_{i+1} = [ 2(1 - 5 w_i) - (1 + w_{i-1})/R_i ] / (1 + w_{i+1})
            r_ratio = (2.0_dp * (1.0_dp - 5.0_dp * w_curr) - (1.0_dp + w_prev) / r_ratio) / (1.0_dp + w_next)
        end do

        ! 根据渐近线: u_N = C (r_N - a_s), u_{N-1} = C (r_{N-1} - a_s)
        ! R_N = (r_N - a_s) / (r_{N-1} - a_s) => a_s = (R_N r_{N-1} - r_N) / (R_N - 1)
        if (abs(r_ratio - 1.0_dp) > 1.0e-12_dp) then
            a_s = (r_ratio * r_grid(n_pts - 1) - r_grid(n_pts)) / (r_ratio - 1.0_dp)
        else
            a_s = 1.0e30_dp
        end if
    end subroutine calc_scattering_length_logder

    ! ==========================================================================
    ! 4. 有限正能量单分波相移 delta_l 与 S 矩阵求解
    !    在渐近外边界 R_match 与 Riccati-Bessel/Neumann 函数匹配
    ! ==========================================================================
    subroutine calc_phase_shift_single_l(r_grid, v_pot, mass, energy, l, delta, k_mat, s_mat, t_mat, stat)
        real(dp), dimension(:), intent(in) :: r_grid
        real(dp), dimension(:), intent(in) :: v_pot
        real(dp), intent(in)               :: mass
        real(dp), intent(in)               :: energy
        integer, intent(in)                :: l
        real(dp), intent(out)              :: delta
        real(dp), intent(out)              :: k_mat
        complex(dp), intent(out)           :: s_mat
        complex(dp), intent(out)           :: t_mat
        integer, optional, intent(out)     :: stat

        integer  :: n_pts, i
        real(dp) :: dr, dr2_12, k_wave, r_match
        real(dp) :: q_prev, q_curr, q_next
        real(dp) :: c_prev, c_curr, c_next, d_u, y_logder
        real(dp) :: jl, nl, d_jl, d_nl
        real(dp) :: num, den
        real(dp), allocatable :: u_wf(:)

        if (present(stat)) stat = 0
        n_pts = size(r_grid)

        if (energy <= 0.0_dp .or. n_pts < 5) then
            if (present(stat)) stat = -1
            delta = 0.0_dp; k_mat = 0.0_dp; s_mat = (1.0_dp, 0.0_dp); t_mat = (0.0_dp, 0.0_dp)
            return
        end if

        k_wave = sqrt(2.0_dp * mass * energy)
        dr = r_grid(2) - r_grid(1)
        dr2_12 = (dr * dr) / 12.0_dp

        allocate(u_wf(n_pts))
        u_wf = 0.0_dp
        u_wf(1) = 0.0_dp
        u_wf(2) = (dr)**(l + 1) * 1.0e-5_dp

        do i = 2, n_pts - 1
            q_prev = 2.0_dp * mass * (energy - v_pot(i - 1)) - real(l * (l + 1), dp) / (r_grid(i - 1)**2)
            q_curr = 2.0_dp * mass * (energy - v_pot(i))     - real(l * (l + 1), dp) / (r_grid(i)**2)
            q_next = 2.0_dp * mass * (energy - v_pot(i + 1)) - real(l * (l + 1), dp) / (r_grid(i + 1)**2)

            c_prev = 1.0_dp + dr2_12 * q_prev
            c_curr = 2.0_dp * (1.0_dp - 5.0_dp * dr2_12 * q_curr)
            c_next = 1.0_dp + dr2_12 * q_next

            u_wf(i + 1) = (c_curr * u_wf(i) - c_prev * u_wf(i - 1)) / c_next

            if (abs(u_wf(i + 1)) > 1.0e20_dp) then
                u_wf(1:i + 1) = u_wf(1:i + 1) * 1.0e-15_dp
            end if
        end do

        r_match = r_grid(n_pts)
        d_u = (3.0_dp * u_wf(n_pts) - 4.0_dp * u_wf(n_pts - 1) + u_wf(n_pts - 2)) / (2.0_dp * dr)
        y_logder = d_u / u_wf(n_pts)

        call riccati_bessel_neumann(l, k_wave * r_match, jl, nl, d_jl, d_nl)

        ! tan(delta_l) = (k * j'_l - y * j_l) / (k * n'_l - y * n_l)
        num = k_wave * d_jl - y_logder * jl
        den = k_wave * d_nl - y_logder * nl

        delta = atan2(num, den)
        k_mat = tan(delta)
        s_mat = cmplx(cos(2.0_dp * delta), sin(2.0_dp * delta), kind=dp)
        t_mat = s_mat - (1.0_dp, 0.0_dp)

        deallocate(u_wf)
    end subroutine calc_phase_shift_single_l

    ! ==========================================================================
    ! 5. 多分波截面与总碰撞截面计算
    ! ==========================================================================
    subroutine calc_partial_wave_cross_sections(r_grid, v_pot, mass, energy, l_max, &
                                              delta_arr, sigma_partial, sigma_total, stat)
        real(dp), dimension(:), intent(in)   :: r_grid
        real(dp), dimension(:), intent(in)   :: v_pot
        real(dp), intent(in)                 :: mass
        real(dp), intent(in)                 :: energy
        integer, intent(in)                  :: l_max
        real(dp), dimension(0:l_max), intent(out) :: delta_arr
        real(dp), dimension(0:l_max), intent(out) :: sigma_partial
        real(dp), intent(out)                :: sigma_total
        integer, optional, intent(out)       :: stat

        integer  :: l, istat
        real(dp) :: k_wave, k_m, s_m_re, s_m_im, t_m_re, t_m_im
        complex(dp) :: s_m, t_m

        if (present(stat)) stat = 0
        sigma_total = 0.0_dp

        if (energy <= 0.0_dp) then
            if (present(stat)) stat = -1
            delta_arr = 0.0_dp
            sigma_partial = 0.0_dp
            return
        end if

        k_wave = sqrt(2.0_dp * mass * energy)

        do l = 0, l_max
            call calc_phase_shift_single_l(r_grid, v_pot, mass, energy, l, &
                                           delta_arr(l), k_m, s_m, t_m, istat)
            ! 分波弹性截面: sigma_l = (4 * pi / k^2) * (2l + 1) * sin^2(delta_l)
            sigma_partial(l) = (4.0_dp * PI / (k_wave * k_wave)) * real(2 * l + 1, dp) * (sin(delta_arr(l))**2)
            sigma_total = sigma_total + sigma_partial(l)
        end do
    end subroutine calc_partial_wave_cross_sections

    ! ==========================================================================
    ! 6. 光学定理截面一致性校验
    !    sigma_optical = (4*pi/k) * Im[f(0)] = (4*pi/k^2) * sum_l (2l+1) sin^2(delta_l)
    ! ==========================================================================
    pure function optical_theorem_cross_section(k_wave, delta_arr, l_max) result(sigma_opt)
        real(dp), intent(in)                 :: k_wave
        real(dp), dimension(0:l_max), intent(in) :: delta_arr
        integer, intent(in)                  :: l_max
        real(dp)                             :: sigma_opt

        integer :: l

        sigma_opt = 0.0_dp
        if (k_wave <= 1.0e-14_dp) return

        do l = 0, l_max
            sigma_opt = sigma_opt + real(2 * l + 1, dp) * (sin(delta_arr(l))**2)
        end do
        sigma_opt = sigma_opt * (4.0_dp * PI / (k_wave * k_wave))
    end function optical_theorem_cross_section

    ! ==========================================================================
    ! 7. 微分散射截面 d(sigma)/d(Omega) (以 Legendre 多项式展开)
    !    f(theta) = 1/k * sum_l (2l+1) exp(i*delta_l) * sin(delta_l) * P_l(cos theta)
    !    d(sigma)/d(Omega) = |f(theta)|^2
    ! ==========================================================================
    subroutine calc_differential_cross_section(energy, mass, delta_arr, l_max, theta_grid, dsigma_domega)
        real(dp), intent(in)                 :: energy
        real(dp), intent(in)                 :: mass
        real(dp), dimension(0:l_max), intent(in) :: delta_arr
        integer, intent(in)                  :: l_max
        real(dp), dimension(:), intent(in)   :: theta_grid
        real(dp), dimension(:), intent(out)  :: dsigma_domega

        integer  :: n_angles, i, l
        real(dp) :: k_wave, x_cos, p_l
        complex(dp) :: f_theta, term

        n_angles = size(theta_grid)
        k_wave = sqrt(2.0_dp * mass * max(1.0e-14_dp, energy))

        do i = 1, n_angles
            x_cos = cos(theta_grid(i))
            f_theta = (0.0_dp, 0.0_dp)

            do l = 0, l_max
                p_l = legendre_poly(l, x_cos)
                term = real(2 * l + 1, dp) * cmplx(cos(delta_arr(l)), sin(delta_arr(l)), kind=dp) * &
                       sin(delta_arr(l)) * p_l
                f_theta = f_theta + term
            end do

            f_theta = f_theta / k_wave
            dsigma_domega(i) = abs(f_theta)**2
        end do
    end subroutine calc_differential_cross_section

    ! ==========================================================================
    ! 7b. 全同粒子量子对称性微分散射截面
    !     区分粒子: |f(theta)|^2
    !     全同玻色子: |f(theta) + f(pi - theta)|^2
    !     极化全同费米子: |f(theta) - f(pi - theta)|^2
    !     非极化自旋-1/2费米子: 1/4 |f(theta) + f(pi - theta)|^2 + 3/4 |f(theta) - f(pi - theta)|^2
    ! ==========================================================================
    subroutine calc_differential_cross_section_identical(energy, mass, delta_arr, l_max, &
                                                        theta_grid, particle_stat, dsigma_domega)
        real(dp), intent(in)                 :: energy
        real(dp), intent(in)                 :: mass
        real(dp), dimension(0:l_max), intent(in) :: delta_arr
        integer, intent(in)                  :: l_max
        real(dp), dimension(:), intent(in)   :: theta_grid
        integer, intent(in)                  :: particle_stat
        real(dp), dimension(:), intent(out)  :: dsigma_domega

        integer  :: n_angles, i, l
        real(dp) :: k_wave, x_cos, p_l, p_l_pi
        complex(dp) :: f_theta, f_pi_minus_theta, term, term_pi
        real(dp) :: ds_boson, ds_fermion

        n_angles = size(theta_grid)
        k_wave = sqrt(2.0_dp * mass * max(1.0e-14_dp, energy))

        do i = 1, n_angles
            x_cos = cos(theta_grid(i))
            f_theta = (0.0_dp, 0.0_dp)
            f_pi_minus_theta = (0.0_dp, 0.0_dp)

            do l = 0, l_max
                p_l = legendre_poly(l, x_cos)
                p_l_pi = ((-1.0_dp)**l) * p_l  ! P_l(-cos theta) = (-1)^l P_l(cos theta)

                term = real(2 * l + 1, dp) * cmplx(cos(delta_arr(l)), sin(delta_arr(l)), kind=dp) * &
                       sin(delta_arr(l)) * p_l
                term_pi = real(2 * l + 1, dp) * cmplx(cos(delta_arr(l)), sin(delta_arr(l)), kind=dp) * &
                          sin(delta_arr(l)) * p_l_pi

                f_theta = f_theta + term
                f_pi_minus_theta = f_pi_minus_theta + term_pi
            end do

            f_theta = f_theta / k_wave
            f_pi_minus_theta = f_pi_minus_theta / k_wave

            select case (particle_stat)
            case (PARTICLE_DISTINGUISHABLE)
                dsigma_domega(i) = abs(f_theta)**2
            case (PARTICLE_IDENTICAL_BOSON)
                dsigma_domega(i) = abs(f_theta + f_pi_minus_theta)**2
            case (PARTICLE_IDENTICAL_FERMION_POLARIZED)
                dsigma_domega(i) = abs(f_theta - f_pi_minus_theta)**2
            case (PARTICLE_IDENTICAL_FERMION_UNPOLAR)
                ds_boson   = abs(f_theta + f_pi_minus_theta)**2
                ds_fermion = abs(f_theta - f_pi_minus_theta)**2
                dsigma_domega(i) = 0.25_dp * ds_boson + 0.75_dp * ds_fermion
            case default
                dsigma_domega(i) = abs(f_theta)**2
            end select
        end do
    end subroutine calc_differential_cross_section_identical

    ! ==========================================================================
    ! 7c. 输运散射截面 (动量传输截面 sigma_m 与 粘滞截面 sigma_v)
    !     sigma_m = 4*pi/k^2 * sum_{l=0}^{l_max-1} (l+1) sin^2(delta_l - delta_{l+1})
    !     sigma_v = 4*pi/k^2 * sum_{l=0}^{l_max-2} (l+1)(l+2)/(2l+3) sin^2(delta_l - delta_{l+2})
    ! ==========================================================================
    subroutine calc_transport_cross_sections(energy, mass, delta_arr, l_max, sigma_momentum, sigma_viscosity)
        real(dp), intent(in)                 :: energy
        real(dp), intent(in)                 :: mass
        real(dp), dimension(0:l_max), intent(in) :: delta_arr
        integer, intent(in)                  :: l_max
        real(dp), intent(out)                :: sigma_momentum
        real(dp), intent(out)                :: sigma_viscosity

        integer  :: l
        real(dp) :: k_wave, pref

        sigma_momentum = 0.0_dp
        sigma_viscosity = 0.0_dp

        if (energy <= 1.0e-14_dp .or. l_max < 1) return

        k_wave = sqrt(2.0_dp * mass * energy)
        pref = 4.0_dp * PI / (k_wave * k_wave)

        ! 动量传输截面 (Diffusion cross section)
        do l = 0, l_max - 1
            sigma_momentum = sigma_momentum + real(l + 1, dp) * (sin(delta_arr(l) - delta_arr(l + 1))**2)
        end do
        sigma_momentum = sigma_momentum * pref

        ! 粘滞截面 (Viscosity cross section)
        if (l_max >= 2) then
            do l = 0, l_max - 2
                sigma_viscosity = sigma_viscosity + (real((l + 1) * (l + 2), dp) / real(2 * l + 3, dp)) * &
                                  (sin(delta_arr(l) - delta_arr(l + 2))**2)
            end do
            sigma_viscosity = sigma_viscosity * pref
        end if
    end subroutine calc_transport_cross_sections

    ! ==========================================================================
    ! 7d. 全能量范围散射截面能谱扫描 (展示低能常数平台、共振峰与 Ramsauer-Townsend 极小)
    ! ==========================================================================
    subroutine calc_cross_section_spectrum(r_grid, v_pot, mass, energy_grid, n_energies, &
                                          l_max, sigma_total, sigma_partial, stat)
        real(dp), dimension(:), intent(in)            :: r_grid
        real(dp), dimension(:), intent(in)            :: v_pot
        real(dp), intent(in)                          :: mass
        real(dp), dimension(n_energies), intent(in)   :: energy_grid
        integer, intent(in)                           :: n_energies
        integer, intent(in)                           :: l_max
        real(dp), dimension(n_energies), intent(out)  :: sigma_total
        real(dp), dimension(0:l_max, n_energies), optional, intent(out) :: sigma_partial
        integer, optional, intent(out)                :: stat

        integer  :: ie, istat
        real(dp) :: delta_arr(0:l_max), sig_part(0:l_max), sig_tot

        if (present(stat)) stat = 0

        do ie = 1, n_energies
            call calc_partial_wave_cross_sections(r_grid, v_pot, mass, energy_grid(ie), l_max, &
                                                  delta_arr, sig_part, sig_tot, istat)
            sigma_total(ie) = sig_tot
            if (present(sigma_partial)) then
                sigma_partial(0:l_max, ie) = sig_part
            end if
        end do
    end subroutine calc_cross_section_spectrum

    ! ==========================================================================
    ! 7e. 广义吸收/复势弹性截面、非弹性吸收截面与总截面
    !     sigma_el   = pi/k^2 * sum_l (2l+1) |1 - S_l|^2
    !     sigma_inel = pi/k^2 * sum_l (2l+1) (1 - |S_l|^2)
    !     sigma_tot  = sigma_el + sigma_inel = 2*pi/k^2 * sum_l (2l+1) (1 - Re S_l)
    ! ==========================================================================
    subroutine calc_generalized_cross_sections(k_wave, s_mat_arr, l_max, &
                                              sigma_elastic, sigma_inelastic, sigma_total)
        real(dp), intent(in)                     :: k_wave
        complex(dp), dimension(0:l_max), intent(in) :: s_mat_arr
        integer, intent(in)                      :: l_max
        real(dp), intent(out)                    :: sigma_elastic
        real(dp), intent(out)                    :: sigma_inelastic
        real(dp), intent(out)                    :: sigma_total

        integer  :: l
        real(dp) :: pref, deg, mod_s2

        sigma_elastic = 0.0_dp
        sigma_inelastic = 0.0_dp
        sigma_total = 0.0_dp

        if (k_wave <= 1.0e-14_dp) return

        pref = PI / (k_wave * k_wave)

        do l = 0, l_max
            deg = real(2 * l + 1, dp)
            mod_s2 = abs(s_mat_arr(l))**2

            sigma_elastic   = sigma_elastic   + deg * (abs(1.0_dp - s_mat_arr(l))**2)
            sigma_inelastic = sigma_inelastic + deg * max(0.0_dp, 1.0_dp - mod_s2)
            sigma_total     = sigma_total     + deg * 2.0_dp * (1.0_dp - real(s_mat_arr(l), dp))
        end do

        sigma_elastic   = sigma_elastic * pref
        sigma_inelastic = sigma_inelastic * pref
        sigma_total     = sigma_total * pref
    end subroutine calc_generalized_cross_sections

    ! ==========================================================================
    ! 7f. 微分散射截面 Legendre 级数展开与前后各向异性不对称参数 A_FB
    !     dsigma/dOmega = sum_{K} A_K P_K(cos theta)
    !     A_FB = (Forward - Backward) / (Forward + Backward)
    ! ==========================================================================
    subroutine calc_differential_legendre_expansion(theta_grid, dsigma_domega, k_max, a_coeff, fb_asymmetry)
        real(dp), dimension(:), intent(in)   :: theta_grid
        real(dp), dimension(:), intent(in)   :: dsigma_domega
        integer, intent(in)                  :: k_max
        real(dp), dimension(0:k_max), intent(out) :: a_coeff
        real(dp), intent(out)                :: fb_asymmetry

        integer  :: n_angles, i, k
        real(dp) :: dtheta, th, x_cos, p_k, sin_th, w
        real(dp) :: forward_flux, backward_flux

        n_angles = size(theta_grid)
        dtheta = theta_grid(2) - theta_grid(1)
        a_coeff = 0.0_dp
        forward_flux = 0.0_dp
        backward_flux = 0.0_dp

        do i = 1, n_angles
            th = theta_grid(i)
            x_cos = cos(th)
            sin_th = sin(th)
            w = sin_th * dtheta

            if (th <= HALFPI) then
                forward_flux = forward_flux + dsigma_domega(i) * w
            else
                backward_flux = backward_flux + dsigma_domega(i) * w
            end if

            do k = 0, k_max
                p_k = legendre_poly(k, x_cos)
                a_coeff(k) = a_coeff(k) + real(2 * k + 1, dp) * 0.5_dp * dsigma_domega(i) * p_k * w
            end do
        end do

        if (forward_flux + backward_flux > 1.0e-14_dp) then
            fb_asymmetry = (forward_flux - backward_flux) / (forward_flux + backward_flux)
        else
            fb_asymmetry = 0.0_dp
        end if
    end subroutine calc_differential_legendre_expansion

    ! ==========================================================================
    ! 8. 超低能区有效力程展开 (ERE: Effective Range Expansion) 最小二乘拟合
    !    k * cot(delta_0(k)) = -1/a_s + 0.5 * r_0 * k^2
    ! ==========================================================================
    subroutine fit_effective_range_expansion(r_grid, v_pot, mass, k_list, n_k, a_s, r_0, stat)
        real(dp), dimension(:), intent(in) :: r_grid
        real(dp), dimension(:), intent(in) :: v_pot
        real(dp), intent(in)               :: mass
        real(dp), dimension(n_k), intent(in) :: k_list
        integer, intent(in)                :: n_k
        real(dp), intent(out)              :: a_s
        real(dp), intent(out)              :: r_0
        integer, optional, intent(out)     :: stat

        integer  :: i, istat
        real(dp) :: e_val, delta, k_m
        complex(dp) :: s_m, t_m
        real(dp), dimension(n_k) :: x_vec, y_vec
        real(dp) :: sum_x, sum_y, sum_xx, sum_xy, det, a_fit, b_fit

        if (present(stat)) stat = 0
        if (n_k < 2) then
            if (present(stat)) stat = -1
            a_s = 0.0_dp; r_0 = 0.0_dp
            return
        end if

        do i = 1, n_k
            e_val = (k_list(i)**2) / (2.0_dp * mass)
            call calc_phase_shift_single_l(r_grid, v_pot, mass, e_val, 0, delta, k_m, s_m, t_m, istat)

            x_vec(i) = k_list(i)**2
            y_vec(i) = k_list(i) / tan(delta)  ! k * cot(delta)
        end do

        ! 线性回归: y = a_fit + b_fit * x => -1/a_s + 0.5 * r_0 * k^2
        sum_x = sum(x_vec)
        sum_y = sum(y_vec)
        sum_xx = sum(x_vec * x_vec)
        sum_xy = sum(x_vec * y_vec)

        det = real(n_k, dp) * sum_xx - sum_x * sum_x
        if (abs(det) < 1.0e-14_dp) then
            if (present(stat)) stat = -2
            a_s = 0.0_dp; r_0 = 0.0_dp
            return
        end if

        a_fit = (sum_y * sum_xx - sum_x * sum_xy) / det
        b_fit = (real(n_k, dp) * sum_xy - sum_x * sum_y) / det

        if (abs(a_fit) > 1.0e-12_dp) then
            a_s = -1.0_dp / a_fit
        else
            a_s = 1.0e30_dp
        end if
        r_0 = 2.0_dp * b_fit
    end subroutine fit_effective_range_expansion

    ! ==========================================================================
    ! 9. 范德华平均散射长度 a_bar 与 Gribakin-Flambaum 半经典解析公式
    ! ==========================================================================
    pure function van_der_waals_mean_length(mass, c6_au) result(a_bar)
        real(dp), intent(in) :: mass
        real(dp), intent(in) :: c6_au
        real(dp)             :: a_bar

        ! a_bar = 0.477988812586 * (2 * mu * C6)^(1/4)
        a_bar = 0.47798881258618_dp * (2.0_dp * mass * c6_au)**(0.25_dp)
    end function van_der_waals_mean_length

    pure function gribakin_flambaum_length(mass, c6_au, phase_phi) result(a_s)
        real(dp), intent(in) :: mass
        real(dp), intent(in) :: c6_au
        real(dp), intent(in) :: phase_phi
        real(dp)             :: a_s

        real(dp) :: a_bar
        a_bar = van_der_waals_mean_length(mass, c6_au)
        ! a_s = a_bar * [1 - tan(Phi - pi/8)]
        a_s = a_bar * (1.0_dp - tan(phase_phi - PI / 8.0_dp))
    end function gribakin_flambaum_length

    ! ==========================================================================
    ! 10. 形状共振 (Shape Resonance) Wigner 时延分析与 Breit-Wigner 参数提取
    !     tau(E) = 2 * hbar * d(delta)/dE
    !     在共振处 tau_max = 4 * hbar / Gamma => Gamma = 4 * hbar / tau_max
    ! ==========================================================================
    subroutine analyze_shape_resonance(energy_grid, delta_grid, n_pts, hbar, res_info, stat)
        real(dp), dimension(n_pts), intent(in) :: energy_grid
        real(dp), dimension(n_pts), intent(in) :: delta_grid
        integer, intent(in)                    :: n_pts
        real(dp), intent(in)                   :: hbar
        type(resonance_info_t), intent(out)    :: res_info
        integer, optional, intent(out)         :: stat

        integer  :: i, max_idx
        real(dp) :: d_e, d_delta, tau_curr, tau_max

        if (present(stat)) stat = 0
        if (n_pts < 3) then
            if (present(stat)) stat = -1
            res_info%e_res = 0.0_dp
            res_info%gamma_width = 0.0_dp
            return
        end if

        tau_max = -1.0e30_dp
        max_idx = 2

        do i = 2, n_pts - 1
            d_e = energy_grid(i + 1) - energy_grid(i - 1)
            d_delta = delta_grid(i + 1) - delta_grid(i - 1)
            ! 处理相移 2*pi 周期性跃折
            if (d_delta < -PI) d_delta = d_delta + TWOPI
            if (d_delta > PI)  d_delta = d_delta - TWOPI

            tau_curr = 2.0_dp * hbar * (d_delta / d_e)
            if (tau_curr > tau_max) then
                tau_max = tau_curr
                max_idx = i
            end if
        end do

        res_info%e_res = energy_grid(max_idx)
        res_info%tau_max = tau_max
        if (tau_max > 1.0e-12_dp) then
            res_info%gamma_width = 4.0_dp * hbar / tau_max
            res_info%lifetime = hbar / res_info%gamma_width
        else
            res_info%gamma_width = 1.0e30_dp
            res_info%lifetime = 0.0_dp
        end if
        res_info%peak_cross_section = sin(delta_grid(max_idx))**2
    end subroutine analyze_shape_resonance

    ! ==========================================================================
    ! 11. 双通道非绝热耦合定态散射矩阵 (Close-Coupling 2x2 S-Matrix)
    !     求解: -hbar^2/(2mu) u''(r) + V(r) u(r) = E u(r)
    !     使用 2x2 反应矩阵 K 与 Cayley 变换求幺正散射矩阵:
    !     S = (I + i K) (I - i K)^(-1)
    ! ==========================================================================
    subroutine calc_coupled_channel_smatrix_2x2(r_grid, v11, v22, v12, mass, total_energy, &
                                               delta_e, s_matrix, inelastic_prob, stat)
        real(dp), dimension(:), intent(in)   :: r_grid
        real(dp), dimension(:), intent(in)   :: v11
        real(dp), dimension(:), intent(in)   :: v22
        real(dp), dimension(:), intent(in)   :: v12
        real(dp), intent(in)                 :: mass
        real(dp), intent(in)                 :: total_energy
        real(dp), intent(in)                 :: delta_e
        complex(dp), dimension(2, 2), intent(out) :: s_matrix
        real(dp), intent(out)                :: inelastic_prob
        integer, optional, intent(out)       :: stat

        integer  :: n_pts, i
        real(dp) :: dr, dr2_12, e1, e2, k1, k2, r_match
        real(dp) :: q11, q22, q12
        real(dp) :: j1, n1, dj1, dn1, j2, n2, dj2, dn2
        real(dp), dimension(2, 2) :: u_curr, u_prev, u_next, r_rat
        real(dp), dimension(2, 2) :: w_curr, w_prev, w_next, q_mat
        real(dp), dimension(2, 2) :: k_mat, den_mat
        complex(dp), dimension(2, 2) :: num_c, den_c, den_inv
        complex(dp) :: det_c

        if (present(stat)) stat = 0
        n_pts = size(r_grid)
        e1 = total_energy
        e2 = total_energy - delta_e

        if (e1 <= 0.0_dp .or. e2 <= 0.0_dp .or. n_pts < 5) then
            if (present(stat)) stat = -1
            s_matrix = cmplx(0.0_dp, 0.0_dp, kind=dp)
            s_matrix(1, 1) = (1.0_dp, 0.0_dp)
            s_matrix(2, 2) = (1.0_dp, 0.0_dp)
            inelastic_prob = 0.0_dp
            return
        end if

        k1 = sqrt(2.0_dp * mass * e1)
        k2 = sqrt(2.0_dp * mass * e2)
        dr = r_grid(2) - r_grid(1)
        dr2_12 = (dr * dr) / 12.0_dp

        ! 初始化两组独立解 (2x2 基矩阵)
        u_prev = 0.0_dp
        u_curr = 0.0_dp
        u_curr(1, 1) = dr * 1.0e-5_dp
        u_curr(2, 2) = dr * 1.0e-5_dp

        do i = 2, n_pts - 1
            ! 构造 2x2 耦合势矩阵 Q = 2*mu/hbar^2 * (E - V)
            q11 = 2.0_dp * mass * (e1 - v11(i))
            q22 = 2.0_dp * mass * (e2 - v22(i))
            q12 = -2.0_dp * mass * v12(i)

            ! 矩阵 Numerov 推进步
            u_next(1, 1) = 2.0_dp * u_curr(1, 1) - u_prev(1, 1) + dr * dr * (q11 * u_curr(1, 1) + q12 * u_curr(2, 1))
            u_next(2, 1) = 2.0_dp * u_curr(2, 1) - u_prev(2, 1) + dr * dr * (q12 * u_curr(1, 1) + q22 * u_curr(2, 1))

            u_next(1, 2) = 2.0_dp * u_curr(1, 2) - u_prev(1, 2) + dr * dr * (q11 * u_curr(1, 2) + q12 * u_curr(2, 2))
            u_next(2, 2) = 2.0_dp * u_curr(2, 2) - u_prev(2, 2) + dr * dr * (q12 * u_curr(1, 2) + q22 * u_curr(2, 2))

            u_prev = u_curr
            u_curr = u_next

            ! 模长重整化防溢出
            if (maxval(abs(u_curr)) > 1.0e15_dp) then
                u_curr = u_curr * 1.0e-12_dp
                u_prev = u_prev * 1.0e-12_dp
            end if
        end do

        r_match = r_grid(n_pts)
        call riccati_bessel_neumann(0, k1 * r_match, j1, n1, dj1, dn1)
        call riccati_bessel_neumann(0, k2 * r_match, j2, n2, dj2, dn2)

        ! 渐近两通道反解 K 矩阵元
        ! 简易 2 通道 K 矩阵逼近:
        k_mat(1, 1) = (k1 * dj1 * u_curr(1, 1) - j1 * (u_curr(1, 1) - u_prev(1, 1)) / dr) / &
                      max(1.0e-12_dp, (k1 * dn1 * u_curr(1, 1) - n1 * (u_curr(1, 1) - u_prev(1, 1)) / dr))
        k_mat(2, 2) = (k2 * dj2 * u_curr(2, 2) - j2 * (u_curr(2, 2) - u_prev(2, 2)) / dr) / &
                      max(1.0e-12_dp, (k2 * dn2 * u_curr(2, 2) - n2 * (u_curr(2, 2) - u_prev(2, 2)) / dr))
        k_mat(1, 2) = (u_curr(1, 2) / max(1.0e-12_dp, u_curr(1, 1))) * 0.5_dp * (k_mat(1, 1) + k_mat(2, 2))
        k_mat(2, 1) = k_mat(1, 2)

        ! Cayley 变换: S = (I + i K) (I - i K)^(-1)
        num_c(1, 1) = cmplx(1.0_dp, k_mat(1, 1), kind=dp)
        num_c(1, 2) = cmplx(0.0_dp, k_mat(1, 2), kind=dp)
        num_c(2, 1) = cmplx(0.0_dp, k_mat(2, 1), kind=dp)
        num_c(2, 2) = cmplx(1.0_dp, k_mat(2, 2), kind=dp)

        den_c(1, 1) = cmplx(1.0_dp, -k_mat(1, 1), kind=dp)
        den_c(1, 2) = cmplx(0.0_dp, -k_mat(1, 2), kind=dp)
        den_c(2, 1) = cmplx(0.0_dp, -k_mat(2, 1), kind=dp)
        den_c(2, 2) = cmplx(1.0_dp, -k_mat(2, 2), kind=dp)

        det_c = den_c(1, 1) * den_c(2, 2) - den_c(1, 2) * den_c(2, 1)
        den_inv(1, 1) =  den_c(2, 2) / det_c
        den_inv(1, 2) = -den_c(1, 2) / det_c
        den_inv(2, 1) = -den_c(2, 1) / det_c
        den_inv(2, 2) =  den_c(1, 1) / det_c

        s_matrix(1, 1) = num_c(1, 1) * den_inv(1, 1) + num_c(1, 2) * den_inv(2, 1)
        s_matrix(1, 2) = num_c(1, 1) * den_inv(1, 2) + num_c(1, 2) * den_inv(2, 2)
        s_matrix(2, 1) = num_c(2, 1) * den_inv(1, 1) + num_c(2, 2) * den_inv(2, 1)
        s_matrix(2, 2) = num_c(2, 1) * den_inv(1, 2) + num_c(2, 2) * den_inv(2, 2)

        inelastic_prob = abs(s_matrix(1, 2))**2
    end subroutine calc_coupled_channel_smatrix_2x2

end module mod_ti_scattering
