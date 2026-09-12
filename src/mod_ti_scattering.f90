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
    use mod_linear_algebra, only: inv_real_matrix, inv_complex_matrix
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

    type, public :: multichannel_result_t
        integer :: n_channels                     ! 总通道数
        integer :: n_open                         ! 开通道数 (E > E_thresh)
        integer :: n_closed                       ! 闭通道数 (E <= E_thresh)
        integer, allocatable :: open_channels(:)   ! 开通道全局索引列表 (1..n_open)
        integer, allocatable :: closed_channels(:) ! 闭通道全局索引列表 (1..n_closed)
        real(dp), allocatable :: k_open(:)         ! 开通道相对动量 k_i = sqrt(2*mu*(E-E_i))/hbar
        real(dp), allocatable :: kappa_closed(:)   ! 闭通道衰减因子 kappa_i = sqrt(2*mu*(E_i-E))/hbar
        real(dp), allocatable :: k_matrix(:, :)    ! 开通道实对称反应矩阵 K_oo (n_open, n_open)
        complex(dp), allocatable :: s_matrix(:, :) ! 开通道严格幺正散射矩阵 S_oo (n_open, n_open)
        complex(dp), allocatable :: t_matrix(:, :) ! 跃迁矩阵 T_oo = S_oo - I
        real(dp), allocatable :: prob_matrix(:, :) ! 态-态跃迁几率 P(i->j) = |S_ij|^2
        real(dp), allocatable :: cross_sections(:, :) ! 态-态部分截面 sigma(i->j)
        real(dp), allocatable :: total_cross_sec(:)   ! 初态总散射截面 sigma_tot(i)
        real(dp) :: eigenphase_sum                 ! 特征相移和 delta_sum (rad)
    end type multichannel_result_t

    ! 多扇区分段径向网格结构 (Segmented / Multi-Sector Radial Grid)
    type, public :: segmented_grid_t
        integer :: n_sectors                        ! 扇区数量
        real(dp), allocatable :: sector_rmin(:)     ! 各扇区起始径向坐标 (n_sectors)
        real(dp), allocatable :: sector_rmax(:)     ! 各扇区终止径向坐标 (n_sectors)
        real(dp), allocatable :: sector_dr(:)       ! 各扇区网格步长 (n_sectors)
        integer, allocatable  :: sector_npts(:)     ! 各扇区内格点数 (n_sectors)
        integer, allocatable  :: sector_offset(:)   ! 各扇区在全局展平网格中的起始索引 (n_sectors)
        integer               :: n_total            ! 全局总格点数
        real(dp), allocatable :: r(:)               ! 展平后的连续单调递增全局径向网格 (1..n_total)
    end type segmented_grid_t

    ! 全同粒子统计与量子对称性参数
    integer, parameter, public :: PARTICLE_DISTINGUISHABLE             = 0
    integer, parameter, public :: PARTICLE_IDENTICAL_BOSON             = 1
    integer, parameter, public :: PARTICLE_IDENTICAL_FERMION_POLARIZED = 2
    integer, parameter, public :: PARTICLE_IDENTICAL_FERMION_UNPOLAR   = 3

    ! 连续态归一化类型常数
    integer, parameter, public :: NORM_ENERGY          = 1 ! delta(E - E') 能量归一化: 渐近振幅 sqrt(2*mu/(pi*hbar^2*k))
    integer, parameter, public :: NORM_MOMENTUM        = 2 ! delta(k - k') 动量归一化: 渐近振幅 sqrt(2/pi)
    integer, parameter, public :: NORM_UNIT_AMPLITUDE  = 3 ! 驻波渐近单位振幅: 渐近振幅 1.0

    ! --------------------------------------------------------------------------
    ! 公共接口导出
    ! --------------------------------------------------------------------------
    public :: riccati_bessel_neumann
    public :: calc_scattering_length_numerov
    public :: calc_scattering_length_logder
    public :: calc_phase_shift_single_l
    public :: calc_scattering_wavefunction_ti
    public :: calc_scattering_wavefunction_1d_cartesian
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
    public :: calc_multichannel_close_coupling_logder
    public :: calc_feshbach_resonance_scan

    ! 多扇区分段网格专用接口导出
    public :: create_segmented_grid
    public :: calc_scattering_length_segmented_numerov
    public :: calc_scattering_wavefunction_segmented_ti
    public :: calc_phase_shift_segmented
    public :: calc_multichannel_close_coupling_segmented_logder

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
    ! 4.1 非含时方法求解定态散射能量本征波函数 u_{l, E}(r)
    !     积分径向薛定谔方程，并使用精确渐近 Riccati 函数与相移进行能量/动量正交归一化
    ! ==========================================================================
    subroutine calc_scattering_wavefunction_ti( &
        r_grid, v_pot, mass, energy, l, norm_type, u_wf, phase_shift, stat)

        real(dp), dimension(:), intent(in)  :: r_grid
        real(dp), dimension(:), intent(in)  :: v_pot
        real(dp), intent(in)                :: mass
        real(dp), intent(in)                :: energy
        integer,  intent(in)                :: l
        integer,  intent(in)                :: norm_type
        real(dp), dimension(:), intent(out) :: u_wf
        real(dp), intent(out)               :: phase_shift
        integer, optional, intent(out)      :: stat

        integer  :: n_pts, i
        real(dp) :: dr, dr2_12, k_wave, r_match
        real(dp) :: q_prev, q_curr, q_next
        real(dp) :: c_prev, c_curr, c_next, d_u, y_logder
        real(dp) :: jl, nl, d_jl, d_nl, num, den
        real(dp) :: asymp_amp, norm_target, scale_factor, target_asymp

        if (present(stat)) stat = 0
        n_pts = size(r_grid)
        u_wf = 0.0_dp
        phase_shift = 0.0_dp

        if (energy <= 0.0_dp .or. n_pts < 5) then
            if (present(stat)) stat = -1
            return
        end if

        k_wave = sqrt(2.0_dp * mass * energy)
        dr = r_grid(2) - r_grid(1)
        dr2_12 = (dr * dr) / 12.0_dp

        ! 1. 从原点正则边界条件出发，Numerov 逐点积分
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

        ! 2. 外边界计算对数导数并匹配 Riccati 函数提取相移 delta_l
        r_match = r_grid(n_pts)
        d_u = (3.0_dp * u_wf(n_pts) - 4.0_dp * u_wf(n_pts - 1) + u_wf(n_pts - 2)) / (2.0_dp * dr)
        y_logder = d_u / u_wf(n_pts)

        call riccati_bessel_neumann(l, k_wave * r_match, jl, nl, d_jl, d_nl)
        num = k_wave * d_jl - y_logder * jl
        den = k_wave * d_nl - y_logder * nl
        phase_shift = atan2(num, den)

        ! 3. 严格连续态物理归一化 (Asymptotic Normalization)
        ! 数值解在渐近区的振幅: A_num = sqrt( u^2 + (u'/k)^2 )
        asymp_amp = sqrt(u_wf(n_pts)**2 + (d_u / k_wave)**2)
        if (asymp_amp < 1.0e-30_dp) asymp_amp = 1.0e-30_dp

        select case (norm_type)
        case (NORM_ENERGY)
            ! delta(E - E') 归一化: 振幅为 sqrt(2*mu / (pi*hbar^2*k))
            norm_target = sqrt(2.0_dp * mass / (PI * k_wave))
        case (NORM_MOMENTUM)
            ! delta(k - k') 归一化: 振幅为 sqrt(2 / pi)
            norm_target = sqrt(2.0_dp / PI)
        case default
            ! 驻波单位振幅: 振幅为 1.0
            norm_target = 1.0_dp
        end select

        ! 相位对齐: 确保正负号与 cos(delta)*jl - sin(delta)*nl 严格同号
        target_asymp = cos(phase_shift) * jl - sin(phase_shift) * nl
        scale_factor = norm_target / asymp_amp
        if (u_wf(n_pts) * target_asymp < 0.0_dp) then
            scale_factor = -scale_factor
        end if

        u_wf = u_wf * scale_factor
    end subroutine calc_scattering_wavefunction_ti

    ! ==========================================================================
    ! 4.2 非含时一维笛卡尔定态散射能量本征波函数 psi_E(x) 求解器
    !     从右边界透射波 psi ~ exp(i*k*x) 逆向 Numerov 积分至左边界，
    !     分解得到入射波振幅 A_inc 与反射波振幅 B_ref，
    !     计算透射几率 T = 1/|A|^2 与反射几率 R = |B|^2/|A|^2，
    !     并进行严格连续能量归一化 psi_E(x) = (psi/A) / sqrt(2*pi) * sqrt(m/(hbar^2*k))
    ! ==========================================================================
    subroutine calc_scattering_wavefunction_1d_cartesian( &
        x_grid, v_pot, mass, energy, norm_type, psi_wf, trans_prob, refl_prob, stat)

        real(dp), dimension(:), intent(in)     :: x_grid
        real(dp), dimension(:), intent(in)     :: v_pot
        real(dp), intent(in)                   :: mass
        real(dp), intent(in)                   :: energy
        integer,  intent(in)                   :: norm_type
        complex(dp), dimension(:), intent(out) :: psi_wf
        real(dp), intent(out)                  :: trans_prob
        real(dp), intent(out)                  :: refl_prob
        integer, optional, intent(out)         :: stat

        integer  :: n_pts, i
        real(dp) :: dx, dx2_12, k_wave
        real(dp), allocatable :: q(:)
        complex(dp) :: c_curr, c_next, c_prev, d_psi_left, a_inc, b_ref
        real(dp) :: norm_factor

        if (present(stat)) stat = 0
        n_pts = size(x_grid)
        psi_wf = (0.0_dp, 0.0_dp)
        trans_prob = 0.0_dp
        refl_prob = 0.0_dp

        if (energy <= 0.0_dp .or. n_pts < 5) then
            if (present(stat)) stat = -1
            return
        end if

        k_wave = sqrt(2.0_dp * mass * energy)
        dx = x_grid(2) - x_grid(1)
        dx2_12 = (dx * dx) / 12.0_dp

        allocate(q(n_pts))
        do i = 1, n_pts
            q(i) = 2.0_dp * mass * (energy - v_pot(i))
        end do

        ! 右边界纯透射波初始条件: psi(x) ~ exp(i * k * x)
        psi_wf(n_pts)     = cmplx(cos(k_wave * x_grid(n_pts)), sin(k_wave * x_grid(n_pts)), kind=dp)
        psi_wf(n_pts - 1) = cmplx(cos(k_wave * x_grid(n_pts - 1)), sin(k_wave * x_grid(n_pts - 1)), kind=dp)

        ! 逆向 Numerov 积分至左边界
        do i = n_pts - 1, 2, -1
            c_curr = 2.0_dp * (1.0_dp - 5.0_dp * dx2_12 * q(i)) * psi_wf(i)
            c_next = (1.0_dp + dx2_12 * q(i + 1)) * psi_wf(i + 1)
            c_prev = cmplx(1.0_dp + dx2_12 * q(i - 1), 0.0_dp, kind=dp)
            psi_wf(i - 1) = (c_curr - c_next) / c_prev

            if (abs(psi_wf(i - 1)) > 1.0e20_dp) then
                psi_wf(i - 1:n_pts) = psi_wf(i - 1:n_pts) * 1.0e-15_dp
            end if
        end do

        ! 左边界数值导数
        d_psi_left = (-3.0_dp * psi_wf(1) + 4.0_dp * psi_wf(2) - psi_wf(3)) / (2.0_dp * dx)

        ! 分解入射与反射波: psi(x_L) = A_inc * exp(i*k*x_L) + B_ref * exp(-i*k*x_L)
        a_inc = 0.5_dp * (psi_wf(1) - (0.0_dp, 1.0_dp) * d_psi_left / k_wave) * &
                cmplx(cos(k_wave * x_grid(1)), -sin(k_wave * x_grid(1)), kind=dp)
        b_ref = 0.5_dp * (psi_wf(1) + (0.0_dp, 1.0_dp) * d_psi_left / k_wave) * &
                cmplx(cos(k_wave * x_grid(1)),  sin(k_wave * x_grid(1)), kind=dp)

        if (abs(a_inc) < 1.0e-30_dp) then
            if (present(stat)) stat = -2
            deallocate(q)
            return
        end if

        trans_prob = 1.0_dp / (abs(a_inc)**2)
        refl_prob  = (abs(b_ref)**2) / (abs(a_inc)**2)

        ! 归一化标定
        select case (norm_type)
        case (NORM_ENERGY)
            ! delta(E - E') 归一化: 入射振幅为 1/sqrt(2*pi) * sqrt(m / (hbar^2 * k))
            norm_factor = (1.0_dp / sqrt(TWOPI)) * sqrt(mass / k_wave)
        case (NORM_MOMENTUM)
            ! delta(k - k') 归一化: 入射振幅为 1/sqrt(2*pi)
            norm_factor = 1.0_dp / sqrt(TWOPI)
        case default
            ! 驻波或单位入射振幅: A_inc = 1.0
            norm_factor = 1.0_dp
        end select

        psi_wf = (psi_wf / a_inc) * norm_factor
        deallocate(q)
    end subroutine calc_scattering_wavefunction_1d_cartesian

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

    ! ==========================================================================
    ! 通用 N 通道定态密耦 Johnson 矩阵对数导数法 (Matrix Log-Derivative Close-Coupling)
    ! 支持任意通道数、开通道/闭通道混合边界以及 Feshbach 共振
    ! 参考文献: B. R. Johnson, J. Comput. Phys. 13, 445 (1973);
    !          D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986).
    ! ==========================================================================
    subroutine calc_multichannel_close_coupling_logder( &
        r_grid, v_mat, mass, total_energy, thresholds, l_channels, &
        res, stat)

        real(dp), intent(in) :: r_grid(:)                    ! 径向格点 (1..n_pts)
        real(dp), intent(in) :: v_mat(:, :, :)               ! 耦合势矩阵 V_ij(r), 形状 (n_chan, n_chan, n_pts)
        real(dp), intent(in) :: mass                         ! 碰撞体系折合质量 mu
        real(dp), intent(in) :: total_energy                 ! 总碰撞能量 E
        real(dp), intent(in) :: thresholds(:)                ! 各通道渐近能级阈值 E_i^thresh (1..n_chan)
        integer,  intent(in) :: l_channels(:)                ! 各通道轨道角动量 l_i (1..n_chan)
        type(multichannel_result_t), intent(out) :: res      ! 输出多通道定态散射物理结果
        integer,  intent(out), optional :: stat              ! 状态码 (0: 正常)

        integer :: n_chan, n_pts, i, j, step, stat_inv
        real(dp) :: dr, dr2_12, r_curr, r_match
        real(dp), allocatable :: w_mat(:, :), q_mat(:, :), q_inv(:, :), m_mat(:, :)
        real(dp), allocatable :: r_curr_mat(:, :), r_next_mat(:, :), r_inv(:, :)
        real(dp), allocatable :: q_prev(:, :), q_curr(:, :)
        real(dp), allocatable :: p1(:, :), p2(:, :), temp_mat(:, :)
        real(dp), allocatable :: y_mat(:, :)
        real(dp), allocatable :: y_oo(:, :), y_oc(:, :), y_co(:, :), y_cc(:, :)
        real(dp), allocatable :: a_cc(:, :), a_cc_inv(:, :)
        real(dp), allocatable :: y_eff(:, :)
        real(dp), allocatable :: j_mat(:, :), n_mat(:, :), dj_mat(:, :), dn_mat(:, :)
        real(dp), allocatable :: mj_mat(:, :), mn_mat(:, :), mn_inv(:, :)
        complex(dp), allocatable :: eye_c(:, :), ik_mat(:, :), den_c(:, :), den_inv(:, :)
        complex(dp), allocatable :: s_mat(:, :)
        real(dp) :: jl, nl, djl, dnl, k_i, e_kin
        complex(dp) :: det_s

        if (present(stat)) stat = 0
        n_pts = size(r_grid)
        n_chan = size(thresholds)

        res%n_channels = n_chan
        res%n_open = 0
        res%n_closed = 0

        if (n_pts < 5 .or. n_chan < 1) then
            if (present(stat)) stat = -1
            return
        end if

        dr = r_grid(2) - r_grid(1)
        dr2_12 = (dr * dr) / 12.0_dp

        ! 1. 统计并分类开通道与闭通道
        allocate(res%open_channels(n_chan))
        allocate(res%closed_channels(n_chan))

        do i = 1, n_chan
            e_kin = total_energy - thresholds(i)
            if (e_kin > 1.0e-13_dp) then
                res%n_open = res%n_open + 1
                res%open_channels(res%n_open) = i
            else
                res%n_closed = res%n_closed + 1
                res%closed_channels(res%n_closed) = i
            end if
        end do

        if (res%n_open == 0) then
            ! 所有通道均为闭通道，体系处于全禁区束缚态区域
            if (present(stat)) stat = 1
            return
        end if

        allocate(res%k_open(res%n_open))
        do i = 1, res%n_open
            res%k_open(i) = sqrt(2.0_dp * mass * (total_energy - thresholds(res%open_channels(i))))
        end do

        if (res%n_closed > 0) then
            allocate(res%kappa_closed(res%n_closed))
            do i = 1, res%n_closed
                res%kappa_closed(i) = sqrt(2.0_dp * mass * max(0.0_dp, thresholds(res%closed_channels(i)) - total_energy))
            end do
        end if

        ! 2. 初始化 Johnson 矩阵对数导数 / 比值矩阵推进器
        allocate(w_mat(n_chan, n_chan))
        allocate(q_mat(n_chan, n_chan))
        allocate(q_inv(n_chan, n_chan))
        allocate(m_mat(n_chan, n_chan))
        allocate(r_curr_mat(n_chan, n_chan))
        allocate(r_next_mat(n_chan, n_chan))
        allocate(r_inv(n_chan, n_chan))
        allocate(q_prev(n_chan, n_chan))
        allocate(q_curr(n_chan, n_chan))
        allocate(p1(n_chan, n_chan))
        allocate(p2(n_chan, n_chan))
        allocate(temp_mat(n_chan, n_chan))
        allocate(y_mat(n_chan, n_chan))

        ! 计算第一格点的 W(r_1) 与 Q(r_1)
        r_curr = r_grid(1)
        do i = 1, n_chan
            do j = 1, n_chan
                w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, 1)
                if (i == j) then
                    w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                  real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                end if
            end do
        end do

        q_mat = -dr2_12 * w_mat
        do i = 1, n_chan
            q_mat(i, i) = q_mat(i, i) + 1.0_dp
        end do

        call inv_real_matrix(n_chan, q_mat, q_inv, stat_inv)
        if (stat_inv /= 0) then
            if (present(stat)) stat = -2
            return
        end if

        ! M_1 = 12 * Q_1^{-1} - 10 * I
        m_mat = 12.0_dp * q_inv
        do i = 1, n_chan
            m_mat(i, i) = m_mat(i, i) - 10.0_dp
        end do

        ! 第一步 R_1^{-1} = 0 (原点波函数 Psi(r_0) = 0)
        ! 从而 R_2 = M_1 - 0 = M_1
        r_curr_mat = m_mat
        q_prev = q_mat

        ! 3. 循环递推推进到边界 r_N
        do step = 2, n_pts - 1
            r_curr = r_grid(step)

            ! 构建 W(r_step)
            do i = 1, n_chan
                do j = 1, n_chan
                    w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, step)
                    if (i == j) then
                        w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                      real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                    end if
                end do
            end do

            q_mat = -dr2_12 * w_mat
            do i = 1, n_chan
                q_mat(i, i) = q_mat(i, i) + 1.0_dp
            end do

            call inv_real_matrix(n_chan, q_mat, q_inv, stat_inv)
            if (stat_inv /= 0) then
                if (present(stat)) stat = -3
                return
            end if

            m_mat = 12.0_dp * q_inv
            do i = 1, n_chan
                m_mat(i, i) = m_mat(i, i) - 10.0_dp
            end do

            ! 求 R_step^{-1}
            call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)
            if (stat_inv /= 0) then
                ! 遇到极点微小微扰正则化
                do i = 1, n_chan
                    r_curr_mat(i, i) = r_curr_mat(i, i) + 1.0e-14_dp
                end do
                call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)
            end if

            ! Johnson 递推: R_{step+1} = M_step - R_step^{-1}
            r_next_mat = m_mat - r_inv
            ! 严格保持对称性
            r_next_mat = 0.5_dp * (r_next_mat + transpose(r_next_mat))

            if (step == n_pts - 2) then
                q_curr = q_mat
            end if

            r_curr_mat = r_next_mat
        end do

        ! 计算外边界 r_N 处的对数导数矩阵 Y(r_N)
        ! 计算 P1 = Q_{N-1}^{-1} * R_N^{-1} * Q_N
        call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)

        ! 计算格点 N 处的 Q_N
        r_curr = r_grid(n_pts)
        do i = 1, n_chan
            do j = 1, n_chan
                w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, n_pts)
                if (i == j) then
                    w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                  real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                end if
            end do
        end do
        q_mat = -dr2_12 * w_mat
        do i = 1, n_chan
            q_mat(i, i) = q_mat(i, i) + 1.0_dp
        end do

        call inv_real_matrix(n_chan, q_curr, q_inv, stat_inv)
        temp_mat = matmul(r_inv, q_mat)
        p1 = matmul(q_inv, temp_mat)

        ! P2 = P1 * P1 (利用两级向后逼近)
        p2 = matmul(p1, p1)

        ! Y = (3*I - 4*P1 + P2) / (2*dr)
        y_mat = p2 - 4.0_dp * p1
        do i = 1, n_chan
            y_mat(i, i) = y_mat(i, i) + 3.0_dp
        end do
        y_mat = y_mat / (2.0_dp * dr)
        y_mat = 0.5_dp * (y_mat + transpose(y_mat))

        ! 4. 开通道/闭通道 Schur 补变换 (Feshbach 投影)
        allocate(y_eff(res%n_open, res%n_open))

        if (res%n_closed == 0) then
            y_eff = y_mat
        else
            allocate(y_oo(res%n_open, res%n_open))
            allocate(y_oc(res%n_open, res%n_closed))
            allocate(y_co(res%n_closed, res%n_open))
            allocate(y_cc(res%n_closed, res%n_closed))
            allocate(a_cc(res%n_closed, res%n_closed))
            allocate(a_cc_inv(res%n_closed, res%n_closed))

            do i = 1, res%n_open
                do j = 1, res%n_open
                    y_oo(j, i) = y_mat(res%open_channels(j), res%open_channels(i))
                end do
            end do

            do i = 1, res%n_closed
                do j = 1, res%n_open
                    y_oc(j, i) = y_mat(res%open_channels(j), res%closed_channels(i))
                    y_co(i, j) = y_mat(res%closed_channels(i), res%open_channels(j))
                end do
            end do

            do i = 1, res%n_closed
                do j = 1, res%n_closed
                    y_cc(j, i) = y_mat(res%closed_channels(j), res%closed_channels(i))
                end do
            end do

            ! A_cc = Y_cc + diag(kappa_closed)
            a_cc = y_cc
            do i = 1, res%n_closed
                a_cc(i, i) = a_cc(i, i) + res%kappa_closed(i)
            end do

            call inv_real_matrix(res%n_closed, a_cc, a_cc_inv, stat_inv)
            if (stat_inv /= 0) then
                ! 遇到 Feshbach 奇异共振点微扰正则化
                do i = 1, res%n_closed
                    a_cc(i, i) = a_cc(i, i) + 1.0e-10_dp
                end do
                call inv_real_matrix(res%n_closed, a_cc, a_cc_inv, stat_inv)
            end if

            ! Y_eff = Y_oo - Y_oc * A_cc^{-1} * Y_co
            y_eff = y_oo - matmul(y_oc, matmul(a_cc_inv, y_co))
            y_eff = 0.5_dp * (y_eff + transpose(y_eff))

            deallocate(y_oo, y_oc, y_co, y_cc, a_cc, a_cc_inv)
        end if

        ! 5. 渐近 Riccati 函数边界匹配提取反应矩阵 K_oo
        r_match = r_grid(n_pts)
        allocate(j_mat(res%n_open, res%n_open))
        allocate(n_mat(res%n_open, res%n_open))
        allocate(dj_mat(res%n_open, res%n_open))
        allocate(dn_mat(res%n_open, res%n_open))
        allocate(mj_mat(res%n_open, res%n_open))
        allocate(mn_mat(res%n_open, res%n_open))
        allocate(mn_inv(res%n_open, res%n_open))
        allocate(res%k_matrix(res%n_open, res%n_open))

        j_mat = 0.0_dp; n_mat = 0.0_dp; dj_mat = 0.0_dp; dn_mat = 0.0_dp

        do i = 1, res%n_open
            k_i = res%k_open(i)
            call riccati_bessel_neumann(l_channels(res%open_channels(i)), k_i * r_match, jl, nl, djl, dnl)
            ! 通量归一化对角元
            j_mat(i, i)  = jl / sqrt(k_i)
            n_mat(i, i)  = nl / sqrt(k_i)
            dj_mat(i, i) = djl * sqrt(k_i)
            dn_mat(i, i) = dnl * sqrt(k_i)
        end do

        ! M_J = J' - Y_eff * J
        ! M_N = N' - Y_eff * N
        mj_mat = dj_mat - matmul(y_eff, j_mat)
        mn_mat = dn_mat - matmul(y_eff, n_mat)

        call inv_real_matrix(res%n_open, mn_mat, mn_inv, stat_inv)
        if (stat_inv /= 0) then
            do i = 1, res%n_open
                mn_mat(i, i) = mn_mat(i, i) + 1.0e-12_dp
            end do
            call inv_real_matrix(res%n_open, mn_mat, mn_inv, stat_inv)
        end if

        ! K_oo = (M_N)^{-1} * M_J
        res%k_matrix = matmul(mn_inv, mj_mat)
        res%k_matrix = 0.5_dp * (res%k_matrix + transpose(res%k_matrix))

        ! 6. Cayley 变换构建严格么正散射矩阵 S_oo = (I + i*K) * (I - i*K)^{-1}
        allocate(eye_c(res%n_open, res%n_open))
        allocate(ik_mat(res%n_open, res%n_open))
        allocate(den_c(res%n_open, res%n_open))
        allocate(den_inv(res%n_open, res%n_open))
        allocate(res%s_matrix(res%n_open, res%n_open))
        allocate(res%t_matrix(res%n_open, res%n_open))

        eye_c = (0.0_dp, 0.0_dp)
        do i = 1, res%n_open
            eye_c(i, i) = (1.0_dp, 0.0_dp)
        end do

        ik_mat = cmplx(0.0_dp, res%k_matrix, kind=dp)
        den_c = eye_c - ik_mat

        call inv_complex_matrix(res%n_open, den_c, den_inv, stat_inv)
        res%s_matrix = matmul(eye_c + ik_mat, den_inv)
        res%t_matrix = res%s_matrix - eye_c

        ! 7. 导出态-态跃迁几率与散射截面
        allocate(res%prob_matrix(res%n_open, res%n_open))
        allocate(res%cross_sections(res%n_open, res%n_open))
        allocate(res%total_cross_sec(res%n_open))

        do i = 1, res%n_open
            res%total_cross_sec(i) = 0.0_dp
            k_i = res%k_open(i)
            do j = 1, res%n_open
                res%prob_matrix(j, i) = abs(res%s_matrix(j, i))**2
                res%cross_sections(j, i) = (PI / (k_i * k_i)) * &
                    real(2 * l_channels(res%open_channels(i)) + 1, dp) * &
                    abs(res%t_matrix(j, i))**2
                res%total_cross_sec(i) = res%total_cross_sec(i) + res%cross_sections(j, i)
            end do
        end do

        ! 8. 计算特征相移和 delta_sum = 0.5 * arg(det(S))
        det_s = (1.0_dp, 0.0_dp)
        if (res%n_open == 1) then
            det_s = res%s_matrix(1, 1)
        else if (res%n_open == 2) then
            det_s = res%s_matrix(1, 1) * res%s_matrix(2, 2) - res%s_matrix(1, 2) * res%s_matrix(2, 1)
        else
            allocate(s_mat(res%n_open, res%n_open))
            s_mat = res%s_matrix
            do i = 1, res%n_open
                det_s = det_s * s_mat(i, i)
            end do
            deallocate(s_mat)
        end if
        res%eigenphase_sum = 0.5_dp * atan2(aimag(det_s), real(det_s, dp))

        ! 清理动态内存
        deallocate(w_mat, q_mat, q_inv, m_mat, r_curr_mat, r_next_mat, r_inv)
        deallocate(q_prev, q_curr, p1, p2, temp_mat, y_mat, y_eff)
        deallocate(j_mat, n_mat, dj_mat, dn_mat, mj_mat, mn_mat, mn_inv)
        deallocate(eye_c, ik_mat, den_c, den_inv)
    end subroutine calc_multichannel_close_coupling_logder

    ! ==========================================================================
    ! Feshbach 共振能谱能量扫描: 提取开通道散射长度 a_s(E) 与特征相移跃升
    ! ==========================================================================
    subroutine calc_feshbach_resonance_scan( &
        r_grid, v_mat, mass, energy_grid, n_energies, thresholds, l_channels, &
        s_wave_length, eigenphase_sums, stat)

        real(dp), intent(in) :: r_grid(:)
        real(dp), intent(in) :: v_mat(:, :, :)
        real(dp), intent(in) :: mass
        real(dp), intent(in) :: energy_grid(:)
        integer,  intent(in) :: n_energies
        real(dp), intent(in) :: thresholds(:)
        integer,  intent(in) :: l_channels(:)
        real(dp), intent(out) :: s_wave_length(n_energies)
        real(dp), intent(out) :: eigenphase_sums(n_energies)
        integer,  intent(out), optional :: stat

        type(multichannel_result_t) :: res
        integer :: ie, s
        real(dp) :: k1

        if (present(stat)) stat = 0
        s_wave_length = 0.0_dp
        eigenphase_sums = 0.0_dp

        do ie = 1, n_energies
            call calc_multichannel_close_coupling_logder( &
                r_grid, v_mat, mass, energy_grid(ie), thresholds, l_channels, res, s)
            if (s == 0 .and. res%n_open >= 1) then
                k1 = res%k_open(1)
                s_wave_length(ie) = -res%k_matrix(1, 1) / max(1.0e-12_dp, k1)
                eigenphase_sums(ie) = res%eigenphase_sum
            end if
        end do
    end subroutine calc_feshbach_resonance_scan

    ! ==========================================================================
    ! 多扇区分段径向网格构造器 (Multi-Sector Segmented Radial Grid Generator)
    ! ==========================================================================
    subroutine create_segmented_grid(r_start, r_bounds, dr_steps, grid, stat)
        real(dp), intent(in)                :: r_start
        real(dp), dimension(:), intent(in)  :: r_bounds
        real(dp), dimension(:), intent(in)  :: dr_steps
        type(segmented_grid_t), intent(out) :: grid
        integer, optional, intent(out)      :: stat

        integer  :: n_sec, k, n_int, curr_idx, i
        real(dp) :: r_curr_start, r_curr_end, dr_target, dr_actual

        if (present(stat)) stat = 0
        n_sec = size(r_bounds)
        if (n_sec < 1 .or. size(dr_steps) /= n_sec) then
            if (present(stat)) stat = -1
            return
        end if

        grid%n_sectors = n_sec
        allocate(grid%sector_rmin(n_sec))
        allocate(grid%sector_rmax(n_sec))
        allocate(grid%sector_dr(n_sec))
        allocate(grid%sector_npts(n_sec))
        allocate(grid%sector_offset(n_sec))

        curr_idx = 1
        do k = 1, n_sec
            if (k == 1) then
                r_curr_start = r_start
            else
                r_curr_start = r_bounds(k - 1)
            end if
            r_curr_end = r_bounds(k)
            if (r_curr_end <= r_curr_start .or. dr_steps(k) <= 0.0_dp) then
                if (present(stat)) stat = -2
                return
            end if

            dr_target = dr_steps(k)
            n_int = max(4, nint((r_curr_end - r_curr_start) / dr_target))
            dr_actual = (r_curr_end - r_curr_start) / real(n_int, dp)

            grid%sector_rmin(k) = r_curr_start
            grid%sector_rmax(k) = r_curr_end
            grid%sector_dr(k)   = dr_actual
            grid%sector_npts(k) = n_int + 1
            grid%sector_offset(k) = curr_idx

            curr_idx = curr_idx + n_int
        end do

        grid%n_total = curr_idx
        allocate(grid%r(grid%n_total))

        ! 填充全局展平单调径向网格点
        do k = 1, n_sec
            curr_idx = grid%sector_offset(k)
            do i = 1, grid%sector_npts(k)
                grid%r(curr_idx + i - 1) = grid%sector_rmin(k) + real(i - 1, dp) * grid%sector_dr(k)
            end do
        end do
    end subroutine create_segmented_grid

    ! ==========================================================================
    ! 分段扇区网格零能散射长度求解器 (Segmented Grid Numerov Scattering Length)
    ! 跨扇区采用 4 阶 Taylor 导数桥接，保持波函数与能量本征导数光滑无缝连续
    ! ==========================================================================
    subroutine calc_scattering_length_segmented_numerov(grid, v_pot, mass, a_s, u_zero, stat)
        type(segmented_grid_t), intent(in)  :: grid
        real(dp), dimension(:), intent(in)  :: v_pot
        real(dp), intent(in)                :: mass
        real(dp), intent(out)               :: a_s
        real(dp), dimension(:), optional, intent(out) :: u_zero
        integer, optional, intent(out)      :: stat

        integer  :: k, i, j_start, j_end, j_bound, n_total
        real(dp) :: h, h2_12, h_prev
        real(dp) :: q_prev, q_curr, q_next
        real(dp) :: c_prev, c_curr, c_next
        real(dp) :: d_u, d2_u, d3_u, d4_u
        real(dp) :: q_j, q_jm1, q_jm2, d_q, d2_q
        real(dp), allocatable :: u_wf(:)

        if (present(stat)) stat = 0
        n_total = grid%n_total
        a_s = 0.0_dp

        if (n_total < 5 .or. size(v_pot) < n_total) then
            if (present(stat)) stat = -1
            return
        end if

        allocate(u_wf(n_total))
        u_wf = 0.0_dp

        ! 扇区 1 启动
        h = grid%sector_dr(1)
        h2_12 = (h * h) / 12.0_dp
        j_start = grid%sector_offset(1)
        j_end   = j_start + grid%sector_npts(1) - 1

        u_wf(1) = 0.0_dp
        u_wf(2) = 1.0e-7_dp

        do i = j_start + 1, j_end - 1
            q_prev = -2.0_dp * mass * v_pot(i - 1)
            q_curr = -2.0_dp * mass * v_pot(i)
            q_next = -2.0_dp * mass * v_pot(i + 1)

            c_prev = 1.0_dp + h2_12 * q_prev
            c_curr = 2.0_dp * (1.0_dp - 5.0_dp * h2_12 * q_curr)
            c_next = 1.0_dp + h2_12 * q_next

            u_wf(i + 1) = (c_curr * u_wf(i) - c_prev * u_wf(i - 1)) / c_next
            if (abs(u_wf(i + 1)) > 1.0e20_dp) then
                u_wf(1:i + 1) = u_wf(1:i + 1) * 1.0e-15_dp
            end if
        end do

        ! 扇区 2 到 n_sectors 跨扇区无缝推进
        do k = 2, grid%n_sectors
            j_bound = grid%sector_offset(k)
            h_prev  = grid%sector_dr(k - 1)
            h       = grid%sector_dr(k)
            h2_12   = (h * h) / 12.0_dp
            j_end   = j_bound + grid%sector_npts(k) - 1

            ! 计算交界面导数 (以扇区 k-1 的步长 h_prev 做高阶向后差分)
            d_u = (11.0_dp * u_wf(j_bound) - 18.0_dp * u_wf(j_bound - 1) + &
                    9.0_dp * u_wf(j_bound - 2) - 2.0_dp * u_wf(j_bound - 3)) / (6.0_dp * h_prev)

            q_j   = -2.0_dp * mass * v_pot(j_bound)
            q_jm1 = -2.0_dp * mass * v_pot(j_bound - 1)
            q_jm2 = -2.0_dp * mass * v_pot(j_bound - 2)

            d_q  = (3.0_dp * q_j - 4.0_dp * q_jm1 + q_jm2) / (2.0_dp * h_prev)
            d2_q = (q_j - 2.0_dp * q_jm1 + q_jm2) / (h_prev * h_prev)

            d2_u = -q_j * u_wf(j_bound)
            d3_u = -d_q * u_wf(j_bound) - q_j * d_u
            d4_u = (-d2_q + q_j * q_j) * u_wf(j_bound) - 2.0_dp * d_q * d_u

            ! 4 阶 Taylor 桥接生成新扇区第二点
            u_wf(j_bound + 1) = u_wf(j_bound) + h * d_u + &
                                0.5_dp * (h * h) * d2_u + &
                                (h**3 / 6.0_dp) * d3_u + &
                                (h**4 / 24.0_dp) * d4_u

            ! 新扇区内部等步长 Numerov 推进
            do i = j_bound + 1, j_end - 1
                q_prev = -2.0_dp * mass * v_pot(i - 1)
                q_curr = -2.0_dp * mass * v_pot(i)
                q_next = -2.0_dp * mass * v_pot(i + 1)

                c_prev = 1.0_dp + h2_12 * q_prev
                c_curr = 2.0_dp * (1.0_dp - 5.0_dp * h2_12 * q_curr)
                c_next = 1.0_dp + h2_12 * q_next

                u_wf(i + 1) = (c_curr * u_wf(i) - c_prev * u_wf(i - 1)) / c_next
                if (abs(u_wf(i + 1)) > 1.0e20_dp) then
                    u_wf(1:i + 1) = u_wf(1:i + 1) * 1.0e-15_dp
                end if
            end do
        end do

        ! 外边界计算散射长度 a_s
        h = grid%sector_dr(grid%n_sectors)
        d_u = (3.0_dp * u_wf(n_total) - 4.0_dp * u_wf(n_total - 1) + u_wf(n_total - 2)) / (2.0_dp * h)

        if (abs(d_u) > 1.0e-14_dp) then
            a_s = grid%r(n_total) - u_wf(n_total) / d_u
        else
            a_s = 1.0e30_dp
        end if

        if (present(u_zero)) then
            if (size(u_zero) == n_total) then
                u_zero = u_wf
            end if
        end if

        deallocate(u_wf)
    end subroutine calc_scattering_length_segmented_numerov

    ! ==========================================================================
    ! 分段扇区网格定态散射波函数 u_{l, E}(r) 与相移求解器
    ! ==========================================================================
    subroutine calc_scattering_wavefunction_segmented_ti( &
        grid, v_pot, mass, energy, l, norm_type, u_wf, phase_shift, stat)

        type(segmented_grid_t), intent(in)  :: grid
        real(dp), dimension(:), intent(in)  :: v_pot
        real(dp), intent(in)                :: mass
        real(dp), intent(in)                :: energy
        integer,  intent(in)                :: l
        integer,  intent(in)                :: norm_type
        real(dp), dimension(:), intent(out) :: u_wf
        real(dp), intent(out)               :: phase_shift
        integer, optional, intent(out)      :: stat

        integer  :: k, i, j_start, j_end, j_bound, n_total
        real(dp) :: h, h2_12, h_prev, k_wave, r_match
        real(dp) :: q_prev, q_curr, q_next
        real(dp) :: c_prev, c_curr, c_next
        real(dp) :: d_u, d2_u, d3_u, d4_u
        real(dp) :: q_j, q_jm1, q_jm2, d_q, d2_q
        real(dp) :: jl, nl, d_jl, d_nl, num, den, y_logder
        real(dp) :: asymp_amp, norm_target, target_asymp, scale_factor

        if (present(stat)) stat = 0
        n_total = grid%n_total
        u_wf = 0.0_dp
        phase_shift = 0.0_dp

        if (energy <= 0.0_dp .or. n_total < 5 .or. size(v_pot) < n_total) then
            if (present(stat)) stat = -1
            return
        end if

        k_wave = sqrt(2.0_dp * mass * energy)

        ! 扇区 1 启动
        h = grid%sector_dr(1)
        h2_12 = (h * h) / 12.0_dp
        j_start = grid%sector_offset(1)
        j_end   = j_start + grid%sector_npts(1) - 1

        u_wf(1) = 0.0_dp
        u_wf(2) = (h)**(l + 1) * 1.0e-5_dp

        do i = j_start + 1, j_end - 1
            q_prev = 2.0_dp * mass * (energy - v_pot(i - 1)) - real(l * (l + 1), dp) / (grid%r(i - 1)**2)
            q_curr = 2.0_dp * mass * (energy - v_pot(i))     - real(l * (l + 1), dp) / (grid%r(i)**2)
            q_next = 2.0_dp * mass * (energy - v_pot(i + 1)) - real(l * (l + 1), dp) / (grid%r(i + 1)**2)

            c_prev = 1.0_dp + h2_12 * q_prev
            c_curr = 2.0_dp * (1.0_dp - 5.0_dp * h2_12 * q_curr)
            c_next = 1.0_dp + h2_12 * q_next

            u_wf(i + 1) = (c_curr * u_wf(i) - c_prev * u_wf(i - 1)) / c_next
            if (abs(u_wf(i + 1)) > 1.0e20_dp) then
                u_wf(1:i + 1) = u_wf(1:i + 1) * 1.0e-15_dp
            end if
        end do

        ! 扇区 2 到 n_sectors
        do k = 2, grid%n_sectors
            j_bound = grid%sector_offset(k)
            h_prev  = grid%sector_dr(k - 1)
            h       = grid%sector_dr(k)
            h2_12   = (h * h) / 12.0_dp
            j_end   = j_bound + grid%sector_npts(k) - 1

            d_u = (11.0_dp * u_wf(j_bound) - 18.0_dp * u_wf(j_bound - 1) + &
                    9.0_dp * u_wf(j_bound - 2) - 2.0_dp * u_wf(j_bound - 3)) / (6.0_dp * h_prev)

            q_j   = 2.0_dp * mass * (energy - v_pot(j_bound))     - real(l * (l + 1), dp) / (grid%r(j_bound)**2)
            q_jm1 = 2.0_dp * mass * (energy - v_pot(j_bound - 1)) - real(l * (l + 1), dp) / (grid%r(j_bound - 1)**2)
            q_jm2 = 2.0_dp * mass * (energy - v_pot(j_bound - 2)) - real(l * (l + 1), dp) / (grid%r(j_bound - 2)**2)

            d_q  = (3.0_dp * q_j - 4.0_dp * q_jm1 + q_jm2) / (2.0_dp * h_prev)
            d2_q = (q_j - 2.0_dp * q_jm1 + q_jm2) / (h_prev * h_prev)

            d2_u = -q_j * u_wf(j_bound)
            d3_u = -d_q * u_wf(j_bound) - q_j * d_u
            d4_u = (-d2_q + q_j * q_j) * u_wf(j_bound) - 2.0_dp * d_q * d_u

            u_wf(j_bound + 1) = u_wf(j_bound) + h * d_u + &
                                0.5_dp * (h * h) * d2_u + &
                                (h**3 / 6.0_dp) * d3_u + &
                                (h**4 / 24.0_dp) * d4_u

            do i = j_bound + 1, j_end - 1
                q_prev = 2.0_dp * mass * (energy - v_pot(i - 1)) - real(l * (l + 1), dp) / (grid%r(i - 1)**2)
                q_curr = 2.0_dp * mass * (energy - v_pot(i))     - real(l * (l + 1), dp) / (grid%r(i)**2)
                q_next = 2.0_dp * mass * (energy - v_pot(i + 1)) - real(l * (l + 1), dp) / (grid%r(i + 1)**2)

                c_prev = 1.0_dp + h2_12 * q_prev
                c_curr = 2.0_dp * (1.0_dp - 5.0_dp * h2_12 * q_curr)
                c_next = 1.0_dp + h2_12 * q_next

                u_wf(i + 1) = (c_curr * u_wf(i) - c_prev * u_wf(i - 1)) / c_next
                if (abs(u_wf(i + 1)) > 1.0e20_dp) then
                    u_wf(1:i + 1) = u_wf(1:i + 1) * 1.0e-15_dp
                end if
            end do
        end do

        ! 外边界计算对数导数与相移
        h = grid%sector_dr(grid%n_sectors)
        r_match = grid%r(n_total)
        d_u = (3.0_dp * u_wf(n_total) - 4.0_dp * u_wf(n_total - 1) + u_wf(n_total - 2)) / (2.0_dp * h)
        y_logder = d_u / u_wf(n_total)

        call riccati_bessel_neumann(l, k_wave * r_match, jl, nl, d_jl, d_nl)
        num = k_wave * d_jl - y_logder * jl
        den = k_wave * d_nl - y_logder * nl
        phase_shift = atan2(num, den)

        ! 渐近物理归一化
        asymp_amp = sqrt(u_wf(n_total)**2 + (d_u / k_wave)**2)
        if (asymp_amp < 1.0e-30_dp) asymp_amp = 1.0e-30_dp

        select case (norm_type)
        case (NORM_ENERGY)
            norm_target = sqrt(2.0_dp * mass / (PI * k_wave))
        case (NORM_MOMENTUM)
            norm_target = sqrt(2.0_dp / PI)
        case default
            norm_target = 1.0_dp
        end select

        target_asymp = cos(phase_shift) * jl - sin(phase_shift) * nl
        scale_factor = norm_target / asymp_amp
        if (u_wf(n_total) * target_asymp < 0.0_dp) then
            scale_factor = -scale_factor
        end if

        u_wf = u_wf * scale_factor
    end subroutine calc_scattering_wavefunction_segmented_ti

    ! ==========================================================================
    ! 分段扇区网格快速计算分波相移与 S-矩阵元
    ! ==========================================================================
    subroutine calc_phase_shift_segmented( &
        grid, v_pot, mass, energy, l, phase_shift, k_mat, s_mat, t_mat, cross_sec, stat)

        type(segmented_grid_t), intent(in) :: grid
        real(dp), dimension(:), intent(in) :: v_pot
        real(dp), intent(in)               :: mass
        real(dp), intent(in)               :: energy
        integer,  intent(in)               :: l
        real(dp), intent(out)              :: phase_shift
        real(dp), intent(out)              :: k_mat
        complex(dp), intent(out)           :: s_mat
        complex(dp), intent(out)           :: t_mat
        real(dp), intent(out)              :: cross_sec
        integer, optional, intent(out)     :: stat

        real(dp), allocatable :: u_wf(:)
        real(dp) :: k_wave

        if (present(stat)) stat = 0
        allocate(u_wf(grid%n_total))

        call calc_scattering_wavefunction_segmented_ti( &
            grid, v_pot, mass, energy, l, NORM_UNIT_AMPLITUDE, u_wf, phase_shift, stat)

        k_wave = sqrt(2.0_dp * mass * energy)
        k_mat = tan(phase_shift)
        s_mat = cmplx(cos(2.0_dp * phase_shift), sin(2.0_dp * phase_shift), kind=dp)
        t_mat = s_mat - (1.0_dp, 0.0_dp)

        if (k_wave > 1.0e-14_dp) then
            cross_sec = (4.0_dp * PI / (k_wave * k_wave)) * real(2 * l + 1, dp) * (sin(phase_shift)**2)
        else
            cross_sec = 0.0_dp
        end if

        deallocate(u_wf)
    end subroutine calc_phase_shift_segmented

    ! ==========================================================================
    ! 分段扇区网格多通道密耦 Johnson Log-Derivative 求解器
    ! ==========================================================================
    subroutine calc_multichannel_close_coupling_segmented_logder( &
        grid, v_mat, mass, total_energy, thresholds, l_channels, res, stat)

        type(segmented_grid_t), intent(in)      :: grid
        real(dp), dimension(:, :, :), intent(in):: v_mat
        real(dp), intent(in)                    :: mass
        real(dp), intent(in)                    :: total_energy
        real(dp), dimension(:), intent(in)      :: thresholds
        integer,  dimension(:), intent(in)      :: l_channels
        type(multichannel_result_t), intent(out):: res
        integer, optional, intent(out)          :: stat

        integer  :: n_chan, n_pts, sec, step, idx, i, j, stat_inv
        integer  :: j_start, j_end, n_pts_sec
        real(dp) :: h, h2_12, r_curr, e_kin, k_i, r_match
        real(dp) :: jl, nl, djl, dnl
        real(dp), allocatable :: w_mat(:, :), q_mat(:, :), q_inv(:, :), m_mat(:, :)
        real(dp), allocatable :: r_curr_mat(:, :), r_next_mat(:, :), r_inv(:, :)
        real(dp), allocatable :: q_prev(:, :), q_curr(:, :), p1(:, :), p2(:, :)
        real(dp), allocatable :: temp_mat(:, :), y_mat(:, :), y_eff(:, :)
        real(dp), allocatable :: y_oo(:, :), y_oc(:, :), y_co(:, :), y_cc(:, :)
        real(dp), allocatable :: a_cc(:, :), a_cc_inv(:, :)
        real(dp), allocatable :: j_mat(:, :), n_mat(:, :), dj_mat(:, :), dn_mat(:, :)
        real(dp), allocatable :: mj_mat(:, :), mn_mat(:, :), mn_inv(:, :)
        complex(dp), allocatable :: eye_c(:, :), ik_mat(:, :), den_c(:, :), den_inv(:, :)

        if (present(stat)) stat = 0
        n_chan = size(thresholds)
        n_pts  = grid%n_total
        res%n_channels = n_chan
        res%n_open = 0
        res%n_closed = 0

        if (n_chan < 1 .or. n_pts < 5) then
            if (present(stat)) stat = -1
            return
        end if

        ! 1. 统计并分类开通道与闭通道
        allocate(res%open_channels(n_chan))
        allocate(res%closed_channels(n_chan))

        do i = 1, n_chan
            e_kin = total_energy - thresholds(i)
            if (e_kin > 1.0e-13_dp) then
                res%n_open = res%n_open + 1
                res%open_channels(res%n_open) = i
            else
                res%n_closed = res%n_closed + 1
                res%closed_channels(res%n_closed) = i
            end if
        end do

        if (res%n_open == 0) then
            if (present(stat)) stat = 1
            return
        end if

        allocate(res%k_open(res%n_open))
        do i = 1, res%n_open
            res%k_open(i) = sqrt(2.0_dp * mass * (total_energy - thresholds(res%open_channels(i))))
        end do

        if (res%n_closed > 0) then
            allocate(res%kappa_closed(res%n_closed))
            do i = 1, res%n_closed
                res%kappa_closed(i) = sqrt(2.0_dp * mass * max(0.0_dp, thresholds(res%closed_channels(i)) - total_energy))
            end do
        end if

        ! 2. 初始化矩阵工作空间
        allocate(w_mat(n_chan, n_chan))
        allocate(q_mat(n_chan, n_chan))
        allocate(q_inv(n_chan, n_chan))
        allocate(m_mat(n_chan, n_chan))
        allocate(r_curr_mat(n_chan, n_chan))
        allocate(r_next_mat(n_chan, n_chan))
        allocate(r_inv(n_chan, n_chan))
        allocate(q_prev(n_chan, n_chan))
        allocate(q_curr(n_chan, n_chan))
        allocate(p1(n_chan, n_chan))
        allocate(p2(n_chan, n_chan))
        allocate(temp_mat(n_chan, n_chan))
        allocate(y_mat(n_chan, n_chan))

        ! 3. 逐扇区递推推进
        do sec = 1, grid%n_sectors
            h = grid%sector_dr(sec)
            h2_12 = (h * h) / 12.0_dp
            j_start   = grid%sector_offset(sec)
            n_pts_sec = grid%sector_npts(sec)
            j_end     = j_start + n_pts_sec - 1

            if (sec == 1) then
                ! 扇区 1 从禁区原点初始化
                r_curr = grid%r(1)
                do i = 1, n_chan
                    do j = 1, n_chan
                        w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, 1)
                        if (i == j) then
                            w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                          real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                        end if
                    end do
                end do

                q_mat = -h2_12 * w_mat
                do i = 1, n_chan
                    q_mat(i, i) = q_mat(i, i) + 1.0_dp
                end do

                call inv_real_matrix(n_chan, q_mat, q_inv, stat_inv)
                if (stat_inv /= 0) then
                    if (present(stat)) stat = -2
                    return
                end if

                m_mat = 12.0_dp * q_inv
                do i = 1, n_chan
                    m_mat(i, i) = m_mat(i, i) - 10.0_dp
                end do

                r_curr_mat = m_mat
                q_prev = q_mat
            else
                ! 扇区 sec > 1: 利用交界面 y_mat 与新步长 h 构建新扇区启动比值矩阵
                r_curr = grid%r(j_start)
                do i = 1, n_chan
                    do j = 1, n_chan
                        w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, j_start)
                        if (i == j) then
                            w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                          real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                        end if
                    end do
                end do

                q_mat = -h2_12 * w_mat
                do i = 1, n_chan
                    q_mat(i, i) = q_mat(i, i) + 1.0_dp
                end do

                call inv_real_matrix(n_chan, q_mat, q_inv, stat_inv)
                if (stat_inv /= 0) then
                    if (present(stat)) stat = -2
                    return
                end if

                m_mat = 12.0_dp * q_inv
                do i = 1, n_chan
                    m_mat(i, i) = m_mat(i, i) - 10.0_dp
                end do

                ! R_start = I + h * Y - 0.5 * h^2 * W
                r_curr_mat = h * y_mat - 0.5_dp * (h * h) * w_mat
                do i = 1, n_chan
                    r_curr_mat(i, i) = r_curr_mat(i, i) + 1.0_dp
                end do

                call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)
                if (stat_inv /= 0) then
                    do i = 1, n_chan
                        r_curr_mat(i, i) = r_curr_mat(i, i) + 1.0e-12_dp
                    end do
                    call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)
                end if

                r_curr_mat = m_mat - r_inv
                r_curr_mat = 0.5_dp * (r_curr_mat + transpose(r_curr_mat))
                q_prev = q_mat
            end if

            ! 扇区内部 Johnson 比值递推
            do step = 2, n_pts_sec - 1
                idx = j_start + step - 1
                r_curr = grid%r(idx)

                do i = 1, n_chan
                    do j = 1, n_chan
                        w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, idx)
                        if (i == j) then
                            w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                          real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                        end if
                    end do
                end do

                q_mat = -h2_12 * w_mat
                do i = 1, n_chan
                    q_mat(i, i) = q_mat(i, i) + 1.0_dp
                end do

                call inv_real_matrix(n_chan, q_mat, q_inv, stat_inv)
                if (stat_inv /= 0) then
                    if (present(stat)) stat = -3
                    return
                end if

                m_mat = 12.0_dp * q_inv
                do i = 1, n_chan
                    m_mat(i, i) = m_mat(i, i) - 10.0_dp
                end do

                call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)
                if (stat_inv /= 0) then
                    do i = 1, n_chan
                        r_curr_mat(i, i) = r_curr_mat(i, i) + 1.0e-14_dp
                    end do
                    call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)
                end if

                r_next_mat = m_mat - r_inv
                r_next_mat = 0.5_dp * (r_next_mat + transpose(r_next_mat))

                if (step == n_pts_sec - 2) then
                    q_curr = q_mat
                end if

                r_curr_mat = r_next_mat
            end do

            ! 扇区终止边界提取局部对数导数矩阵 Y(r_{j_end})
            call inv_real_matrix(n_chan, r_curr_mat, r_inv, stat_inv)

            r_curr = grid%r(j_end)
            do i = 1, n_chan
                do j = 1, n_chan
                    w_mat(j, i) = 2.0_dp * mass * v_mat(j, i, j_end)
                    if (i == j) then
                        w_mat(i, i) = w_mat(i, i) - 2.0_dp * mass * (total_energy - thresholds(i)) + &
                                      real(l_channels(i) * (l_channels(i) + 1), dp) / (r_curr * r_curr)
                    end if
                end do
            end do

            q_mat = -h2_12 * w_mat
            do i = 1, n_chan
                q_mat(i, i) = q_mat(i, i) + 1.0_dp
            end do

            call inv_real_matrix(n_chan, q_curr, q_inv, stat_inv)
            temp_mat = matmul(r_inv, q_mat)
            p1 = matmul(q_inv, temp_mat)
            p2 = matmul(p1, p1)

            y_mat = p2 - 4.0_dp * p1
            do i = 1, n_chan
                y_mat(i, i) = y_mat(i, i) + 3.0_dp
            end do
            y_mat = y_mat / (2.0_dp * h)
            y_mat = 0.5_dp * (y_mat + transpose(y_mat))
        end do

        ! 4. 开/闭通道 Schur 补变换与 Feshbach 投影
        allocate(y_eff(res%n_open, res%n_open))
        if (res%n_closed == 0) then
            y_eff = y_mat
        else
            allocate(y_oo(res%n_open, res%n_open))
            allocate(y_oc(res%n_open, res%n_closed))
            allocate(y_co(res%n_closed, res%n_open))
            allocate(y_cc(res%n_closed, res%n_closed))
            allocate(a_cc(res%n_closed, res%n_closed))
            allocate(a_cc_inv(res%n_closed, res%n_closed))

            do i = 1, res%n_open
                do j = 1, res%n_open
                    y_oo(j, i) = y_mat(res%open_channels(j), res%open_channels(i))
                end do
            end do

            do i = 1, res%n_closed
                do j = 1, res%n_open
                    y_oc(j, i) = y_mat(res%open_channels(j), res%closed_channels(i))
                    y_co(i, j) = y_mat(res%closed_channels(i), res%open_channels(j))
                end do
            end do

            do i = 1, res%n_closed
                do j = 1, res%n_closed
                    y_cc(j, i) = y_mat(res%closed_channels(j), res%closed_channels(i))
                end do
            end do

            a_cc = y_cc
            do i = 1, res%n_closed
                a_cc(i, i) = a_cc(i, i) + res%kappa_closed(i)
            end do

            call inv_real_matrix(res%n_closed, a_cc, a_cc_inv, stat_inv)
            if (stat_inv /= 0) then
                do i = 1, res%n_closed
                    a_cc(i, i) = a_cc(i, i) + 1.0e-10_dp
                end do
                call inv_real_matrix(res%n_closed, a_cc, a_cc_inv, stat_inv)
            end if

            y_eff = y_oo - matmul(y_oc, matmul(a_cc_inv, y_co))
            y_eff = 0.5_dp * (y_eff + transpose(y_eff))

            deallocate(y_oo, y_oc, y_co, y_cc, a_cc, a_cc_inv)
        end if

        ! 5. 渐近外边界 Riccati 函数匹配提取开通道反应矩阵 K_oo
        r_match = grid%r(n_pts)
        allocate(j_mat(res%n_open, res%n_open))
        allocate(n_mat(res%n_open, res%n_open))
        allocate(dj_mat(res%n_open, res%n_open))
        allocate(dn_mat(res%n_open, res%n_open))
        allocate(mj_mat(res%n_open, res%n_open))
        allocate(mn_mat(res%n_open, res%n_open))
        allocate(mn_inv(res%n_open, res%n_open))
        allocate(res%k_matrix(res%n_open, res%n_open))

        j_mat = 0.0_dp; n_mat = 0.0_dp; dj_mat = 0.0_dp; dn_mat = 0.0_dp

        do i = 1, res%n_open
            k_i = res%k_open(i)
            call riccati_bessel_neumann(l_channels(res%open_channels(i)), k_i * r_match, jl, nl, djl, dnl)
            j_mat(i, i)  = jl / sqrt(k_i)
            n_mat(i, i)  = nl / sqrt(k_i)
            dj_mat(i, i) = sqrt(k_i) * djl
            dn_mat(i, i) = sqrt(k_i) * dnl
        end do

        mj_mat = dj_mat - matmul(y_eff, j_mat)
        mn_mat = dn_mat - matmul(y_eff, n_mat)

        call inv_real_matrix(res%n_open, mn_mat, mn_inv, stat_inv)
        if (stat_inv /= 0) then
            do i = 1, res%n_open
                mn_mat(i, i) = mn_mat(i, i) + 1.0e-12_dp
            end do
            call inv_real_matrix(res%n_open, mn_mat, mn_inv, stat_inv)
        end if

        res%k_matrix = matmul(mn_inv, mj_mat)
        res%k_matrix = 0.5_dp * (res%k_matrix + transpose(res%k_matrix))

        ! 6. Cayley 变换求解幺正散射矩阵 S = (I + i*K)(I - i*K)^{-1}
        allocate(eye_c(res%n_open, res%n_open))
        allocate(ik_mat(res%n_open, res%n_open))
        allocate(den_c(res%n_open, res%n_open))
        allocate(den_inv(res%n_open, res%n_open))
        allocate(res%s_matrix(res%n_open, res%n_open))
        allocate(res%t_matrix(res%n_open, res%n_open))
        allocate(res%prob_matrix(res%n_open, res%n_open))
        allocate(res%cross_sections(res%n_open, res%n_open))
        allocate(res%total_cross_sec(res%n_open))

        eye_c = (0.0_dp, 0.0_dp)
        do i = 1, res%n_open
            eye_c(i, i) = (1.0_dp, 0.0_dp)
        end do

        do i = 1, res%n_open
            do j = 1, res%n_open
                ik_mat(j, i) = cmplx(0.0_dp, res%k_matrix(j, i), kind=dp)
            end do
        end do

        den_c = eye_c - ik_mat
        call inv_complex_matrix(res%n_open, den_c, den_inv, stat_inv)
        res%s_matrix = matmul(eye_c + ik_mat, den_inv)
        res%t_matrix = res%s_matrix - eye_c

        ! 7. 跃迁几率、态-态截面与特征相移和
        res%total_cross_sec = 0.0_dp
        do i = 1, res%n_open
            k_i = res%k_open(i)
            do j = 1, res%n_open
                res%prob_matrix(j, i) = abs(res%s_matrix(j, i))**2
                res%cross_sections(j, i) = (PI / (k_i * k_i)) * abs(res%t_matrix(j, i))**2
                res%total_cross_sec(i) = res%total_cross_sec(i) + res%cross_sections(j, i)
            end do
        end do

        res%eigenphase_sum = 0.0_dp
        do i = 1, res%n_open
            res%eigenphase_sum = res%eigenphase_sum + atan(res%k_matrix(i, i))
        end do

        ! 释放内存
        deallocate(w_mat, q_mat, q_inv, m_mat, r_curr_mat, r_next_mat, r_inv)
        deallocate(q_prev, q_curr, p1, p2, temp_mat, y_mat, y_eff)
        deallocate(j_mat, n_mat, dj_mat, dn_mat, mj_mat, mn_mat, mn_inv)
        deallocate(eye_c, ik_mat, den_c, den_inv)
    end subroutine calc_multichannel_close_coupling_segmented_logder

end module mod_ti_scattering
