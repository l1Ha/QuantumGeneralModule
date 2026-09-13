!> \brief 强场非顺序双电离 (NSDI) 与电子重碰撞相关动量动力学模块
!> \details 基于 Corkum 三步模型 (Three-Step Model)、ADK 准静态隧穿电离理论、
!>          Lotz 碰撞电离截面与两电子经典轨道数值积分，
!>          计算最大回碰动能 (3.17 Up 截断)、(e,2e) 直接碰撞与 RESI 激发机制、
!>          平行双电子二维相关动量谱 P(pz1, pz2) 以及双电离产率“膝盖结构” (Knee Structure)。
!> \author LiHao
!> \date 2026-09-13
module mod_strong_field_nsdi
    use mod_constants, only: dp, PI, TWOPI, HALFPI, C_LIGHT, AU2W_CM2, W_CM2AU, NM2AU, EV2AU
    implicit none
    private

    public :: nsdi_laser_t, nsdi_target_t, nsdi_result_t
    public :: init_nsdi_laser, init_nsdi_target
    public :: calc_ponderomotive_energy, calc_keldysh_gamma, calc_adk_rate
    public :: calc_recollision_trajectory, calc_lotz_cross_section
    public :: calc_nsdi_drift_momenta, calc_nsdi_2d_momentum_dist
    public :: calc_double_ion_yield_curve

    !> \brief 强激光场参数类型
    type :: nsdi_laser_t
        real(dp) :: wavelength_nm   = 800.0_dp          !< 激光波长 (nm)
        real(dp) :: intensity_w_cm2 = 2.0e14_dp         !< 激光峰值光强 (W/cm^2)
        real(dp) :: omega_au        = 0.05696_dp        !< 载波中心角频率 (a.u.)
        real(dp) :: field_peak_au   = 0.07548_dp        !< 峰值电场振幅 F0 (a.u.)
        real(dp) :: up_au           = 0.4397_dp         !< 有质动力势 Up = F0^2 / (4*omega^2) (a.u.)
    end type nsdi_laser_t

    !> \brief 目标原子参数类型 (He, Ne, Ar, Xe 等)
    type :: nsdi_target_t
        character(len=16) :: name       = "He"          !< 目标名称
        real(dp)          :: ip1_au     = 0.90355_dp    !< 第一电离势 Ip1 (a.u.), He: 24.59 eV
        real(dp)          :: ip2_au     = 2.00000_dp    !< 第二电离势 Ip2 (a.u.), He+: 54.42 eV
        real(dp)          :: ip_exc_au  = 1.50000_dp    !< 最低激发态能级差 (用于 RESI 激发通道)
        real(dp)          :: z_core     = 1.0_dp        !< 母体离子有效核电荷
    end type nsdi_target_t

    !> \brief NSDI 计算结果汇总
    type :: nsdi_result_t
        real(dp) :: yield_nsdi     = 0.0_dp             !< 非顺序双电离产率
        real(dp) :: yield_sdi      = 0.0_dp             !< 顺序双电离产率
        real(dp) :: ratio_nsdi_sdi = 0.0_dp             !< NSDI / SDI 产率比
        real(dp) :: e_rec_max_au   = 0.0_dp             !< 经典最大回碰动能 (3.17 Up)
        real(dp) :: p_rec_max_au   = 0.0_dp             !< 最大回碰动量 sqrt(2 * E_rec_max)
        real(dp) :: corr_coeff     = 0.0_dp             !< 双电子相关系数 <p1*p2>/sqrt(<p1^2><p2^2>)
    end type nsdi_result_t

contains

    ! ==========================================================================
    ! init_nsdi_laser: 初始化激光场参数
    ! ==========================================================================
    subroutine init_nsdi_laser(wavelength_nm, intensity_w_cm2, laser, stat)
        real(dp), intent(in)            :: wavelength_nm
        real(dp), intent(in)            :: intensity_w_cm2
        type(nsdi_laser_t), intent(out) :: laser
        integer, intent(out)            :: stat

        real(dp) :: lambda_au, intensity_au

        stat = 0
        if (wavelength_nm <= 0.0_dp .or. intensity_w_cm2 <= 0.0_dp) then
            stat = 1
            return
        end if

        laser%wavelength_nm   = wavelength_nm
        laser%intensity_w_cm2 = intensity_w_cm2

        ! 波长转换为原子单位: lambda (Bohr)
        lambda_au = wavelength_nm * NM2AU

        ! 载波角频率 omega = 2 * pi * c / lambda (a.u.)
        laser%omega_au = TWOPI * 137.035999084_dp / lambda_au

        ! 激光光强在原子单位下: I_au = I (W/cm^2) * W_CM2AU
        intensity_au = intensity_w_cm2 * W_CM2AU
        laser%field_peak_au = sqrt(intensity_au)

        ! 激光有质动力能量 Up = F0^2 / (4 * omega^2)
        laser%up_au = calc_ponderomotive_energy(laser)
    end subroutine init_nsdi_laser

    ! ==========================================================================
    ! init_nsdi_target: 初始化目标原子参数
    ! ==========================================================================
    subroutine init_nsdi_target(target_name, target, stat)
        character(len=*), intent(in)     :: target_name
        type(nsdi_target_t), intent(out) :: target
        integer, intent(out)             :: stat

        stat = 0
        target%name = trim(target_name)

        if (trim(target_name) == "He" .or. trim(target_name) == "he") then
            target%ip1_au    = 24.587_dp * EV2AU   ! ~ 0.90355 a.u.
            target%ip2_au    = 54.418_dp * EV2AU   ! ~ 2.00000 a.u.
            target%ip_exc_au = 40.814_dp * EV2AU   ! 1s -> 2p 激发能
            target%z_core    = 1.0_dp
        else if (trim(target_name) == "Ne" .or. trim(target_name) == "ne") then
            target%ip1_au    = 21.564_dp * EV2AU   ! ~ 0.79246 a.u.
            target%ip2_au    = 40.963_dp * EV2AU   ! ~ 1.50536 a.u.
            target%ip_exc_au = 26.910_dp * EV2AU
            target%z_core    = 1.0_dp
        else if (trim(target_name) == "Ar" .or. trim(target_name) == "ar") then
            target%ip1_au    = 15.760_dp * EV2AU   ! ~ 0.57917 a.u.
            target%ip2_au    = 27.630_dp * EV2AU   ! ~ 1.01538 a.u.
            target%ip_exc_au = 16.420_dp * EV2AU
            target%z_core    = 1.0_dp
        else if (trim(target_name) == "Xe" .or. trim(target_name) == "xe") then
            target%ip1_au    = 12.130_dp * EV2AU   ! ~ 0.44577 a.u.
            target%ip2_au    = 20.975_dp * EV2AU   ! ~ 0.77082 a.u.
            target%ip_exc_au = 11.270_dp * EV2AU
            target%z_core    = 1.0_dp
        else
            stat = -1  ! 未知原子，采用默认 Helium 参数
            target%ip1_au    = 24.587_dp * EV2AU
            target%ip2_au    = 54.418_dp * EV2AU
            target%ip_exc_au = 40.814_dp * EV2AU
            target%z_core    = 1.0_dp
        end if
    end subroutine init_nsdi_target

    ! ==========================================================================
    ! calc_ponderomotive_energy: 计算有质动力能量 Up = F0^2 / (4 * omega^2)
    ! ==========================================================================
    pure function calc_ponderomotive_energy(laser) result(up_au)
        type(nsdi_laser_t), intent(in) :: laser
        real(dp) :: up_au

        if (laser%omega_au > 0.0_dp) then
            up_au = (laser%field_peak_au**2) / (4.0_dp * (laser%omega_au**2))
        else
            up_au = 0.0_dp
        end if
    end function calc_ponderomotive_energy

    ! ==========================================================================
    ! calc_keldysh_gamma: 计算 Keldysh 参数 gamma = omega * sqrt(2 * Ip) / F0
    ! gamma < 1 时隧穿机制主导，gamma > 1 时多光子电离主导
    ! ==========================================================================
    pure function calc_keldysh_gamma(laser, ip_au) result(gamma_k)
        type(nsdi_laser_t), intent(in) :: laser
        real(dp), intent(in)           :: ip_au
        real(dp) :: gamma_k

        if (laser%field_peak_au > 0.0_dp .and. ip_au > 0.0_dp) then
            gamma_k = (laser%omega_au * sqrt(2.0_dp * ip_au)) / laser%field_peak_au
        else
            gamma_k = 1.0e10_dp
        end if
    end function calc_keldysh_gamma

    ! ==========================================================================
    ! calc_adk_rate: 计算准静态隧穿电离率 (Ammosov-Delone-Krainov 理论)
    ! w(F) = C_n*^2 * (2 * F_c / F)^(2*n* - 1) * exp(-2 * F_c / (3 * F))
    ! ==========================================================================
    pure function calc_adk_rate(field_au, ip_au, z_core) result(rate)
        real(dp), intent(in) :: field_au
        real(dp), intent(in) :: ip_au
        real(dp), intent(in) :: z_core
        real(dp) :: rate

        real(dp) :: f_abs, n_star, f_c, prefactor, expo

        f_abs = abs(field_au)
        if (f_abs < 1.0e-12_dp .or. ip_au <= 0.0_dp) then
            rate = 0.0_dp
            return
        end if

        ! 有效主量子数 n* = Z / sqrt(2 * Ip)
        n_star = z_core / sqrt(2.0_dp * ip_au)

        ! 特征库仑电场 F_c = (2 * Ip)^(3/2)
        f_c = (2.0_dp * ip_au)**1.5_dp

        expo = -2.0_dp * f_c / (3.0_dp * f_abs)
        if (expo < -200.0_dp) then
            rate = 0.0_dp
            return
        end if

        ! ADK 前因子 (取 s 态近似 l = 0, m = 0)
        prefactor = ((2.0_dp * f_c / f_abs)**(2.0_dp * n_star - 1.0_dp)) / &
                    (n_star**4.0_dp)

        rate = prefactor * exp(expo)
    end function calc_adk_rate

    ! ==========================================================================
    ! calc_recollision_trajectory: 求解经典电子轨道与母核回碰时刻
    ! 电子在出生相位 phi_0 = omega * t_0 处以初速度零脱离原子核
    ! x(phi) = (F0 / omega^2) * [ cos(phi) - cos(phi_0) + sin(phi_0) * (phi - phi_0) ]
    ! 求解非平凡根 x(phi_r) = 0 (phi_r > phi_0)
    ! 计算回碰动能 E_rec(phi_0) = 2 * Up * (sin(phi_r) - sin(phi_0))^2
    ! ==========================================================================
    pure subroutine calc_recollision_trajectory(phi_0, up_au, phi_r, e_rec_au, recollided)
        real(dp), intent(in)  :: phi_0
        real(dp), intent(in)  :: up_au
        real(dp), intent(out) :: phi_r
        real(dp), intent(out) :: e_rec_au
        logical, intent(out)  :: recollided

        real(dp) :: phi_curr, phi_next, f_val, df_val
        integer  :: iter

        recollided = .false.
        phi_r = 0.0_dp
        e_rec_au = 0.0_dp

        ! 只有在激光电场极性反转前出生的电子才可能在反向电场加速下返回母核
        ! 在一个周期内，回碰电子主要出生在相位 phi_0 ∈ (0, pi/2)
        if (phi_0 <= 0.0_dp .or. phi_0 >= HALFPI) return

        ! 牛顿-拉夫逊法求解回碰相位 phi_r
        ! 初值选取: 经典回碰通常发生在大半个光周期后 phi ~ 4.0 - 4.5 rad
        phi_curr = max(3.5_dp, phi_0 + PI)

        do iter = 1, 50
            f_val  = cos(phi_curr) - cos(phi_0) + sin(phi_0) * (phi_curr - phi_0)
            df_val = -sin(phi_curr) + sin(phi_0)

            if (abs(df_val) < 1.0e-10_dp) df_val = 1.0e-10_dp
            phi_next = phi_curr - f_val / df_val

            ! 限制在合理物理区间内
            if (phi_next <= phi_0) phi_next = phi_0 + 0.5_dp * PI
            if (phi_next > TWOPI + phi_0) phi_next = TWOPI + phi_0

            if (abs(phi_next - phi_curr) < 1.0e-9_dp) then
                phi_curr = phi_next
                exit
            end if
            phi_curr = phi_next
        end do

        ! 验证是否为回碰根
        f_val = cos(phi_curr) - cos(phi_0) + sin(phi_0) * (phi_curr - phi_0)
        if (abs(f_val) < 1.0e-5_dp .and. phi_curr > phi_0) then
            recollided = .true.
            phi_r = phi_curr
            ! 回碰动能: E_rec = 1/2 * v^2 = 2 * Up * (sin(phi_r) - sin(phi_0))^2
            e_rec_au = 2.0_dp * up_au * ((sin(phi_r) - sin(phi_0))**2)
        end if
    end subroutine calc_recollision_trajectory

    ! ==========================================================================
    ! calc_lotz_cross_section: 电子碰撞电离 (e,2e) Lotz 经验截面公式
    ! sigma(E) = a * q * ln(E / Ip2) / (E * Ip2), 当 E > Ip2; 否则为 0
    ! ==========================================================================
    pure function calc_lotz_cross_section(e_rec_au, ip2_au) result(sigma)
        real(dp), intent(in) :: e_rec_au
        real(dp), intent(in) :: ip2_au
        real(dp) :: sigma

        real(dp) :: e_ratio

        if (e_rec_au <= ip2_au .or. ip2_au <= 0.0_dp) then
            sigma = 0.0_dp
            return
        end if

        e_ratio = e_rec_au / ip2_au
        ! 无量纲/原子单位标度碰撞电离截面
        sigma = log(e_ratio) / (e_rec_au * ip2_au)
    end function calc_lotz_cross_section

    ! ==========================================================================
    ! calc_nsdi_drift_momenta: 计算双电子在强激光场中的最终漂移动量
    ! 在碰撞时刻 t_r 发生直接 (e,2e) 碰撞电离，
    ! 多余动能 Delta_E = E_rec - Ip2 被双电子按比例 alpha 分享
    ! 碰撞后电场继续做功产生漂移动量 p_drift = -A(t_r) = (F0 / omega) * sin(phi_r)
    ! pz1 = p_share1 + p_drift,  pz2 = p_share2 + p_drift
    ! ==========================================================================
    pure subroutine calc_nsdi_drift_momenta(phi_r, e_rec_au, ip2_au, f0_au, omega_au, &
                                            alpha_share, pz1, pz2)
        real(dp), intent(in)  :: phi_r
        real(dp), intent(in)  :: e_rec_au
        real(dp), intent(in)  :: ip2_au
        real(dp), intent(in)  :: f0_au
        real(dp), intent(in)  :: omega_au
        real(dp), intent(in)  :: alpha_share  ! 能量分配系数 0 <= alpha <= 1
        real(dp), intent(out) :: pz1, pz2

        real(dp) :: delta_e, p_drift, k1, k2

        delta_e = max(0.0_dp, e_rec_au - ip2_au)
        ! 强场矢势赋予的共同漂移动量
        p_drift = (f0_au / omega_au) * sin(phi_r)

        ! 双电子分享出射初动量
        k1 = sqrt(2.0_dp * alpha_share * delta_e)
        k2 = sqrt(2.0_dp * (1.0_dp - alpha_share) * delta_e)

        ! 激光偏振轴方向总末动量
        pz1 = p_drift + k1
        pz2 = p_drift + k2
    end subroutine calc_nsdi_drift_momenta

    ! ==========================================================================
    ! calc_nsdi_2d_momentum_dist:
    ! 数值积分生成双电子二维纵向相关动量谱 P(pz1, pz2) 并计算关联系数
    ! ==========================================================================
    subroutine calc_nsdi_2d_momentum_dist(laser, target, n_pts, p_max, grid_p, dist_2d, corr_coeff)
        type(nsdi_laser_t), intent(in)     :: laser
        type(nsdi_target_t), intent(in)    :: target
        integer, intent(in)                :: n_pts
        real(dp), intent(in)               :: p_max
        real(dp), intent(out)              :: grid_p(n_pts)
        real(dp), intent(out)              :: dist_2d(n_pts, n_pts)
        real(dp), intent(out)              :: corr_coeff

        integer  :: i, j, i_phi, i_share
        real(dp) :: dp_step, phi_0, field_inst, w_tunnel, phi_r, e_rec
        real(dp) :: sigma_coll, alpha, pz1, pz2, weight
        real(dp) :: sum_weight, sum_p1p2, sum_p1_sq, sum_p2_sq
        real(dp) :: sigma_p, dist_sum
        logical  :: ok

        ! 构建动量网格
        dp_step = 2.0_dp * p_max / real(n_pts - 1, dp)
        do i = 1, n_pts
            grid_p(i) = -p_max + real(i - 1, dp) * dp_step
        end do
        dist_2d = 0.0_dp

        sum_weight = 0.0_dp
        sum_p1p2   = 0.0_dp
        sum_p1_sq  = 0.0_dp
        sum_p2_sq  = 0.0_dp

        sigma_p = max(0.15_dp, 0.08_dp * p_max)

        ! 遍历隧穿出生相位 phi_0 ∈ (0.02, pi/2 - 0.02)
        do i_phi = 1, 100
            phi_0 = 0.02_dp + real(i_phi - 1, dp) * (HALFPI - 0.04_dp) / 99.0_dp
            field_inst = laser%field_peak_au * cos(phi_0)
            w_tunnel   = calc_adk_rate(field_inst, target%ip1_au, target%z_core)

            if (w_tunnel <= 1.0e-30_dp) cycle

            call calc_recollision_trajectory(phi_0, laser%up_au, phi_r, e_rec, ok)
            if (.not. ok) cycle

            ! 碰撞电离截面
            sigma_coll = calc_lotz_cross_section(e_rec, target%ip2_au)
            if (sigma_coll <= 1.0e-30_dp) cycle

            ! 遍历能量共享比例 alpha ∈ [0, 1]
            do i_share = 1, 11
                alpha = real(i_share - 1, dp) / 10.0_dp
                call calc_nsdi_drift_momenta(phi_r, e_rec, target%ip2_au, &
                                             laser%field_peak_au, laser%omega_au, &
                                             alpha, pz1, pz2)

                weight = w_tunnel * sigma_coll

                ! 对称性: 激光反向半周期贡献 (-pz1, -pz2)
                ! 将高斯展宽沉积到二维网格
                do i = 1, n_pts
                    if (abs(grid_p(i) - pz1) > 3.0_dp * sigma_p .and. &
                        abs(grid_p(i) + pz1) > 3.0_dp * sigma_p) cycle

                    do j = 1, n_pts
                        dist_2d(i, j) = dist_2d(i, j) + weight * ( &
                            exp(-((grid_p(i) - pz1)**2 + (grid_p(j) - pz2)**2) / (2.0_dp * sigma_p**2)) + &
                            exp(-((grid_p(i) + pz1)**2 + (grid_p(j) + pz2)**2) / (2.0_dp * sigma_p**2)) &
                        )
                    end do
                end do

                sum_weight = sum_weight + 2.0_dp * weight
                sum_p1p2   = sum_p1p2   + 2.0_dp * weight * (pz1 * pz2)
                sum_p1_sq  = sum_p1_sq  + 2.0_dp * weight * (pz1**2)
                sum_p2_sq  = sum_p2_sq  + 2.0_dp * weight * (pz2**2)
            end do
        end do

        ! 归一化二维谱
        dist_sum = sum(dist_2d) * (dp_step**2)
        if (dist_sum > 1.0e-30_dp) then
            dist_2d = dist_2d / dist_sum
        end if

        ! 计算统计关联系数 C = <p1 * p2> / sqrt(<p1^2> * <p2^2>)
        if (sum_weight > 1.0e-30_dp .and. sum_p1_sq > 0.0_dp .and. sum_p2_sq > 0.0_dp) then
            corr_coeff = (sum_p1p2 / sum_weight) / sqrt((sum_p1_sq / sum_weight) * (sum_p2_sq / sum_weight))
        else
            corr_coeff = 0.0_dp
        end if
    end subroutine calc_nsdi_2d_momentum_dist

    ! ==========================================================================
    ! calc_double_ion_yield_curve:
    ! 计算光强依赖的非顺序 (NSDI) 与顺序 (SDI) 双电离产率曲线，重现著名的“膝盖结构”
    ! ==========================================================================
    subroutine calc_double_ion_yield_curve(wavelength_nm, target, n_int, &
                                           int_min, int_max, intensities, &
                                           yields_nsdi, yields_sdi)
        real(dp), intent(in)            :: wavelength_nm
        type(nsdi_target_t), intent(in) :: target
        integer, intent(in)             :: n_int
        real(dp), intent(in)            :: int_min
        real(dp), intent(in)            :: int_max
        real(dp), intent(out)           :: intensities(n_int)
        real(dp), intent(out)           :: yields_nsdi(n_int)
        real(dp), intent(out)           :: yields_sdi(n_int)

        integer           :: k, i_phi, stat
        real(dp)          :: log_i_min, log_i_max, log_i
        type(nsdi_laser_t):: laser
        real(dp)          :: phi_0, field_inst, w1, w2, phi_r, e_rec, sigma_coll
        real(dp)          :: nsdi_integral, sdi_integral
        logical           :: ok

        log_i_min = log10(int_min)
        log_i_max = log10(int_max)

        do k = 1, n_int
            log_i = log_i_min + real(k - 1, dp) * (log_i_max - log_i_min) / real(n_int - 1, dp)
            intensities(k) = 10.0_dp**log_i

            call init_nsdi_laser(wavelength_nm, intensities(k), laser, stat)

            nsdi_integral = 0.0_dp
            sdi_integral  = 0.0_dp

            ! 积分周期内的电离几率
            do i_phi = 1, 100
                phi_0 = 0.02_dp + real(i_phi - 1, dp) * (HALFPI - 0.04_dp) / 99.0_dp
                field_inst = laser%field_peak_au * cos(phi_0)

                w1 = calc_adk_rate(field_inst, target%ip1_au, 1.0_dp)
                w2 = calc_adk_rate(field_inst, target%ip2_au, 2.0_dp)

                ! 顺序双电离: 第二电子直接隧穿电离 (乘积积分)
                sdi_integral = sdi_integral + w1 * w2

                ! 非顺序双电离: 第一电子隧穿后重碰撞诱导第二电子碰撞电离
                call calc_recollision_trajectory(phi_0, laser%up_au, phi_r, e_rec, ok)
                if (ok) then
                    sigma_coll = calc_lotz_cross_section(e_rec, target%ip2_au)
                    ! 若低于直接电离阈值，考虑 RESI (碰撞激发-后续场致电离) 通道
                    if (sigma_coll <= 0.0_dp .and. e_rec > (target%ip2_au - target%ip_exc_au)) then
                        sigma_coll = 0.25_dp * log(e_rec / (target%ip2_au - target%ip_exc_au)) / &
                                     (e_rec * (target%ip2_au - target%ip_exc_au))
                    end if
                    nsdi_integral = nsdi_integral + w1 * sigma_coll
                end if
            end do

            ! 产率标度
            yields_sdi(k)  = sdi_integral * 1.0e-3_dp
            yields_nsdi(k) = nsdi_integral * 0.05_dp
        end do
    end subroutine calc_double_ion_yield_curve

end module mod_strong_field_nsdi
