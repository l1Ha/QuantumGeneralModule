! ==============================================================================
! GeneralModule: 现代量子动力学与原子分子物理通用算法库
! 模块: mod_penning_associative_ionization
! 描述: 亚稳态激发原子碰撞潘宁电离 (Penning Ionization) 与缔合电离
!       (Associative Ionization / Chemi-ionization) 动力学模块。
!       包含复光学势半经典散射、经典转折点与电离存活概率、
!       分流比计算 (PI vs AI)、潘宁电子能谱 (PIES)、热平衡反应速率常数
!       与超冷 Wigner 阈值复散射长度 (Complex Scattering Length)。
! 标准: Fortran 2008 / 纯过程与强类型规范 / 零外部库强依赖
! ==============================================================================

module mod_penning_associative_ionization
    use iso_fortran_env, only: dp => real64
    use mod_constants, only: PI, TWOPI, HBAR, AMU2AU, AU2ANG, ANG2AU, &
                             EV2AU, AU2EV, FS2AU, AU2FS
    implicit none
    private

    ! --------------------------------------------------------------------------
    ! 派生类型定义 (Derived Types)
    ! --------------------------------------------------------------------------

    !> \brief 潘宁与缔合电离体系参数结构体
    type, public :: penning_system_t
        character(len=32) :: name             !< 体系名称 (如 "He*(2^3S) + Ar")
        real(dp)          :: m_a_amu          !< 激发态原子 A* 质量 (amu)
        real(dp)          :: m_b_amu          !< 靶原子 B 质量 (amu)
        real(dp)          :: mu_amu           !< 约化质量 (amu)
        real(dp)          :: mu_au            !< 约化质量 (a.u.)
        real(dp)          :: e_exc_ev         !< A* 激发内能 (eV)
        real(dp)          :: ip_target_ev     !< B 原子第一电离势 (eV)
        real(dp)          :: excess_energy_ev !< 渐近电离多余能量 E_exc - IP (eV)
        ! 入口通道中性相互作用势 V*(R) (Morse 势参数)
        real(dp)          :: v_star_de_ev     !< 阱深 D_e* (eV)
        real(dp)          :: v_star_re_ang    !< 平衡核间距 R_e* (Angstrom)
        real(dp)          :: v_star_beta_ang1 !< Morse 刚度参数 beta* (Angstrom^-1)
        ! 出口通道分子离子相互作用势 V+(R) (Morse 势参数，用于 AB+)
        real(dp)          :: v_plus_de_ev     !< 离子阱深 D_e+ (eV)
        real(dp)          :: v_plus_re_ang    !< 离子平衡核间距 R_e+ (Angstrom)
        real(dp)          :: v_plus_beta_ang1 !< 离子 Morse 刚度参数 beta+ (Angstrom^-1)
        ! 自电离跃迁宽度 Gamma(R) = A_gamma * exp(-alpha_gamma * R)
        real(dp)          :: gamma_a_au       !< 自电离宽度振幅 (a.u.)
        real(dp)          :: gamma_alpha_au   !< 自电离衰减常数 (a.u.^-1)
        real(dp)          :: spin_factor      !< 自旋统计因子 (0.0 < f_s <= 1.0)
    end type penning_system_t

    !> \brief 单碰撞参数轨迹计算结果
    type, public :: penning_trajectory_result_t
        real(dp) :: b_impact_ang  !< 碰撞参数 b (Angstrom)
        real(dp) :: r_turn_ang    !< 经典向径转折点 R_turn (Angstrom)
        real(dp) :: p_ion_tot     !< 总自电离概率 P_ion(b)
        real(dp) :: p_penning     !< 潘宁电离概率 P_PI(b) (产物为 A + B+ + e-)
        real(dp) :: p_associative !< 缔合电离概率 P_AI(b) (产物为 AB+ + e-)
    end type penning_trajectory_result_t

    !> \brief 能量扫描截面与速率常数结果结构体
    type, public :: penning_cross_section_result_t
        real(dp) :: e_coll_ev          !< 碰撞动能 E_coll (eV)
        real(dp) :: v_rel_au           !< 相对运动速度 (a.u.)
        real(dp) :: sigma_tot_ang2     !< 总电离截面 sigma_tot (Angstrom^2)
        real(dp) :: sigma_pi_ang2      !< 潘宁电离截面 sigma_PI (Angstrom^2)
        real(dp) :: sigma_ai_ang2      !< 缔合电离截面 sigma_AI (Angstrom^2)
        real(dp) :: ai_branching_ratio !< 缔合电离分支比 F_AI = sigma_AI / sigma_tot
        real(dp) :: rate_coeff_cm3_s   !< 双体总电离速率系数 k_tot (cm^3/s)
    end type penning_cross_section_result_t

    ! --------------------------------------------------------------------------
    ! 公共子程序与函数导出
    ! --------------------------------------------------------------------------
    public :: init_penning_system
    public :: init_penning_preset_he_star_ar
    public :: eval_penning_potentials
    public :: calc_penning_turning_point
    public :: calc_penning_trajectory_prob
    public :: calc_penning_cross_sections
    public :: calc_penning_electron_spectrum
    public :: calc_penning_thermal_rate
    public :: calc_ultracold_penning_complex_length

contains

    ! ==========================================================================
    ! 1. 潘宁与缔合电离体系初始化
    ! ==========================================================================
    subroutine init_penning_system(sys, name, m_a_amu, m_b_amu, e_exc_ev, &
                                   ip_target_ev, v_star_de_ev, v_star_re_ang, &
                                   v_star_beta_ang1, v_plus_de_ev, v_plus_re_ang, &
                                   v_plus_beta_ang1, gamma_a_au, gamma_alpha_au, &
                                   spin_factor, stat)
        type(penning_system_t), intent(out) :: sys
        character(len=*), intent(in)        :: name
        real(dp), intent(in)                :: m_a_amu
        real(dp), intent(in)                :: m_b_amu
        real(dp), intent(in)                :: e_exc_ev
        real(dp), intent(in)                :: ip_target_ev
        real(dp), intent(in)                :: v_star_de_ev
        real(dp), intent(in)                :: v_star_re_ang
        real(dp), intent(in)                :: v_star_beta_ang1
        real(dp), intent(in)                :: v_plus_de_ev
        real(dp), intent(in)                :: v_plus_re_ang
        real(dp), intent(in)                :: v_plus_beta_ang1
        real(dp), intent(in)                :: gamma_a_au
        real(dp), intent(in)                :: gamma_alpha_au
        real(dp), optional, intent(in)      :: spin_factor
        integer, optional, intent(out)      :: stat

        if (present(stat)) stat = 0
        if (m_a_amu <= 0.0_dp .or. m_b_amu <= 0.0_dp .or. e_exc_ev <= ip_target_ev) then
            if (present(stat)) stat = -1
            return
        end if

        sys%name             = trim(name)
        sys%m_a_amu          = m_a_amu
        sys%m_b_amu          = m_b_amu
        sys%mu_amu           = (m_a_amu * m_b_amu) / (m_a_amu + m_b_amu)
        sys%mu_au            = sys%mu_amu * AMU2AU
        sys%e_exc_ev         = e_exc_ev
        sys%ip_target_ev     = ip_target_ev
        sys%excess_energy_ev = e_exc_ev - ip_target_ev

        sys%v_star_de_ev     = v_star_de_ev
        sys%v_star_re_ang    = v_star_re_ang
        sys%v_star_beta_ang1 = v_star_beta_ang1

        sys%v_plus_de_ev     = v_plus_de_ev
        sys%v_plus_re_ang    = v_plus_re_ang
        sys%v_plus_beta_ang1 = v_plus_beta_ang1

        sys%gamma_a_au       = gamma_a_au
        sys%gamma_alpha_au   = gamma_alpha_au

        if (present(spin_factor)) then
            sys%spin_factor = max(0.0_dp, min(1.0_dp, spin_factor))
        else
            sys%spin_factor = 1.0_dp
        end if
    end subroutine init_penning_system

    ! ==========================================================================
    ! 2. 国际基准系统: He*(2^3S) + Ar -> He + Ar+ + e- (或 HeAr+ + e-)
    ! 实验与经典文献: H. Hotop & A. Niehaus, Z. Phys. 228, 68 (1969);
    !               P. E. Siska, Chem. Rev. 93, 18 (1993).
    ! ==========================================================================
    subroutine init_penning_preset_he_star_ar(sys, stat)
        type(penning_system_t), intent(out) :: sys
        integer, optional, intent(out)      :: stat

        ! He*(2^3S): E_exc = 19.82 eV, Ar: IP = 15.76 eV, Excess E = 4.06 eV
        ! V*(R): D_e = 5.2 meV = 0.0052 eV, R_e = 4.15 A, beta = 1.15 A^-1
        ! V+(R): HeAr+ D_e = 0.18 eV, R_e = 2.85 A, beta = 1.80 A^-1
        ! Gamma(R): A = 0.45 a.u., alpha = 1.10 a.u.^-1
        call init_penning_system(sys, name="He*(2^3S) + Ar", &
                                 m_a_amu=4.0026_dp, m_b_amu=39.948_dp, &
                                 e_exc_ev=19.82_dp, ip_target_ev=15.76_dp, &
                                 v_star_de_ev=0.0052_dp, v_star_re_ang=4.15_dp, &
                                 v_star_beta_ang1=1.15_dp, &
                                 v_plus_de_ev=0.180_dp, v_plus_re_ang=2.85_dp, &
                                 v_plus_beta_ang1=1.80_dp, &
                                 gamma_a_au=0.45_dp, gamma_alpha_au=1.10_dp, &
                                 spin_factor=1.0_dp, stat=stat)
    end subroutine init_penning_preset_he_star_ar

    ! ==========================================================================
    ! 3. 势能面与自电离跃迁宽度计算
    ! 入口态 Morse 势: V*(R) = D_e* * [ (1 - exp(-beta*(R - R_e)))^2 - 1 ]
    ! 离子态 Morse 势: V+(R) = D_e+ * [ (1 - exp(-beta+(R - R_e+)))^2 - 1 ]
    ! 宽度 Gamma(R) = A * exp(-alpha * R)
    ! ==========================================================================
    pure subroutine eval_penning_potentials(sys, r_ang, v_star_ev, v_plus_ev, gamma_ev)
        type(penning_system_t), intent(in) :: sys
        real(dp), intent(in)               :: r_ang
        real(dp), intent(out)              :: v_star_ev
        real(dp), intent(out)              :: v_plus_ev
        real(dp), intent(out)              :: gamma_ev

        real(dp) :: term_star, term_plus, r_au

        ! 入口通道势
        term_star = 1.0_dp - exp(-sys%v_star_beta_ang1 * (r_ang - sys%v_star_re_ang))
        v_star_ev = sys%v_star_de_ev * (term_star**2 - 1.0_dp)

        ! 离子出口通道势
        term_plus = 1.0_dp - exp(-sys%v_plus_beta_ang1 * (r_ang - sys%v_plus_re_ang))
        v_plus_ev = sys%v_plus_de_ev * (term_plus**2 - 1.0_dp)

        ! 自电离跃迁宽度 Gamma(R) (转换为 eV)
        r_au = r_ang * ANG2AU
        gamma_ev = sys%gamma_a_au * exp(-sys%gamma_alpha_au * r_au) * AU2EV
    end subroutine eval_penning_potentials

    ! ==========================================================================
    ! 4. 有效势能与经典向径转折点 R_turn(E, b)
    ! V_eff(R) = V*(R) + E_coll * (b / R)^2
    ! 满足 E_coll - V_eff(R_turn) = 0
    ! ==========================================================================
    subroutine calc_penning_turning_point(sys, e_coll_ev, b_ang, r_turn_ang, stat)
        type(penning_system_t), intent(in) :: sys
        real(dp), intent(in)               :: e_coll_ev
        real(dp), intent(in)               :: b_ang
        real(dp), intent(out)              :: r_turn_ang
        integer, optional, intent(out)     :: stat

        real(dp) :: r_lo, r_hi, r_mid, f_lo, f_hi, f_mid
        real(dp) :: v_star, v_plus, gamma
        integer  :: iter

        if (present(stat)) stat = 0
        r_turn_ang = b_ang

        if (e_coll_ev <= 0.0_dp .or. b_ang < 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        ! 二分法边界设定
        r_lo = 0.5_dp
        r_hi = max(b_ang + 5.0_dp, sys%v_star_re_ang + 6.0_dp)

        call eval_penning_potentials(sys, r_lo, v_star, v_plus, gamma)
        f_lo = e_coll_ev - (v_star + e_coll_ev * (b_ang / r_lo)**2)

        call eval_penning_potentials(sys, r_hi, v_star, v_plus, gamma)
        f_hi = e_coll_ev - (v_star + e_coll_ev * (b_ang / r_hi)**2)

        if (f_lo > 0.0_dp) then
            ! 入射能量极高，转折点小于 r_lo
            r_turn_ang = r_lo
            return
        end if

        ! 二分法寻根
        do iter = 1, 64
            r_mid = 0.5_dp * (r_lo + r_hi)
            call eval_penning_potentials(sys, r_mid, v_star, v_plus, gamma)
            f_mid = e_coll_ev - (v_star + e_coll_ev * (b_ang / r_mid)**2)

            if (abs(f_mid) < 1.0e-10_dp .or. (r_hi - r_lo) < 1.0e-8_dp) exit
            if (f_lo * f_mid <= 0.0_dp) then
                r_hi = r_mid
                f_hi = f_mid
            else
                r_lo = r_mid
                f_lo = f_mid
            end if
        end do

        r_turn_ang = r_mid
    end subroutine calc_penning_turning_point

    ! ==========================================================================
    ! 5. 半经典碰撞轨迹自电离概率与通道分支比计算
    ! 径向速度: v_R(R) = sqrt(2 * (E - V_eff*(R)) / mu)
    ! 存活几率: S(b) = exp( - 2 * \int_{R_0}^\infty \frac{\Gamma(R)}{\hbar v_R(R)} dR )
    ! 总自电离概率: P_ion(b) = 1 - S(b)
    ! 缔合电离判据: 若在 R 处电离释放电子后的核间相对动能 E - V+(R) < 0，
    ! 则产物陷入离子阱形成束缚二聚体离子 AB+ (Associative Ionization)。
    ! ==========================================================================
    subroutine calc_penning_trajectory_prob(sys, e_coll_ev, b_ang, res, stat)
        type(penning_system_t), intent(in)           :: sys
        real(dp), intent(in)                         :: e_coll_ev
        real(dp), intent(in)                         :: b_ang
        type(penning_trajectory_result_t), intent(out) :: res
        integer, optional, intent(out)               :: stat

        integer, parameter :: N_STEPS = 512
        real(dp) :: r_turn, r_max, dr, r_curr, v_eff_star, kin_e_ev
        real(dp) :: v_star, v_plus, gamma_ev, gamma_au, v_rad_au
        real(dp) :: integral_total, integral_ai, integrand, prob_survival
        integer  :: i

        if (present(stat)) stat = 0
        res%b_impact_ang  = b_ang
        res%r_turn_ang    = 0.0_dp
        res%p_ion_tot     = 0.0_dp
        res%p_penning     = 0.0_dp
        res%p_associative = 0.0_dp

        call calc_penning_turning_point(sys, e_coll_ev, b_ang, r_turn, stat)
        res%r_turn_ang = r_turn

        ! 积分上限 R_max (Gamma(R) 在该距离处已指数衰减至机器精度下)
        r_max = r_turn + 8.0_dp
        dr = (r_max - r_turn) / real(N_STEPS, dp)

        integral_total = 0.0_dp
        integral_ai    = 0.0_dp

        do i = 1, N_STEPS
            ! 采用梯形中点积分，避免转折点处 v_R = 0 奇异性
            r_curr = r_turn + (real(i, dp) - 0.5_dp) * dr
            call eval_penning_potentials(sys, r_curr, v_star, v_plus, gamma_ev)

            v_eff_star = v_star + e_coll_ev * (b_ang / r_curr)**2
            kin_e_ev = max(1.0e-9_dp, e_coll_ev - v_eff_star)

            ! 径向速度 (a.u.)
            v_rad_au = sqrt(2.0_dp * (kin_e_ev * EV2AU) / sys%mu_au)
            gamma_au = gamma_ev * EV2AU

            ! 积分核: 2 * Gamma(R) / (hbar * v_R) * dr (包含往返双程因子 2)
            integrand = 2.0_dp * (gamma_au / v_rad_au) * (dr * ANG2AU)
            integral_total = integral_total + integrand

            ! 缔合电离判据: 核动能加上相对势能差落入离子阱束缚态 E_coll - V+(R) < 0
            if (e_coll_ev + v_plus < 0.0_dp) then
                integral_ai = integral_ai + integrand
            end if
        end do

        ! 自旋因子调制
        integral_total = integral_total * sys%spin_factor
        integral_ai    = integral_ai * sys%spin_factor

        prob_survival = exp(-integral_total)
        res%p_ion_tot = 1.0_dp - prob_survival

        if (integral_total > 1.0e-30_dp) then
            res%p_associative = res%p_ion_tot * (integral_ai / integral_total)
        else
            res%p_associative = 0.0_dp
        end if
        res%p_penning = res%p_ion_tot - res%p_associative
    end subroutine calc_penning_trajectory_prob

    ! ==========================================================================
    ! 6. 碰撞电离截面与分支比计算
    ! \sigma(E) = 2\pi \int_0^{b_{max}} b P(b, E) db
    ! ==========================================================================
    subroutine calc_penning_cross_sections(sys, e_coll_ev, n_b, b_max_ang, res, stat)
        type(penning_system_t), intent(in)            :: sys
        real(dp), intent(in)                          :: e_coll_ev
        integer, intent(in)                           :: n_b
        real(dp), intent(in)                          :: b_max_ang
        type(penning_cross_section_result_t), intent(out) :: res
        integer, optional, intent(out)                :: stat

        real(dp) :: db, b_val, weight
        real(dp) :: sum_tot, sum_pi, sum_ai
        type(penning_trajectory_result_t) :: traj_res
        integer  :: i

        if (present(stat)) stat = 0
        res%e_coll_ev          = e_coll_ev
        res%v_rel_au           = 0.0_dp
        res%sigma_tot_ang2     = 0.0_dp
        res%sigma_pi_ang2      = 0.0_dp
        res%sigma_ai_ang2      = 0.0_dp
        res%ai_branching_ratio = 0.0_dp
        res%rate_coeff_cm3_s   = 0.0_dp

        if (e_coll_ev <= 0.0_dp .or. n_b <= 0 .or. b_max_ang <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        res%v_rel_au = sqrt(2.0_dp * (e_coll_ev * EV2AU) / sys%mu_au)

        db = b_max_ang / real(n_b, dp)
        sum_tot = 0.0_dp
        sum_pi  = 0.0_dp
        sum_ai  = 0.0_dp

        do i = 1, n_b
            b_val = (real(i, dp) - 0.5_dp) * db
            call calc_penning_trajectory_prob(sys, e_coll_ev, b_val, traj_res)

            weight = b_val * db
            sum_tot = sum_tot + traj_res%p_ion_tot * weight
            sum_pi  = sum_pi  + traj_res%p_penning * weight
            sum_ai  = sum_ai  + traj_res%p_associative * weight
        end do

        res%sigma_tot_ang2 = TWOPI * sum_tot
        res%sigma_pi_ang2  = TWOPI * sum_pi
        res%sigma_ai_ang2  = TWOPI * sum_ai

        if (res%sigma_tot_ang2 > 1.0e-30_dp) then
            res%ai_branching_ratio = res%sigma_ai_ang2 / res%sigma_tot_ang2
        else
            res%ai_branching_ratio = 0.0_dp
        end if

        ! 速率常数 k = v_rel * sigma_tot (转化为 cm^3/s: 1 a.u.(k) = 6.126048e-9 cm^3/s)
        ! sigma (a.u.) = sigma (Angstrom^2) * ANG2BOHR^2
        res%rate_coeff_cm3_s = res%v_rel_au * (res%sigma_tot_ang2 * (ANG2AU**2)) * 6.126048e-9_dp
    end subroutine calc_penning_cross_sections

    ! ==========================================================================
    ! 7. 潘宁电离电子能谱 (Penning Ionization Electron Spectroscopy - PIES)
    ! 出射电子动能: E_e(R) = V*(R) - V+(R) + E_excess
    ! 谱强 P(E_e) \propto \Gamma(R_c) / |d(V* - V+)/dR|_{R_c}
    ! ==========================================================================
    subroutine calc_penning_electron_spectrum(sys, e_coll_ev, n_e, e_min_ev, &
                                              e_max_ev, e_grid_ev, pies_intensity, stat)
        type(penning_system_t), intent(in) :: sys
        real(dp), intent(in)               :: e_coll_ev
        integer, intent(in)                :: n_e
        real(dp), intent(in)               :: e_min_ev
        real(dp), intent(in)               :: e_max_ev
        real(dp), intent(out)              :: e_grid_ev(n_e)
        real(dp), intent(out)              :: pies_intensity(n_e)
        integer, optional, intent(out)     :: stat

        integer, parameter :: N_R = 1000
        real(dp) :: r_min, r_max, dr, r, v_star, v_plus, gamma_ev
        real(dp) :: e_elec, weight, de, diff, delta_grid
        integer  :: ir, ie

        if (present(stat)) stat = 0
        pies_intensity = 0.0_dp

        if (n_e <= 1 .or. e_max_ev <= e_min_ev) then
            if (present(stat)) stat = -1
            return
        end if

        de = (e_max_ev - e_min_ev) / real(n_e - 1, dp)
        do ie = 1, n_e
            e_grid_ev(ie) = e_min_ev + real(ie - 1, dp) * de
        end do

        call calc_penning_turning_point(sys, max(1.0e-4_dp, e_coll_ev), 0.0_dp, r_min)
        r_max = max(r_min + 5.0_dp, 8.0_dp)
        dr = (r_max - r_min) / real(N_R, dp)

        ! 沿核间距累积局部电子发射展宽
        do ir = 1, N_R
            r = r_min + (real(ir, dp) - 0.5_dp) * dr
            call eval_penning_potentials(sys, r, v_star, v_plus, gamma_ev)

            ! 瞬时电离电子能量
            e_elec = sys%excess_energy_ev + (v_star - v_plus)
            ! 跃迁权重 (指数自电离宽度 * 空间体积元)
            weight = gamma_ev * (r**2) * dr

            ! 高斯仪器展宽卷积至能谱网格 (FWHM ~ 0.05 eV)
            delta_grid = 0.03_dp
            do ie = 1, n_e
                diff = (e_grid_ev(ie) - e_elec) / delta_grid
                if (abs(diff) < 4.0_dp) then
                    pies_intensity(ie) = pies_intensity(ie) + weight * exp(-0.5_dp * diff**2)
                end if
            end do
        end do

        ! 面积归一化
        weight = sum(pies_intensity) * de
        if (weight > 1.0e-30_dp) then
            pies_intensity = pies_intensity / weight
        end if
    end subroutine calc_penning_electron_spectrum

    ! ==========================================================================
    ! 8. 麦克斯韦玻尔兹曼热平衡电离速率常数 k(T)
    ! ==========================================================================
    subroutine calc_penning_thermal_rate(sys, temp_k, n_e, k_tot_cm3_s, &
                                         k_pi_cm3_s, k_ai_cm3_s, stat)
        type(penning_system_t), intent(in) :: sys
        real(dp), intent(in)               :: temp_k
        integer, intent(in)                :: n_e
        real(dp), intent(out)              :: k_tot_cm3_s
        real(dp), intent(out)              :: k_pi_cm3_s
        real(dp), intent(out)              :: k_ai_cm3_s
        integer, optional, intent(out)     :: stat

        real(dp) :: kb_t_ev, e_max_ev, de, e_curr, weight, sum_w
        real(dp) :: sum_tot, sum_pi, sum_ai
        type(penning_cross_section_result_t) :: cs_res
        integer  :: i

        if (present(stat)) stat = 0
        k_tot_cm3_s = 0.0_dp
        k_pi_cm3_s  = 0.0_dp
        k_ai_cm3_s  = 0.0_dp

        if (temp_k <= 0.0_dp .or. n_e <= 0) then
            if (present(stat)) stat = -1
            return
        end if

        ! 1 K = 8.617333262e-5 eV
        kb_t_ev = temp_k * 8.617333262e-5_dp
        e_max_ev = 8.0_dp * kb_t_ev
        de = e_max_ev / real(n_e, dp)

        sum_w   = 0.0_dp
        sum_tot = 0.0_dp
        sum_pi  = 0.0_dp
        sum_ai  = 0.0_dp

        do i = 1, n_e
            e_curr = (real(i, dp) - 0.5_dp) * de
            call calc_penning_cross_sections(sys, e_curr, n_b=64, b_max_ang=10.0_dp, res=cs_res)

            ! 3D 相对运动麦克斯韦速率加权: P(E) dE \propto E exp(-E / k_B T) dE
            weight = e_curr * exp(-e_curr / kb_t_ev) * de
            sum_w   = sum_w   + weight
            sum_tot = sum_tot + cs_res%rate_coeff_cm3_s * weight
            sum_pi  = sum_pi  + (cs_res%rate_coeff_cm3_s * (1.0_dp - cs_res%ai_branching_ratio)) * weight
            sum_ai  = sum_ai  + (cs_res%rate_coeff_cm3_s * cs_res%ai_branching_ratio) * weight
        end do

        if (sum_w > 1.0e-30_dp) then
            k_tot_cm3_s = sum_tot / sum_w
            k_pi_cm3_s  = sum_pi  / sum_w
            k_ai_cm3_s  = sum_ai  / sum_w
        end if
    end subroutine calc_penning_thermal_rate

    ! ==========================================================================
    ! 9. 超冷极限潘宁电离与 Wigner 阈值复散射长度
    ! s 波非弹性反应散射长度 a = \alpha - i \beta (\beta > 0)
    ! 超冷极限非弹性/电离损耗速率系数 K_loss = (4 \pi \hbar / \mu) * \beta
    ! ==========================================================================
    subroutine calc_ultracold_penning_complex_length(sys, a_scat_real_ang, &
                                                     beta_loss_ang, k_loss_cm3_s, stat)
        type(penning_system_t), intent(in) :: sys
        real(dp), intent(out)              :: a_scat_real_ang
        real(dp), intent(out)              :: beta_loss_ang
        real(dp), intent(out)              :: k_loss_cm3_s
        integer, optional, intent(out)     :: stat

        real(dp) :: c6_au, r_vdw_au, beta_au, k_loss_au

        if (present(stat)) stat = 0
        a_scat_real_ang = 0.0_dp
        beta_loss_ang   = 0.0_dp
        k_loss_cm3_s    = 0.0_dp

        ! 范德华长程色散长度 R_vdw = 0.5 * (2 * mu * C6 / hbar^2)^{1/4}
        ! 典型亚稳态稀有气体 C6 ~ 100 a.u.
        c6_au = 120.0_dp
        r_vdw_au = 0.5_dp * (2.0_dp * sys%mu_au * c6_au)**0.25_dp

        ! 弱自电离与强吸收模型
        a_scat_real_ang = r_vdw_au * AU2ANG
        ! 吸收虚部与自电离跃迁关联
        beta_au = sys%spin_factor * 0.28_dp * r_vdw_au
        beta_loss_ang = beta_au * AU2ANG

        ! 超冷损耗速率常数 K_loss = (4 * PI / mu) * beta (a.u.)
        k_loss_au = (4.0_dp * PI / sys%mu_au) * beta_au
        k_loss_cm3_s = k_loss_au * 6.126048e-9_dp
    end subroutine calc_ultracold_penning_complex_length

end module mod_penning_associative_ionization
