! ==============================================================================
! GeneralModule: mod_dipolar_scattering.f90
! ------------------------------------------------------------------------------
! 现代超冷量子散射理论：各向异性磁偶极与电偶极碰撞及自旋弛豫模块
! 权威理论参考文献 (详见 LITERATURE.md):
!   - H. T. C. Stoof et al., Phys. Rev. B 38, 4688 (1988) [自旋弛豫与偶极相互作用]
!   - A. J. Moerdijk & B. J. Verhaar, Phys. Rev. A 53, 4341 (1996) [冷原子磁阱碰撞损耗]
!   - F. H. Mies, C. J. Williams, P. S. Julienne, J. Res. NIST 101, 521 (1996)
!   - J. L. Bohn, M. Cavagnero, C. Ticknor, New J. Phys. 11, 055039 (2009) [电偶极散射]
!   - G. Quemener & J. L. Bohn, Phys. Rev. A 81, 022702 (2010) [外加电场下极性分子碰撞]
!   - K.-K. Ni et al., Science 322, 231 (2008) [超冷极性分子量子气体 KRb]
! 核心特性：
! 1. 电子磁偶极-偶极相互作用 (MDDI) 2阶球张量算符及轨道-自旋矩阵元解析求值
! 2. s-d 分波各向异性耦合 (L=0 <-> L=2) 与宇称/总投影严格守恒定则
! 3. 弱场寻优态冷原子在静磁场下的非弹性偶极自旋弛豫 (Dipolar Relaxation) 截面与热速率
! 4. 极性分子外加直流电场 (Stark 效应) 态混合与实验室系诱导电偶极矩 d_ind(E) 求解
! 5. 各向异性电偶极相互作用 (EDDI) 矩阵元、偶极长度尺度 a_d 与多分波密耦势矩阵装配
! ==============================================================================
module mod_dipolar_scattering
    use mod_constants, only: dp, PI, TWOPI, HBAR, AMU2AU, AU2CM, CM2AU
    use mod_special_functions, only: wigner_3j, clebsch_gordan, wigner_3j_half, &
                                     clebsch_gordan_half, wigner_6j_half, wigner_9j_half
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    ! --------------------------------------------------------------------------
    ! 基础微观物理常量 (原子单位 a.u.)
    ! --------------------------------------------------------------------------
    real(dp), parameter, public :: FINE_STRUCT_ALPHA = 1.0_dp / 137.035999084_dp !< 精细结构常数 alpha
    real(dp), parameter, public :: ELECTRON_GS       = 2.00231930436256_dp       !< 电子自旋 g 因子
    real(dp), parameter, public :: BOHR_MAGNETON_AU  = 0.5_dp                    !< 玻尔磁子 mu_B (a.u.)

    ! --------------------------------------------------------------------------
    ! 派生类型定义
    ! --------------------------------------------------------------------------
    !> \brief 极性双原子分子参数结构体
    type, public :: polar_molecule_t
        character(len=16) :: name          !< 分子名称 (如 "KRb", "LiCs", "NaK")
        real(dp)          :: mass_amu      !< 分子质量 (amu)
        real(dp)          :: b_rot_cm1     !< 刚体转动常数 B_rot (cm^-1)
        real(dp)          :: dipole_debye  !< 永久电偶极矩 d_0 (Debye)
    end type polar_molecule_t

    !> \brief 各向异性散射分波通道结构体 (包含空间轨道与总自旋自由度)
    type, public :: dipolar_channel_t
        integer  :: l       !< 轨道角动量量子数 L (0=s, 2=d, 4=g)
        integer  :: ml      !< 轨道磁量子数 M_L (-L..L)
        integer  :: s_tot   !< 双电子总自旋 S (0=单重态, 1=三重态)
        integer  :: ms_tot  !< 总自旋磁量子数 M_S (-S..S)
        real(dp) :: energy_thresh !< 渐近通道阈值能量 (a.u.)
    end type dipolar_channel_t

    !> \brief 偶极自旋弛豫计算结果结构体
    type, public :: dipolar_relaxation_result_t
        real(dp) :: b_field_gauss          !< 外加磁场 (Gauss)
        real(dp) :: e_incident             !< 入射碰撞动能 (a.u.)
        real(dp) :: e_released             !< 释放的 Zeeman 内能 (a.u.)
        real(dp) :: k_incident             !< 入射动量
        real(dp) :: k_exit                 !< 出射通道相对动量
        real(dp) :: cross_section_au       !< 弛豫截面 sigma_rel (a.u.)
        real(dp) :: rate_coeff_cm3_s       !< 双体速率系数 K_rel (cm^3/s)
    end type dipolar_relaxation_result_t

    ! --------------------------------------------------------------------------
    ! 公共子程序与函数导出
    ! --------------------------------------------------------------------------
    public :: calc_mddi_spatial_matrix_element
    public :: calc_mddi_spin_matrix_element
    public :: calc_mddi_total_matrix_element
    public :: calc_dipolar_relaxation_cross_section
    public :: calc_dipolar_relaxation_thermal_rate
    public :: calc_stark_induced_dipole
    public :: calc_eddi_matrix_element
    public :: calc_dipole_length_scale
    public :: build_dipolar_coupled_potential_matrix

contains

    ! ==========================================================================
    ! 1. 磁偶极-偶极相互作用 (MDDI) 空间部分张量矩阵元
    ! <L', M_L' | C_{2, q}(\hat{r}) | L, M_L>
    ! 其中 C_{2, q} = sqrt(4*pi/5) * Y_{2, q} 为归一化 Racah 2 阶球谐张量
    ! 满足严格宇称守恒 (L + L' + 2 为偶数) 与磁量子数守恒 (q = M_L' - M_L)
    ! ==========================================================================
    pure function calc_mddi_spatial_matrix_element(l_bra, ml_bra, l_ket, ml_ket, q) result(mat_elem)
        integer, intent(in) :: l_bra, ml_bra
        integer, intent(in) :: l_ket, ml_ket
        integer, intent(in) :: q
        real(dp)            :: mat_elem

        real(dp) :: w3j_zero, w3j_m, pref
        integer  :: phase

        mat_elem = 0.0_dp

        ! 磁量子数选择定则: q + M_L = M_L' => q = M_L' - M_L
        if (ml_bra - ml_ket /= q) return

        ! 宇称定则: L + L' 必须为偶数
        if (mod(l_bra + l_ket, 2) /= 0) return

        ! 三角定则: |L - 2| <= L' <= L + 2
        if (l_bra < abs(l_ket - 2) .or. l_bra > (l_ket + 2)) return
        if (abs(ml_bra) > l_bra .or. abs(ml_ket) > l_ket .or. abs(q) > 2) return

        ! (2L+1)(2L'+1) 前因子
        pref = sqrt(real((2 * l_bra + 1) * (2 * l_ket + 1), dp))

        ! Wigner 3j 系数 (L' 2 L; 0 0 0)
        w3j_zero = wigner_3j(l_bra, 2, l_ket, 0, 0, 0)
        if (abs(w3j_zero) < 1.0e-15_dp) return

        ! Wigner 3j 系数 (L' 2 L; -M_L' q M_L)
        w3j_m = wigner_3j(l_bra, 2, l_ket, -ml_bra, q, ml_ket)

        ! 相位因子 (-1)^{M_L'}
        if (mod(abs(ml_bra), 2) == 1) then
            phase = -1
        else
            phase = 1
        end if

        mat_elem = real(phase, dp) * pref * w3j_zero * w3j_m
    end function calc_mddi_spatial_matrix_element

    ! ==========================================================================
    ! 2. 磁偶极-偶极相互作用 (MDDI) 自旋部分 2 阶不可约张量矩阵元
    ! <S', M_S' | [s1 (x) s2]^{(2)}_q | S, M_S>
    ! 双电子自旋 s1 = s2 = 1/2.
    ! 由角动量耦合定则，该 2 阶张量仅在 S = S' = 1 (三重态) 之间存在非零矩阵元！
    ! 约化矩阵元为 <1 || [s1 (x) s2]^{(2)} || 1> = 3/2.
    ! ==========================================================================
    pure function calc_mddi_spin_matrix_element(s_bra, ms_bra, s_ket, ms_ket, q) result(mat_elem)
        integer, intent(in) :: s_bra, ms_bra
        integer, intent(in) :: s_ket, ms_ket
        integer, intent(in) :: q
        real(dp)            :: mat_elem

        real(dp) :: w3j_val, red_mat
        integer  :: phase

        mat_elem = 0.0_dp

        ! 仅在三重态之间耦合 (S = S' = 1)
        if (s_bra /= 1 .or. s_ket /= 1) return
        if (abs(ms_bra) > 1 .or. abs(ms_ket) > 1 .or. abs(q) > 2) return

        ! 投影定则: q = M_S' - M_S
        if (ms_bra - ms_ket /= q) return

        ! 约化矩阵元 <1 || [s1 x s2]^(2) || 1> = 3/2 = 1.5
        red_mat = 1.5_dp

        ! Wigner 3j 符号 (1 2 1; -M_S' q M_S)
        w3j_val = wigner_3j(1, 2, 1, -ms_bra, q, ms_ket)

        ! 相位 (-1)^{1 - M_S'}
        if (mod(abs(1 - ms_bra), 2) == 1) then
            phase = -1
        else
            phase = 1
        end if

        mat_elem = real(phase, dp) * w3j_val * red_mat
    end function calc_mddi_spin_matrix_element

    ! ==========================================================================
    ! 3. 磁偶极-偶极相互作用总角向各向异性耦合强度系数 (无量纲角向因子)
    ! V_{dd}(r) = (alpha^2 / r^3) * C_{angle}
    ! 其中 C_{angle} = -sqrt(6) * sum_q (-1)^q <L' M_L'|C_{2, -q}|L M_L> *
    !                                         <S' M_S'|[s1 x s2]^{(2)}_q|S M_S>
    ! 严格守恒总磁量子数: M_tot = M_L + M_S = M_L' + M_S'
    ! ==========================================================================
    pure function calc_mddi_total_matrix_element(s_bra, ms_bra, l_bra, ml_bra, &
                                                 s_ket, ms_ket, l_ket, ml_ket) result(c_angle)
        integer, intent(in) :: s_bra, ms_bra, l_bra, ml_bra
        integer, intent(in) :: s_ket, ms_ket, l_ket, ml_ket
        real(dp)            :: c_angle

        integer  :: q, q_spatial, phase_q
        real(dp) :: spin_part, spat_part

        c_angle = 0.0_dp

        ! 检验全局总磁量子数守恒 M_tot = M_L + M_S
        if ((ml_bra + ms_bra) /= (ml_ket + ms_ket)) return

        ! 确定张量分量 q = M_S' - M_S
        q = ms_bra - ms_ket
        if (abs(q) > 2) return

        ! 对应的空间球张量分量为 -q (使得 -q = M_L' - M_L)
        q_spatial = -q
        if (ml_bra - ml_ket /= q_spatial) return

        spin_part = calc_mddi_spin_matrix_element(s_bra, ms_bra, s_ket, ms_ket, q)
        if (abs(spin_part) < 1.0e-15_dp) return

        spat_part = calc_mddi_spatial_matrix_element(l_bra, ml_bra, l_ket, ml_ket, q_spatial)
        if (abs(spat_part) < 1.0e-15_dp) return

        ! 相位 (-1)^q
        if (mod(abs(q), 2) == 1) then
            phase_q = -1
        else
            phase_q = 1
        end if

        ! V_{dd}(r) 前置因子为 -sqrt(6)
        c_angle = -sqrt(6.0_dp) * real(phase_q, dp) * spat_part * spin_part
    end function calc_mddi_total_matrix_element

    ! ==========================================================================
    ! 4. 磁阱中弱场寻优态冷原子的非弹性偶极自旋弛豫 (Dipolar Relaxation) 截面
    ! 入射初态: |S=1, M_S=1, L=0, M_L=0> (s-波)
    ! 出射末态: |S=1, M_S'=0, L'=2, M_L'=1> 或 |S=1, M_S'=-1, L'=2, M_L'=2> (d-波)
    ! 释放的动能: Delta_E = Delta_M_S * g_s * mu_B * B
    ! 基于 Born 近似与短程截止截面公式 (Stoof et al. 1988, Moerdijk 1996)
    ! ==========================================================================
    subroutine calc_dipolar_relaxation_cross_section(mass_amu, b_field_gauss, e_incident_au, &
                                                     r_cutoff_au, res, stat)
        real(dp), intent(in)                           :: mass_amu
        real(dp), intent(in)                           :: b_field_gauss
        real(dp), intent(in)                           :: e_incident_au
        real(dp), intent(in)                           :: r_cutoff_au
        type(dipolar_relaxation_result_t), intent(out) :: res
        integer, optional, intent(out)                 :: stat

        real(dp) :: mu_mass, b_au, delta_zeeman, e_exit, ki, kf
        real(dp) :: c_angle_dm1, c_angle_dm2, coupling_sq
        real(dp) :: matrix_r_integral, v_rel

        if (present(stat)) stat = 0
        if (mass_amu <= 0.0_dp .or. b_field_gauss < 0.0_dp .or. e_incident_au <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        ! 质心约化质量 mu = m / 2
        mu_mass = 0.5_dp * mass_amu * AMU2AU

        ! 外磁场转换为原子单位 (1 Gauss = 4.254382e-10 a.u.)
        b_au = b_field_gauss * 4.254382e-10_dp

        ! 单自旋翻转释放的 Zeeman 能量 Delta_E = g_s * mu_B * B
        delta_zeeman = ELECTRON_GS * BOHR_MAGNETON_AU * b_au
        e_exit = e_incident_au + delta_zeeman

        ki = sqrt(2.0_dp * mu_mass * e_incident_au)
        kf = sqrt(2.0_dp * mu_mass * e_exit)

        ! 计算两路主要弛豫通道的角向耦合强度:
        ! 通道 A: |1, 1; 0, 0> -> |1, 0; 2, 1> (Delta_M_S = -1, Delta_M_L = +1)
        c_angle_dm1 = calc_mddi_total_matrix_element(1, 0, 2, 1, 1, 1, 0, 0)

        ! 通道 B: |1, 1; 0, 0> -> |1, -1; 2, 2> (Delta_M_S = -2, Delta_M_L = +2)
        c_angle_dm2 = calc_mddi_total_matrix_element(1, -1, 2, 2, 1, 1, 0, 0)

        ! 总角向矩阵元平方和
        coupling_sq = c_angle_dm1**2 + c_angle_dm2**2

        ! 径向积分估计 (从硬核截止半径 r_cutoff 到无穷远处积分 1/r^3)
        matrix_r_integral = (FINE_STRUCT_ALPHA**2) / (r_cutoff_au**2)

        res%b_field_gauss    = b_field_gauss
        res%e_incident       = e_incident_au
        res%e_released       = delta_zeeman
        res%k_incident       = ki
        res%k_exit           = kf

        ! 构造 Born 截面标度
        res%cross_section_au = (4.0_dp * PI / (ki**2)) * (kf / ki) * &
                               (mu_mass * matrix_r_integral)**2 * coupling_sq * 0.1_dp

        ! 相对速度 v_rel = ki / mu_mass (a.u.)
        v_rel = ki / mu_mass
        ! 速率系数 K_rel = v_rel * sigma_rel, 转换为 cm^3 / s
        res%rate_coeff_cm3_s = v_rel * res%cross_section_au * 6.126048e-9_dp
    end subroutine calc_dipolar_relaxation_cross_section

    ! ==========================================================================
    ! 5. 玻尔兹曼热系综平均自旋弛豫速率系数 <K_{rel}>(T)
    ! ==========================================================================
    subroutine calc_dipolar_relaxation_thermal_rate(mass_amu, b_field_gauss, temp_k, &
                                                    r_cutoff_au, k_thermal, stat)
        real(dp), intent(in)           :: mass_amu
        real(dp), intent(in)           :: b_field_gauss
        real(dp), intent(in)           :: temp_k
        real(dp), intent(in)           :: r_cutoff_au
        real(dp), intent(out)          :: k_thermal
        integer, optional, intent(out) :: stat

        integer, parameter :: N_PTS = 64
        integer  :: i
        real(dp) :: kb_t_au, e_max, de, e_curr, weight, sum_k, sum_w
        type(dipolar_relaxation_result_t) :: res

        if (present(stat)) stat = 0
        if (temp_k <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        kb_t_au = temp_k * 3.166811563e-6_dp
        e_max = 8.0_dp * kb_t_au
        de = e_max / real(N_PTS, dp)

        sum_k = 0.0_dp
        sum_w = 0.0_dp

        do i = 1, N_PTS
            e_curr = (real(i, dp) - 0.5_dp) * de
            call calc_dipolar_relaxation_cross_section(mass_amu, b_field_gauss, e_curr, &
                                                       r_cutoff_au, res)
            weight = sqrt(e_curr) * exp(-e_curr / kb_t_au) * de
            sum_k = sum_k + res%rate_coeff_cm3_s * weight
            sum_w = sum_w + weight
        end do

        if (sum_w > 1.0e-30_dp) then
            k_thermal = sum_k / sum_w
        else
            k_thermal = 0.0_dp
        end if
    end subroutine calc_dipolar_relaxation_thermal_rate

    ! ==========================================================================
    ! 6. 极性分子外加直流电场 (Stark 效应) 实验室系诱导电偶极矩 d_ind(E)
    ! ==========================================================================
    subroutine calc_stark_induced_dipole(mol, e_field_kv_cm, d_ind_debye, stat)
        type(polar_molecule_t), intent(in) :: mol
        real(dp), intent(in)               :: e_field_kv_cm
        real(dp), intent(out)              :: d_ind_debye
        integer, optional, intent(out)     :: stat

        integer, parameter :: N_BASIS = 4
        real(dp) :: h_mat(N_BASIS, N_BASIS), eig_vals(N_BASIS), eig_vecs(N_BASIS, N_BASIS)
        real(dp) :: b_rot_au, d0_au, e_field_au
        integer  :: j1, j2, st
        real(dp) :: cos_mat(N_BASIS, N_BASIS)

        if (present(stat)) stat = 0
        if (mol%b_rot_cm1 <= 0.0_dp .or. mol%dipole_debye <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        b_rot_au = mol%b_rot_cm1 * CM2AU
        d0_au = mol%dipole_debye * 0.393430307_dp
        e_field_au = e_field_kv_cm / 5.14220674763e6_dp

        h_mat = 0.0_dp
        cos_mat = 0.0_dp

        do j1 = 1, N_BASIS
            h_mat(j1, j1) = real((j1 - 1) * j1, dp) * b_rot_au
        end do

        do j1 = 2, N_BASIS
            j2 = j1 - 1
            cos_mat(j1, j2) = real(j2, dp) / sqrt(real((2 * j2 - 1) * (2 * j2 + 1), dp))
            cos_mat(j2, j1) = cos_mat(j1, j2)
            h_mat(j1, j2) = -d0_au * e_field_au * cos_mat(j1, j2)
            h_mat(j2, j1) = h_mat(j1, j2)
        end do

        call diag_symmetric_matrix(N_BASIS, h_mat, eig_vals, eig_vecs, st)
        if (st /= 0) then
            if (present(stat)) stat = st
            return
        end if

        d_ind_debye = 0.0_dp
        do j1 = 1, N_BASIS
            do j2 = 1, N_BASIS
                d_ind_debye = d_ind_debye + eig_vecs(j1, 1) * eig_vecs(j2, 1) * &
                              mol%dipole_debye * cos_mat(j1, j2)
            end do
        end do
    end subroutine calc_stark_induced_dipole

    ! ==========================================================================
    ! 7. 各向异性电偶极相互作用 (EDDI) 空间张量矩阵元
    ! ==========================================================================
    pure function calc_eddi_matrix_element(l_bra, ml_bra, l_ket, ml_ket, d_ind_au) result(v_mat)
        integer, intent(in) :: l_bra, ml_bra
        integer, intent(in) :: l_ket, ml_ket
        real(dp), intent(in):: d_ind_au
        real(dp)            :: v_mat

        real(dp) :: c20_elem

        v_mat = 0.0_dp
        if (ml_bra /= ml_ket) return

        c20_elem = calc_mddi_spatial_matrix_element(l_bra, ml_bra, l_ket, ml_ket, 0)
        v_mat = -2.0_dp * (d_ind_au**2) * c20_elem
    end function calc_eddi_matrix_element

    ! ==========================================================================
    ! 8. 极性分子电偶极特征长度尺度 (Dipole Length Scale) a_d
    ! ==========================================================================
    pure function calc_dipole_length_scale(mass_amu, d_ind_debye) result(a_d_au)
        real(dp), intent(in) :: mass_amu
        real(dp), intent(in) :: d_ind_debye
        real(dp)            :: a_d_au

        real(dp) :: mu_au, d_au

        mu_au = 0.5_dp * mass_amu * AMU2AU
        d_au  = d_ind_debye * 0.393430307_dp
        a_d_au = 0.5_dp * mu_au * (d_au**2)
    end function calc_dipole_length_scale

    ! ==========================================================================
    ! 9. 组装多分波偶极耦合势能矩阵 V_{ij}(r)
    ! ==========================================================================
    subroutine build_dipolar_coupled_potential_matrix(channels, n_channels, mass_amu, c6_au, &
                                                      d_ind_debye, r, v_mat, stat)
        type(dipolar_channel_t), intent(in) :: channels(:)
        integer, intent(in)                 :: n_channels
        real(dp), intent(in)                :: mass_amu
        real(dp), intent(in)                :: c6_au
        real(dp), intent(in)                :: d_ind_debye
        real(dp), intent(in)                :: r
        real(dp), intent(out)               :: v_mat(:, :)
        integer, optional, intent(out)      :: stat

        integer  :: i, j
        real(dp) :: mu_au, d_ind_au, r3, r6, v_cent, v_c6, v_dip

        if (present(stat)) stat = 0
        if (r <= 0.0_dp .or. n_channels < 1) then
            if (present(stat)) stat = -1
            return
        end if

        mu_au = 0.5_dp * mass_amu * AMU2AU
        d_ind_au = d_ind_debye * 0.393430307_dp
        r3 = r**3
        r6 = r**6

        v_mat(1:n_channels, 1:n_channels) = 0.0_dp

        do i = 1, n_channels
            do j = 1, n_channels
                if (i == j) then
                    v_cent = real(channels(i)%l * (channels(i)%l + 1), dp) / (2.0_dp * mu_au * r**2)
                    v_c6   = -c6_au / r6
                    v_mat(i, i) = channels(i)%energy_thresh + v_cent + v_c6
                end if

                v_dip = calc_eddi_matrix_element(channels(i)%l, channels(i)%ml, &
                                                 channels(j)%l, channels(j)%ml, d_ind_au) / r3
                v_mat(i, j) = v_mat(i, j) + v_dip
            end do
        end do
    end subroutine build_dipolar_coupled_potential_matrix

end module mod_dipolar_scattering
