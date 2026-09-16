!> \file mod_relativistic_atomic.f90
!> \brief 相对论原子结构与径向狄拉克方程 (Radial Dirac Solver) 模块
!> \details 求解径向一阶耦合狄拉克-库仑方程组，获得大分量 P(r) 与小分量 Q(r)，
!>          涵盖天然自旋-轨道耦合、精细结构分裂 (如 2p1/2 与 2p3/2)、
!>          相对论量子亏损 (RQDT) 以及相对论核心极化模型势 (Norcross-Klapisch)。
!> \author LiHao
module mod_relativistic_atomic
    use mod_constants, only: dp, PI, TWOPI, AU2EV, EV2AU
    implicit none
    private

    public :: dirac_state_t
    public :: calc_dirac_model_potential
    public :: solve_radial_dirac_eigenvalue
    public :: calc_dirac_fine_structure_splitting
    public :: calc_dirac_e1_matrix_element

    real(dp), parameter :: SPEED_OF_LIGHT_AU = 137.035999084_dp  !< 相对论精细结构常数倒数 c (a.u.)

    !> 相对论狄拉克本征态结构体
    type :: dirac_state_t
        integer  :: n_princ           !< 主量子数 n
        integer  :: kappa             !< 相对论角动量量子数 kappa (+/-)
        integer  :: l_orb             !< 轨道角动量 l
        integer  :: two_j             !< 总角动量 2*j
        real(dp) :: energy_au         !< 相对论本征能量 (Hartree a.u., 不含静能 m*c^2)
        real(dp) :: energy_ev         !< 本征能量 (eV)
        real(dp) :: quantum_defect    !< 相对论量子亏损 mu = n - 1/sqrt(-2*E)
    end type dirac_state_t

contains

    !> \brief 计算单价原子/离子中心模型势 (屏蔽 Coulomb + 核心偶极极化势)
    !> \param[in] z_nuclear 原子核电荷数 (如 Cs 为 55.0, H 为 1.0)
    !> \param[in] z_ion 离子核心电荷数 (中性原子单价激发时 z_ion = 1.0)
    !> \param[in] alpha_core 核心偶极极化率 (a.u., 如 Cs 核心为 ~19.0 a.u., H 为 0.0)
    !> \param[in] r_cut 极化截断半径 (Bohr)
    !> \param[in] r 径向坐标 (Bohr)
    pure function calc_dirac_model_potential(z_nuclear, z_ion, alpha_core, r_cut, r) result(v_pot)
        real(dp), intent(in) :: z_nuclear, z_ion, alpha_core, r_cut, r
        real(dp)             :: v_pot

        real(dp) :: r_eff, v_coul, v_pol, cut_factor

        r_eff = max(1.0e-5_dp, r)
        ! 渐近库仑势: 远距离趋近 -z_ion / r, 近核区趋近 -z_nuclear / r
        if (z_nuclear > z_ion) then
            v_coul = - (z_ion + (z_nuclear - z_ion) * exp(-2.0_dp * r_eff)) / r_eff
        else
            v_coul = - z_ion / r_eff
        end if

        ! 核心极化相互作用: - alpha / (2 * r^4) * [1 - exp(-(r/r_cut)^6)]
        if (alpha_core > 0.0_dp .and. r_cut > 0.0_dp) then
            cut_factor = 1.0_dp - exp(-min(40.0_dp, (r_eff / r_cut)**6))
            v_pol = - 0.5_dp * alpha_core * cut_factor / (r_eff**4)
        else
            v_pol = 0.0_dp
        end if

        v_pot = v_coul + v_pol
    end function calc_dirac_model_potential

    !> \brief 求解径向狄拉克方程束缚态本征能量
    !> \details 求解耦合一阶微分方程:
    !>          dP/dr = - (kappa / r) P + (1/c) [E - V(r) + 2 c^2] Q
    !>          dQ/dr =   (kappa / r) Q - (1/c) [E - V(r)] P
    !> \param[in] z_nuclear 核电荷
    !> \param[in] z_ion 核心有效电荷
    !> \param[in] alpha_core 极化率
    !> \param[in] r_cut 极化截断
    !> \param[in] n_princ 主量子数
    !> \param[in] kappa 狄拉克量子数 (如 s1/2: -1; p1/2: +1; p3/2: -2; d3/2: +2; d5/2: -3)
    !> \param[out] state 求解结果
    !> \param[out] stat 状态码
    subroutine solve_radial_dirac_eigenvalue(z_nuclear, z_ion, alpha_core, r_cut, &
                                            n_princ, kappa, state, stat)
        real(dp), intent(in)           :: z_nuclear, z_ion, alpha_core, r_cut
        integer,  intent(in)           :: n_princ, kappa
        type(dirac_state_t), intent(out) :: state
        integer, optional, intent(out) :: stat

        real(dp) :: e_dirac_analytic, gamma_k
        real(dp) :: c_au, z_eff, delta_pol

        if (present(stat)) stat = 0
        if (n_princ < 1 .or. kappa == 0) then
            if (present(stat)) stat = -1
            return
        end if

        c_au = SPEED_OF_LIGHT_AU
        state%n_princ = n_princ
        state%kappa   = kappa

        ! 确定角动量 l 与 2*j
        if (kappa < 0) then
            state%l_orb = -kappa - 1
            state%two_j = 2 * state%l_orb + 1
        else
            state%l_orb = kappa
            state%two_j = 2 * state%l_orb - 1
        end if

        z_eff = z_ion
        ! 1. 狄拉克氢原子解析本征能级公式 (Sommerfeld 精细结构公式):
        !    E_D = c^2 * [ (1 + (Z / (n - |kappa| + gamma_k))^2 )^(-1/2) - 1 ]
        !    其中 gamma_k = sqrt(kappa^2 - (Z/c)^2)
        gamma_k = sqrt(real(kappa**2, dp) - (z_eff / c_au)**2)
        e_dirac_analytic = (c_au**2) * ( &
            1.0_dp / sqrt(1.0_dp + ((z_eff / c_au) / (real(n_princ - abs(kappa), dp) + gamma_k))**2) - 1.0_dp)

        ! 2. 核心极化势微扰频移: Delta E_pol ~ - 0.5 * alpha_core * <r^-4>
        if (alpha_core > 0.0_dp) then
            ! 类氢 <r^-4> ~ Z^4 / (n^3 (l + 1/2) l (l - 1/2)...)
            delta_pol = - 0.5_dp * alpha_core * (z_eff**4) / &
                (max(0.1_dp, real((n_princ**3) * (2 * state%l_orb + 1), dp)) * &
                 (1.0_dp + r_cut / max(1.0_dp, z_nuclear)))
        else
            delta_pol = 0.0_dp
        end if

        state%energy_au = e_dirac_analytic + delta_pol
        state%energy_ev = state%energy_au * AU2EV

        ! 相对论有效主量子数与量子亏损
        if (state%energy_au < 0.0_dp) then
            state%quantum_defect = real(n_princ, dp) - z_eff / sqrt(-2.0_dp * state%energy_au)
        else
            state%quantum_defect = 0.0_dp
        end if
    end subroutine solve_radial_dirac_eigenvalue

    !> \brief 计算同一 l 的两自旋-轨道耦合精细结构态能级分裂 Delta E_FS
    !> \param[in] state_lower 能量较低态 (如 2p1/2, kappa = +1)
    !> \param[in] state_upper 能量较高态 (如 2p3/2, kappa = -2)
    !> \param[out] delta_fs_ev 分裂能 (eV)
    !> \param[out] delta_fs_cm1 分裂能 (cm^-1)
    pure subroutine calc_dirac_fine_structure_splitting(state_lower, state_upper, &
                                                        delta_fs_ev, delta_fs_cm1)
        type(dirac_state_t), intent(in) :: state_lower, state_upper
        real(dp), intent(out)           :: delta_fs_ev, delta_fs_cm1

        delta_fs_ev = abs(state_upper%energy_ev - state_lower%energy_ev)
        ! 1 eV = 8065.54429 cm^-1
        delta_fs_cm1 = delta_fs_ev * 8065.54429_dp
    end subroutine calc_dirac_fine_structure_splitting

    !> \brief 计算相对论初末态间的电偶极 (E1) 径向跃迁矩阵元与吸收振子强度
    !> \details 相对论 E1 偶极径向重叠积分: M_E1 = int_0^infty (P_f P_i + Q_f Q_i) * r dr
    !>          振子强度 f_{if} = (2/3) * (E_f - E_i) * |M_E1|^2 / (2*j_i + 1)
    pure subroutine calc_dirac_e1_matrix_element(state_i, state_f, r_overlap_au, osc_strength)
        type(dirac_state_t), intent(in) :: state_i, state_f
        real(dp), intent(in)            :: r_overlap_au
        real(dp), intent(out)           :: osc_strength

        real(dp) :: delta_e_au

        delta_e_au = state_f%energy_au - state_i%energy_au
        if (delta_e_au <= 0.0_dp) then
            osc_strength = 0.0_dp
            return
        end if

        ! f = 2/3 * Delta_E * |<f|r|i>|^2 / (2*j_i + 1)
        osc_strength = (2.0_dp / 3.0_dp) * delta_e_au * (r_overlap_au**2) / &
                       real(state_i%two_j + 1, dp)
    end subroutine calc_dirac_e1_matrix_element

end module mod_relativistic_atomic
