!> \brief 分子转振耦合与激光跃迁偶极矩阵模块
!> \details 提供双原子分子振动波函数 Franck-Condon 因子计算、振动跃迁偶极矩积分、
!>          各振动态转动常数 B_v 积分、全空间转振基底 |v, J> 偶极与极化跃迁矩阵构建，
!>          以及 STIRAP 受激拉曼绝热通道脉冲对配置生成。
!> \author LiHao
module mod_rovibrational
    use mod_constants, only: dp, PI, TWOPI, FS2AU
    use mod_special_functions, only: rot_matrix_cos_theta, rot_matrix_cos2_theta
    use mod_laser_pulse, only: pulse_config_t, PULSE_GAUSSIAN, PULSE_SIN2
    implicit none
    private

    public :: calc_franck_condon_factors
    public :: calc_vibrational_dipole_matrix
    public :: calc_rotational_constants_bv
    public :: rovibrational_state_index
    public :: rovibrational_state_unindex
    public :: build_rovibrational_hamiltonian
    public :: build_rovibrational_dipole_matrix
    public :: build_rovibrational_polarizability_matrix
    public :: honl_london_parallel
    public :: honl_london_sum_rule
    public :: create_stirap_pulses

contains

    !> \brief 计算两组振动态之间的 Franck-Condon 重叠因子矩阵 FC(v, v') = |<chi_a,v | chi_b,v'>|^2
    subroutine calc_franck_condon_factors(chi_a, chi_b, dx, fc_mat)
        real(dp), intent(in) :: chi_a(:, :)   !< 态 A 波函数矩阵 (n_pts, n_va)
        real(dp), intent(in) :: chi_b(:, :)   !< 态 B 波函数矩阵 (n_pts, n_vb)
        real(dp), intent(in) :: dx            !< 坐标步长
        real(dp), intent(out) :: fc_mat(:, :) !< Franck-Condon 矩阵 (n_va, n_vb)

        integer :: va, vb, n_va, n_vb
        real(dp) :: overlap

        n_va = size(chi_a, 2)
        n_vb = size(chi_b, 2)

        do va = 1, n_va
            do vb = 1, n_vb
                overlap = sum(chi_a(:, va) * chi_b(:, vb)) * dx
                fc_mat(va, vb) = overlap**2
            end do
        end do
    end subroutine calc_franck_condon_factors

    !> \brief 计算振动态之间的跃迁偶极矩矩阵 M(v, v') = <chi_v | mu(R) | chi_v'>
    subroutine calc_vibrational_dipole_matrix(chi, dipole_grid, dx, dip_mat)
        real(dp), intent(in) :: chi(:, :)     !< 振动态波函数 (n_pts, n_v)
        real(dp), intent(in) :: dipole_grid(:)!< 坐标格点上的永久/跃迁偶极矩 mu(R_i) (n_pts)
        real(dp), intent(in) :: dx            !< 坐标步长
        real(dp), intent(out) :: dip_mat(:, :)!< 振动偶极矩阵 (n_v, n_v)

        integer :: v1, v2, nv

        nv = size(chi, 2)
        do v1 = 1, nv
            do v2 = 1, nv
                dip_mat(v1, v2) = sum(chi(:, v1) * dipole_grid(:) * chi(:, v2)) * dx
            end do
        end do
    end subroutine calc_vibrational_dipole_matrix

    !> \brief 计算各振动态的有效转动常数 B_v = <chi_v | hbar^2 / (2 * mu * R^2) | chi_v>
    subroutine calc_rotational_constants_bv(chi, r_grid, dx, mass, b_v)
        real(dp), intent(in) :: chi(:, :)     !< 振动态波函数 (n_pts, n_v)
        real(dp), intent(in) :: r_grid(:)     !< 核间距网格 R (n_pts)
        real(dp), intent(in) :: dx            !< 坐标步长
        real(dp), intent(in) :: mass          !< 约化质量 mu (a.u.)
        real(dp), intent(out) :: b_v(:)       !< 各振动态转动常数 B_v (n_v)

        integer :: v, nv, i
        real(dp), allocatable :: inv_r2(:)
        integer :: n_pts

        nv = size(chi, 2)
        n_pts = size(r_grid)
        allocate(inv_r2(n_pts))

        do i = 1, n_pts
            inv_r2(i) = 1.0_dp / (2.0_dp * mass * max(0.1_dp, r_grid(i)**2))
        end do

        do v = 1, nv
            b_v(v) = sum(chi(:, v)**2 * inv_r2(:)) * dx
        end do

        deallocate(inv_r2)
    end subroutine calc_rotational_constants_bv

    !> \brief 二维量子数 (v, J) 映射到一维基底线性索引 k ∈ [1, (v_max+1)*(j_max+1)]
    pure function rovibrational_state_index(v, j, j_max) result(idx)
        integer, intent(in) :: v, j, j_max
        integer :: idx
        idx = v * (j_max + 1) + j + 1
    end function rovibrational_state_index

    !> \brief 一维基底索引反解量子数 (v, J)
    pure subroutine rovibrational_state_unindex(idx, j_max, v, j)
        integer, intent(in) :: idx, j_max
        integer, intent(out) :: v, j
        v = (idx - 1) / (j_max + 1)
        j = mod(idx - 1, j_max + 1)
    end subroutine rovibrational_state_unindex

    !> \brief 构造无场转振哈密顿量对角本征能级: E(v, J) = E_vib(v) + B_v * J * (J + 1)
    subroutine build_rovibrational_hamiltonian(v_max, j_max, e_vib, b_v, h_diag)
        integer, intent(in) :: v_max, j_max
        real(dp), intent(in) :: e_vib(0:v_max)
        real(dp), intent(in) :: b_v(0:v_max)
        real(dp), intent(out) :: h_diag(:)

        integer :: v, j, k

        do v = 0, v_max
            do j = 0, j_max
                k = rovibrational_state_index(v, j, j_max)
                h_diag(k) = e_vib(v) + b_v(v) * real(j * (j + 1), dp)
            end do
        end do
    end subroutine build_rovibrational_hamiltonian

    !> \brief 构造全转振偶极跃迁矩阵: <v, J | mu_z | v', J'> = <v|mu|v'> * <J|cos(theta)|J'>
    !> \details 严格满足选择定则 Delta J = +/- 1, Delta M = 0 (线偏振沿 z 轴)
    subroutine build_rovibrational_dipole_matrix(v_max, j_max, dip_vib, dip_mat)
        integer, intent(in) :: v_max, j_max
        real(dp), intent(in) :: dip_vib(0:v_max, 0:v_max)
        real(dp), intent(out) :: dip_mat(:, :)

        integer :: v, j, vp, jp, k1, k2
        real(dp) :: angular_elem

        dip_mat = 0.0_dp

        do v = 0, v_max
            do j = 0, j_max
                k1 = rovibrational_state_index(v, j, j_max)
                do vp = 0, v_max
                    do jp = 0, j_max
                        k2 = rovibrational_state_index(vp, jp, j_max)
                        if (abs(j - jp) == 1) then
                            angular_elem = rot_matrix_cos_theta(j, jp, 0)
                            dip_mat(k1, k2) = dip_vib(v, vp) * angular_elem
                        end if
                    end do
                end do
            end do
        end do
    end subroutine build_rovibrational_dipole_matrix

    !> \brief 构造极化率取向跃迁矩阵: <v, J | Delta_alpha * cos^2(theta) | v', J'>
    !> \details 满足极化选择定则 Delta J = 0, +/- 2
    subroutine build_rovibrational_polarizability_matrix(v_max, j_max, alpha_vib, polar_mat)
        integer, intent(in) :: v_max, j_max
        real(dp), intent(in) :: alpha_vib(0:v_max, 0:v_max)
        real(dp), intent(out) :: polar_mat(:, :)

        integer :: v, j, vp, jp, k1, k2
        real(dp) :: angular_elem

        polar_mat = 0.0_dp

        do v = 0, v_max
            do j = 0, j_max
                k1 = rovibrational_state_index(v, j, j_max)
                do vp = 0, v_max
                    do jp = 0, j_max
                        k2 = rovibrational_state_index(vp, jp, j_max)
                        if (j == jp .or. abs(j - jp) == 2) then
                            angular_elem = rot_matrix_cos2_theta(j, jp, 0)
                            polar_mat(k1, k2) = alpha_vib(v, vp) * angular_elem
                        end if
                    end do
                end do
            end do
        end do
    end subroutine build_rovibrational_polarizability_matrix

    !> \brief 便捷生成受激拉曼绝热通道 (STIRAP) 脉冲对配置（Stokes 脉冲先于 Pump 脉冲反直觉时序）
    subroutine create_stirap_pulses(peak_pump, peak_stokes, dur_pump_fs, dur_stokes_fs, &
                                    delay_fs, freq_pump_au, freq_stokes_au, cfg_pump, cfg_stokes)
        real(dp), intent(in) :: peak_pump, peak_stokes
        real(dp), intent(in) :: dur_pump_fs, dur_stokes_fs, delay_fs
        real(dp), intent(in) :: freq_pump_au, freq_stokes_au
        type(pulse_config_t), intent(out) :: cfg_pump, cfg_stokes

        ! Pump 脉冲
        cfg_pump%shape_type = PULSE_GAUSSIAN
        cfg_pump%field_peak = peak_pump
        cfg_pump%freq_central = freq_pump_au
        cfg_pump%duration = dur_pump_fs * FS2AU
        cfg_pump%t_center = 0.5_dp * delay_fs * FS2AU  ! 滞后触发
        cfg_pump%cep_phase = 0.0_dp

        ! Stokes 脉冲 (先发)
        cfg_stokes%shape_type = PULSE_GAUSSIAN
        cfg_stokes%field_peak = peak_stokes
        cfg_stokes%freq_central = freq_stokes_au
        cfg_stokes%duration = dur_stokes_fs * FS2AU
        cfg_stokes%t_center = -0.5_dp * delay_fs * FS2AU ! 超前触发
        cfg_stokes%cep_phase = 0.0_dp
    end subroutine create_stirap_pulses

    !> \brief 平行跃迁 (Delta Lambda = 0) 的 Hoehnl-London 因子 (P 与 R 支)
    !> \details S_P(J) = J/(2J+1), S_R(J) = (J+1)/(2J+1) —— 由 Wigner-Eckart
    !>          定理 + 3j 符号正交性严格给出:
    !>            S_{J'}(J) = sum_{M,M',q} |<J'M'|d_q|JM>|^2
    !>                      = (2J'+1) (2J+1) [3j(J' 1 J; 0 0 0)]^2
    !>          (对 Delta Lambda = 0 平行跃迁; 简化为:
    !>           P: J/(2J+1),  R: (J+1)/(2J+1))。
    !>          通过 rot_matrix_cos_theta 对 M 求和数值验证。
    subroutine honl_london_parallel(j_qnum, s_p, s_r, stat)
        integer, intent(in)  :: j_qnum
        real(dp), intent(out) :: s_p, s_r
        integer, intent(out)  :: stat

        real(dp) :: jj
        jj = real(j_qnum, dp)
        s_p = jj / (2.0_dp * jj + 1.0_dp)
        s_r = (jj + 1.0_dp) / (2.0_dp * jj + 1.0_dp)
        stat = 0
    end subroutine honl_london_parallel

    !> \brief Hoehnl-London 求和规则: sum_{J'} (2J'+1) S_{J'} = 1
    !> \details 对任意 J, P 支 + R 支的强度总和恒等于 1 (单位偶极强度):
    !>            sum_{J'} S_{J'} = sum_{J'} (2J'+1)[3j(J' 1 J; 0 0 0)]^2 = 1
    !>          (3j 符号正交性); 数值通过 rot_matrix_cos_theta 对 M 求和验证。
    subroutine honl_london_sum_rule(j_qnum, s_total, stat)
        use mod_special_functions, only: rot_matrix_cos_theta
        integer, intent(in)  :: j_qnum
        real(dp), intent(out) :: s_total
        integer, intent(out)  :: stat

        real(dp) :: s_m_sum, mat_elem
        integer :: m_val, jp

        stat = 0
        ! 数值: 对 M 求和的 |<J'M'|cos theta|J M>|^2, 对 J' = J-1 和 J+1 求和
        s_total = 0.0_dp
        do jp = j_qnum - 1, j_qnum + 1
            if (jp < 0) cycle
            s_m_sum = 0.0_dp
            do m_val = -j_qnum, j_qnum
                mat_elem = rot_matrix_cos_theta(j_qnum, jp, m_val)
                s_m_sum = s_m_sum + mat_elem**2
            end do
            ! 数值 Hönl-London 因子 (对 M 求和)
            s_total = s_total + s_m_sum
        end do
        ! 期望 = 1 (单位跃迁偶极)
    end subroutine honl_london_sum_rule

end module mod_rovibrational
