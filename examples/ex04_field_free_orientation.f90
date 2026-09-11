!> \brief 示例 4: 刚性转子分子外场定向/排列动力学与玻尔兹曼热系综平均
!> \details 模拟超快脉冲触发 CO 分子转动波包形成，计算无场量子拍与定向度 <cos(theta)>(t)，演示多初态热加权平均。
program ex04_field_free_orientation
    use general_module
    implicit none

    integer, parameter :: j_max = 15
    integer, parameter :: n_states = j_max + 1
    integer, parameter :: nt = 500

    real(dp) :: b_rot, dipole, alpha_par, alpha_perp, delta_alpha
    real(dp) :: temp_k, t_rot, dt, t
    real(dp) :: w_boltz(0:j_max), z_part
    real(dp) :: cos_mat(0:j_max, 0:j_max)
    real(dp) :: cos2_mat(0:j_max, 0:j_max)
    complex(dp) :: c_coeff(0:j_max)
    complex(dp) :: c_all(nt, 0:j_max)
    real(dp) :: cos_series(nt, 0:j_max)
    real(dp) :: cos_avg(nt)

    type(pulse_config_t) :: p_kick
    integer :: j, jp, j0, it, file_unit
    real(dp) :: kick_phase

    print '(A)', "=========================================================="
    print '(A)', "  Example 04: Field-Free Orientation & Thermal Ensemble   "
    print '(A)', "=========================================================="

    ! 1. 分子参数设置 (CO 分子参数，原子单位 a.u.)
    b_rot = 1.931_dp * CM2AU             ! 转动常数 B_e ~ 1.93 cm^-1
    dipole = 0.112_dp * DEBYE2AU         ! 永久偶极矩 mu_0 ~ 0.11 Debye
    alpha_par = 15.4_dp                  ! 平行极化率 (a.u.)
    alpha_perp = 11.8_dp                 ! 垂直极化率 (a.u.)
    delta_alpha = alpha_par - alpha_perp ! 极化率各向异性

    t_rot = PI / b_rot                   ! 转动周期 T_rot = pi / B (a.u.)
    print '(A, F10.3, A)', " CO Rotational Period: ", t_rot * AU2PS, " ps"

    ! 2. 预计算转动跃迁矩阵元 (M=0 子空间)
    cos_mat = 0.0_dp
    cos2_mat = 0.0_dp
    do j = 0, j_max
        do jp = 0, j_max
            cos_mat(j, jp) = rot_matrix_cos_theta(j, jp, 0)
            cos2_mat(j, jp) = rot_matrix_cos2_theta(j, jp, 0)
        end do
    end do

    ! 3. 玻尔兹曼转动权重计算 (T = 30 K)
    temp_k = 30.0_dp
    call boltzmann_rotational_weights(temp_k, b_rot, j_max, w_boltz, z_part)
    print '(A, F6.1, A, F10.4)', " Boltzmann Partition Z at T = ", temp_k, " K: ", z_part

    ! 4. 设置飞秒 THz 触发定向脉冲
    p_kick%shape_type = PULSE_SIN2
    p_kick%field_peak = 0.005_dp          ! ~25 MV/cm 强 THz 场
    p_kick%duration = 200.0_dp * FS2AU   ! 半高宽 200 fs (小于转动周期)
    p_kick%t_center = 0.0_dp

    ! 5. 对每个初态 J0 模拟转动跃迁并在无场自由演化
    dt = t_rot / real(nt, dp)

    do j0 = 0, j_max
        ! 初始态 |J0>
        c_coeff = (0.0_dp, 0.0_dp)
        c_coeff(j0) = (1.0_dp, 0.0_dp)

        ! 脉冲作用瞬时冲量近似 (Impulsive Kick):
        ! delta_P = int mu * E(t) dt ~ mu * E_0 * duration
        kick_phase = dipole * p_kick%field_peak * p_kick%duration
        ! 激发形成相干叠加态 c(J)
        do j = 0, j_max
            if (abs(j - j0) == 1) then
                c_coeff(j) = -EYE * 0.5_dp * kick_phase * cos_mat(j, j0)
            end if
        end do
        ! 重新归一化
        c_coeff = c_coeff / sqrt(sum(abs(c_coeff)**2))

        ! 自由演化含时展开: c_J(t) = c_J(0) * exp(-i * E_J * t)
        do it = 1, nt
            t = real(it - 1, dp) * dt
            do j = 0, j_max
                c_all(it, j) = c_coeff(j) * exp(-EYE * b_rot * real(j * (j + 1), dp) * t)
            end do

            ! 计算期望值 <cos(theta)>(t) = sum_{j, jp} c_j* c_jp <j|cos|jp>
            cos_series(it, j0) = 0.0_dp
            do j = 0, j_max
                do jp = 0, j_max
                    cos_series(it, j0) = cos_series(it, j0) + &
                        real(conjg(c_all(it, j)) * c_all(it, jp), dp) * cos_mat(j, jp)
                end do
            end do
        end do
    end do

    ! 6. 热系综平均: <cos(theta)>_T = sum_J0 w_{J0} * <cos>_{J0}(t)
    call thermal_average_2d(cos_series, w_boltz, cos_avg)

    ! 7. 保存结果数据
    open(newunit=file_unit, file="field_free_alignment.dat", status="replace", action="write")
    write(file_unit, '(A)') "# Time(ps)  <cos(theta)>_T(30K)  <cos(theta)>_J0  <cos(theta)>_J1"
    do it = 1, nt
        t = real(it - 1, dp) * dt
        write(file_unit, '(4ES16.8)') t * AU2PS, cos_avg(it), cos_series(it, 0), cos_series(it, 1)
    end do
    close(file_unit)
    print '(A)', "Thermal orientation dynamics saved to: field_free_alignment.dat"

end program ex04_field_free_orientation
