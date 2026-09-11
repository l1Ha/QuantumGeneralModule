!> \brief 完整激光脉冲调控分子转振态分布测试与综合算例
!> \details 综合利用 Sinc-DVR、FGH 束缚态求解、转动常数积分、全转振偶极跃迁矩阵构建、
!>          超快红外激光脉冲合成与 TDSE 动力学推进，模拟激光脉冲驱动双原子分子
!>          从基态 |v=0, J=0> 到激发转振态 |v=1, J=1> 及梯级跃迁的相干布居调控全过程。
program test_laser_rovibrational_control
    use general_module
    implicit none

    ! ----------------------------------------------------
    ! 1. 物理模型与参数定义 (HF 双原子分子模型)
    ! ----------------------------------------------------
    integer, parameter :: v_max = 2        ! 振动能级截断 v = 0, 1, 2
    integer, parameter :: j_max = 3        ! 转动能级截断 J = 0, 1, 2, 3
    integer, parameter :: n_states = (v_max + 1) * (j_max + 1) ! 总转振态维度 = 12
    integer, parameter :: n_pts = 256      ! 空间网格点数
    integer, parameter :: nt = 4000        ! 时间演化步数

    type(dvr_1d_t) :: dvr
    real(dp) :: r_min, r_max, mu_mass
    real(dp) :: d_e, r_e, beta
    real(dp), allocatable :: v_morse(:), dip_r(:), inv_r2(:)
    real(dp), allocatable :: eig_vals(:), eig_vecs(:, :), chi_wavefuncs(:, :)
    integer :: stat, v, j, k, step, file_unit

    ! 转振基底参数
    real(dp) :: e_vib(0:v_max)
    real(dp) :: b_v(0:v_max)
    real(dp) :: dip_vib(0:v_max, 0:v_max)
    real(dp) :: fc_mat(v_max+1, v_max+1)
    real(dp) :: h_diag(n_states)
    real(dp) :: dip_mat(n_states, n_states)

    ! 激光脉冲与动力学演化变量
    type(pulse_config_t) :: ir_pulse
    real(dp) :: omega_res, pulse_fwhm, pulse_peak
    real(dp) :: t_start, t_end, dt, t, e_field, norm_tot
    complex(dp) :: c_state(n_states)
    real(dp) :: pop_matrix(0:v_max, 0:j_max)
    real(dp) :: pop_v0, pop_v1, pop_v2
    real(dp) :: time_history(nt), pop_history(nt, n_states)

    ! 测试断言变量
    integer :: n_tests = 0, n_passed = 0
    real(dp) :: tol = 1.0e-5_dp

    call print_banner("GeneralModule: Laser Rovibrational State Control Test", 68)

    ! ----------------------------------------------------
    ! 2. Sinc-DVR 求解双原子分子振动态波函数与本征能级
    ! ----------------------------------------------------
    d_e = 0.225_dp           ! 势阱深度 (Hartree)
    r_e = 1.733_dp           ! 平衡核间距 (Bohr)
    beta = 1.174_dp          ! Morse 刚度参数
    mu_mass = 1744.5_dp      ! 约化质量 (a.u.)

    r_min = 0.8_dp
    r_max = 5.0_dp
    call dvr_sinc_init(r_min, r_max, n_pts, mu_mass, dvr)

    allocate(v_morse(n_pts), dip_r(n_pts), inv_r2(n_pts))
    allocate(eig_vals(n_pts), eig_vecs(n_pts, n_pts))

    ! 构造 Morse 势能与偶极函数 mu(R) = mu0 + mu1*(R - Re)
    do k = 1, n_pts
        v_morse(k) = d_e * (1.0_dp - exp(-beta * (dvr%x(k) - r_e)))**2
        dip_r(k) = 0.716_dp + 0.40_dp * (dvr%x(k) - r_e)
    end do

    call fgh_solve_bound_states(dvr, v_morse, eig_vals, eig_vecs, stat)
    call assert_true("FGH bound states solved successfully", stat == 0)

    ! 截取前 (v_max+1) 个束缚态并正确归一化为空间连续波函数: chi(x) = eig_vecs / sqrt(dx)
    allocate(chi_wavefuncs(n_pts, 0:v_max))
    do v = 0, v_max
        e_vib(v) = eig_vals(v + 1)
        chi_wavefuncs(:, v) = eig_vecs(:, v + 1) / sqrt(dvr%dx)
    end do

    ! ----------------------------------------------------
    ! 3. 计算各振动态转动常数 B_v 与跃迁偶极矩 M(v, v')
    ! ----------------------------------------------------
    call calc_rotational_constants_bv(chi_wavefuncs, dvr%x, dvr%dx, mu_mass, b_v)
    call calc_vibrational_dipole_matrix(chi_wavefuncs, dip_r, dvr%dx, dip_vib)
    call calc_franck_condon_factors(chi_wavefuncs, chi_wavefuncs, dvr%dx, fc_mat)

    print '(A, 3F12.5, A)', " Vibrational Energies E_v (cm^-1):   ", (e_vib(v) * AU2CM, v = 0, v_max)
    print '(A, 3F12.5, A)', " Rotational Constants B_v (cm^-1):   ", (b_v(v) * AU2CM, v = 0, v_max)
    print '(A, F10.6, A)',  " v=0 -> v=1 Transition Dipole (a.u.):", dip_vib(0, 1), ""

    ! 验证 Franck-Condon 对角归一化
    call assert_close("FC diagonal normalization FC(0,0) = 1.0", 1.0_dp, fc_mat(1, 1), 1.0e-5_dp)
    call assert_close("FC diagonal normalization FC(1,1) = 1.0", 1.0_dp, fc_mat(2, 2), 1.0e-5_dp)

    ! ----------------------------------------------------
    ! 4. 构造全转振基底 |v, J> 哈密顿量与偶极矩阵
    ! ----------------------------------------------------
    call build_rovibrational_hamiltonian(v_max, j_max, e_vib, b_v, h_diag)
    call build_rovibrational_dipole_matrix(v_max, j_max, dip_vib, dip_mat)

    ! 验证转动选择定则: <0,0|mu|0,0> 必须严格为 0
    k = rovibrational_state_index(0, 0, j_max)
    call assert_close("Dipole selection rule <0,0|mu|0,0> = 0", 0.0_dp, dip_mat(k, k), 1.0e-12_dp)

    ! 验证 R 支跃迁矩阵元 |<0,0|mu|1,1>| > 0
    k = rovibrational_state_index(1, 1, j_max)
    call assert_true("Dipole matrix element |<0,0|mu|1,1>| > 0", abs(dip_mat(1, k)) > 0.0_dp)

    ! ----------------------------------------------------
    ! 5. 激光脉冲参数定制 (调谐至 v=0,J=0 -> v=1,J=1 R支跃迁)
    ! ----------------------------------------------------
    ! 跃迁能量: Delta E = E(v=1, J=1) - E(v=0, J=0)
    omega_res = h_diag(rovibrational_state_index(1, 1, j_max)) - &
                h_diag(rovibrational_state_index(0, 0, j_max))

    pulse_fwhm = 80.0_dp    ! 80 fs 脉冲宽度
    pulse_peak = 0.015_dp   ! 强红外脉冲峰值场强 ~8e12 W/cm^2

    ! 便捷脉冲构造函数
    call create_sin2_pulse(peak=pulse_peak, dur_fs=pulse_fwhm, &
                           freq_ev=omega_res * AU2EV, t_center_fs=120.0_dp, &
                           cfg=ir_pulse)

    print '(A)', "------------------------------------------------------------"
    print '(A, F10.2, A)', " Resonant IR Pulse Frequency:    ", omega_res * AU2CM, " cm^-1"
    print '(A, F10.2, A)', " Resonant Photon Energy:         ", omega_res * AU2EV, " eV"
    print '(A, F10.2, A)', " Laser Duration (FWHM):          ", pulse_fwhm, " fs"
    print '(A, F10.4, A)', " Laser Peak Electric Field:      ", pulse_peak, " a.u."
    print '(A)', "------------------------------------------------------------"

    ! ----------------------------------------------------
    ! 6. TDSE 时间演化推进 (初态处于 |v=0, J=0>)
    ! ----------------------------------------------------
    c_state = (0.0_dp, 0.0_dp)
    c_state(rovibrational_state_index(0, 0, j_max)) = (1.0_dp, 0.0_dp)

    t_start = 0.0_dp
    t_end = 250.0_dp * FS2AU
    dt = (t_end - t_start) / real(nt, dp)

    do step = 1, nt
        t = t_start + real(step - 1, dp) * dt
        e_field = pulse_electric_field(t, ir_pulse)

        ! 4阶 Runge-Kutta 推进转振态系数向量
        call rk4_rovibrational_step(c_state, h_diag, dip_mat, e_field, dt)

        ! 记录时序数据
        time_history(step) = t * AU2FS
        do k = 1, n_states
            pop_history(step, k) = abs(c_state(k))**2
        end do
    end do

    ! ----------------------------------------------------
    ! 7. 演化结果分析与转振分布统计
    ! ----------------------------------------------------
    norm_tot = sum(abs(c_state)**2)
    call assert_close("Total probability norm conservation = 1.0", 1.0_dp, norm_tot, 5.0e-5_dp)

    ! 统计末态振动态与转动能级二维分布
    pop_v0 = 0.0_dp
    pop_v1 = 0.0_dp
    pop_v2 = 0.0_dp
    do v = 0, v_max
        do j = 0, j_max
            k = rovibrational_state_index(v, j, j_max)
            pop_matrix(v, j) = abs(c_state(k))**2
            if (v == 0) pop_v0 = pop_v0 + pop_matrix(v, j)
            if (v == 1) pop_v1 = pop_v1 + pop_matrix(v, j)
            if (v == 2) pop_v2 = pop_v2 + pop_matrix(v, j)
        end do
    end do

    print '(A)', " Final Rovibrational State Populations P(v, J):"
    print '(A10, 4A12)', "v \ J", "J = 0", "J = 1", "J = 2", "J = 3"
    print '(A)', "------------------------------------------------------------"
    do v = 0, v_max
        print '(A6, I2, 4F12.5)', "  v = ", v, (pop_matrix(v, j), j = 0, j_max)
    end do
    print '(A)', "------------------------------------------------------------"
    print '(A, F10.4)', " Total v=0 Population: ", pop_v0
    print '(A, F10.4)', " Total v=1 Population: ", pop_v1
    print '(A, F10.4)', " Total v=2 Population: ", pop_v2
    print '(A)', "------------------------------------------------------------"

    ! 断言检验布居转移成功发生:
    ! 1. 目标态 |v=1, J=1> 布居显著激发 (> 0.20)
    call assert_true("Target rovibrational state |1,1> excited (> 20%)", pop_matrix(1, 1) > 0.20_dp)
    ! 2. 总 v=1 振动布居数显著升高
    call assert_true("Total v=1 population significantly populated (> 25%)", pop_v1 > 0.25_dp)

    ! ----------------------------------------------------
    ! 8. 保存时序数据到标准文件 (测试 I/O 工具模块)
    ! ----------------------------------------------------
    call save_data_table_2d("test_rovibrational_dynamics.dat", time_history, pop_history, &
                            header="Time(fs) vs Rovibrational State Populations P(v, J)")
    print '(A)', ">> Rovibrational dynamics history saved to: test_rovibrational_dynamics.dat"

    ! 最终测试汇总
    print '(A)', "============================================================"
    print '(A, I2, A, I2, A)', "Rovibrational Control Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: Laser rovibrational control test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> \brief 4 阶 Runge-Kutta 推进含外场耦合的转振展开系数向量:
    !>        i * dc/dt = H_diag * c - E(t) * Dip_mat * c
    subroutine rk4_rovibrational_step(c_vec, diag_e, dip_m, e_t, dt_step)
        complex(dp), intent(inout) :: c_vec(:)
        real(dp), intent(in) :: diag_e(:), dip_m(:, :), e_t, dt_step

        complex(dp) :: k1(size(c_vec)), k2(size(c_vec)), k3(size(c_vec)), k4(size(c_vec))
        complex(dp) :: c_tmp(size(c_vec))

        k1 = deriv_c(c_vec, diag_e, dip_m, e_t)
        c_tmp = c_vec + 0.5_dp * dt_step * k1
        k2 = deriv_c(c_tmp, diag_e, dip_m, e_t)
        c_tmp = c_vec + 0.5_dp * dt_step * k2
        k3 = deriv_c(c_tmp, diag_e, dip_m, e_t)
        c_tmp = c_vec + dt_step * k3
        k4 = deriv_c(c_tmp, diag_e, dip_m, e_t)

        c_vec = c_vec + (dt_step / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end subroutine rk4_rovibrational_step

    !> \brief 计算系数导数: dc/dt = -i * (H_0 - E(t)*mu) * c
    pure function deriv_c(c_in, diag_e, dip_m, e_t) result(dcdt)
        complex(dp), intent(in) :: c_in(:)
        real(dp), intent(in) :: diag_e(:), dip_m(:, :), e_t
        complex(dp) :: dcdt(size(c_in))

        integer :: i_idx, j_idx, dim_n
        complex(dp) :: h_c

        dim_n = size(c_in)
        do i_idx = 1, dim_n
            h_c = diag_e(i_idx) * c_in(i_idx)
            do j_idx = 1, dim_n
                if (dip_m(i_idx, j_idx) /= 0.0_dp) then
                    h_c = h_c - e_t * dip_m(i_idx, j_idx) * c_in(j_idx)
                end if
            end do
            dcdt(i_idx) = -EYE * h_c
        end do
    end function deriv_c

    subroutine assert_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        n_tests = n_tests + 1
        if (condition) then
            n_passed = n_passed + 1
            print '(A, A45, A)', " [PASS] ", name, ""
        else
            print '(A, A45, A)', " [FAIL] ", name, " (Condition was FALSE)"
        end if
    end subroutine assert_true

    subroutine assert_close(name, expected, actual, eps)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, eps
        real(dp) :: diff

        n_tests = n_tests + 1
        diff = abs(expected - actual)
        if (diff <= eps) then
            n_passed = n_passed + 1
            print '(A, A45, A)', " [PASS] ", name, ""
        else
            print '(A, A45, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_laser_rovibrational_control
