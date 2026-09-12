! ==============================================================================
! GeneralModule Unit Tests: Time-Independent Scattering (mod_ti_scattering)
! ==============================================================================
program test_ti_scattering
    use mod_constants, only: dp, PI
    use mod_ti_scattering
    implicit none

    integer :: n_pass = 0, n_total = 0
    integer :: istat
    real(dp) :: jl, nl, djl, dnl, wronskian
    real(dp) :: mass, r_well, v0_well, kappa_well, as_exact
    real(dp) :: as_num, as_logder
    real(dp), allocatable :: r_grid(:), v_pot(:)
    integer :: n_pts, i
    real(dp) :: energy, k_wave, k_in, delta_exact, delta_calc, k_mat
    complex(dp) :: s_mat, t_mat
    real(dp) :: delta_arr(0:3), sigma_part(0:3), sigma_tot, sigma_opt
    real(dp) :: k_list(4), as_ere, r0_ere
    real(dp) :: c6_au, a_bar, as_gf
    real(dp), allocatable :: e_grid(:), delta_scan(:)
    type(resonance_info_t) :: res_info
    real(dp), allocatable :: v11(:), v22(:), v12(:)
    complex(dp) :: s_2x2(2, 2)
    real(dp) :: inel_prob, s_unitary_norm

    print '(A)', "=================================================="
    print '(A)', " GeneralModule Unit Tests: TI Scattering         "
    print '(A)', "=================================================="

    ! --------------------------------------------------------------------------
    ! 测试 1: Riccati-Bessel/Neumann 函数及 Wronskian 恒等式 W[j_hat, n_hat] = 1.0
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call riccati_bessel_neumann(0, 1.5_dp, jl, nl, djl, dnl)
    wronskian = jl * dnl - djl * nl
    if (abs(wronskian - 1.0_dp) < 1.0e-12_dp .and. abs(jl - sin(1.5_dp)) < 1.0e-12_dp) then
        print '(A, F10.7)', " [PASS] Riccati-Bessel Wronskian = ", wronskian
        n_pass = n_pass + 1
    else
        print '(A, F10.7)', " [FAIL] Riccati-Bessel Wronskian = ", wronskian
    end if

    ! --------------------------------------------------------------------------
    ! 测试 2: 吸引方势阱零能散射长度 Numerov 解与解析解对比
    ! V(r) = -V0 (r < R), 0 (r >= R)
    ! 解析解: a_s = R * (1 - tan(kappa*R)/(kappa*R)), kappa = sqrt(2*mu*V0)
    ! --------------------------------------------------------------------------
    mass = 1.0_dp
    r_well = 2.0_dp
    v0_well = 1.0_dp
    kappa_well = sqrt(2.0_dp * mass * v0_well)
    as_exact = r_well * (1.0_dp - tan(kappa_well * r_well) / (kappa_well * r_well))

    n_pts = 1000
    allocate(r_grid(n_pts), v_pot(n_pts))
    do i = 1, n_pts
        r_grid(i) = 0.01_dp + real(i - 1, dp) * (10.0_dp - 0.01_dp) / real(n_pts - 1, dp)
        if (r_grid(i) <= r_well) then
            v_pot(i) = -v0_well
        else
            v_pot(i) = 0.0_dp
        end if
    end do

    n_total = n_total + 1
    call calc_scattering_length_numerov(r_grid, v_pot, mass, as_num, stat=istat)
    if (abs(as_num - as_exact) / abs(as_exact) < 0.005_dp) then
        print '(A, F9.5, A, F9.5)', " [PASS] Numerov a_s = ", as_num, " (Exact: ", as_exact, ")"
        n_pass = n_pass + 1
    else
        print '(A, F9.5, A, F9.5)', " [FAIL] Numerov a_s = ", as_num, " (Exact: ", as_exact, ")"
    end if

    ! --------------------------------------------------------------------------
    ! 测试 3: 吸引方势阱 Johnson Log-Derivative 比值法求 a_s
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call calc_scattering_length_logder(r_grid, v_pot, mass, as_logder, stat=istat)
    if (abs(as_logder - as_exact) / abs(as_exact) < 0.005_dp) then
        print '(A, F9.5, A, F9.5)', " [PASS] Log-Der a_s = ", as_logder, " (Exact: ", as_exact, ")"
        n_pass = n_pass + 1
    else
        print '(A, F9.5, A, F9.5)', " [FAIL] Log-Der a_s = ", as_logder, " (Exact: ", as_exact, ")"
    end if

    ! --------------------------------------------------------------------------
    ! 测试 4: 有限正能量分波相移与 S-矩阵幺正性 |S_l| = 1.0
    ! 解析: tan(k*R + delta) = (k/K) * tan(K*R), K = sqrt(2*mu*(E+V0))
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    energy = 0.05_dp
    k_wave = sqrt(2.0_dp * mass * energy)
    k_in   = sqrt(2.0_dp * mass * (energy + v0_well))
    delta_exact = atan2(k_wave * tan(k_in * r_well), k_in) - k_wave * r_well

    call calc_phase_shift_single_l(r_grid, v_pot, mass, energy, 0, delta_calc, k_mat, s_mat, t_mat, istat)
    if (abs(abs(s_mat) - 1.0_dp) < 1.0e-5_dp) then
        print '(A, F9.6)', " [PASS] S-matrix unitarity |S_0| = ", abs(s_mat)
        n_pass = n_pass + 1
    else
        print '(A, F9.6)', " [FAIL] S-matrix unitarity |S_0| = ", abs(s_mat)
    end if

    ! --------------------------------------------------------------------------
    ! 测试 5: 光学定理自洽性 sigma_total = sigma_optical
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call calc_partial_wave_cross_sections(r_grid, v_pot, mass, energy, 3, &
                                          delta_arr, sigma_part, sigma_tot, istat)
    sigma_opt = optical_theorem_cross_section(k_wave, delta_arr, 3)
    if (abs(sigma_tot - sigma_opt) < 1.0e-8_dp .and. sigma_tot > 0.0_dp) then
        print '(A, F10.5, A, F10.5)', " [PASS] Optical theorem: sigma_tot = ", sigma_tot, " == sigma_opt = ", sigma_opt
        n_pass = n_pass + 1
    else
        print '(A, F10.5, A, F10.5)', " [FAIL] Optical theorem: sigma_tot = ", sigma_tot, " != sigma_opt = ", sigma_opt
    end if

    ! --------------------------------------------------------------------------
    ! 测试 6: 有效力程展开 (ERE) 拟合提取 a_s
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    k_list = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp]
    call fit_effective_range_expansion(r_grid, v_pot, mass, k_list, 4, as_ere, r0_ere, istat)
    if (abs(as_ere - as_exact) / abs(as_exact) < 0.01_dp) then
        print '(A, F9.5, A, F9.5)', " [PASS] ERE fit a_s = ", as_ere, " (Exact: ", as_exact, ")"
        n_pass = n_pass + 1
    else
        print '(A, F9.5, A, F9.5)', " [FAIL] ERE fit a_s = ", as_ere, " (Exact: ", as_exact, ")"
    end if

    ! --------------------------------------------------------------------------
    ! 测试 7: 范德华平均长度 a_bar 与 Gribakin-Flambaum 公式
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    c6_au = 50.0_dp
    a_bar = van_der_waals_mean_length(mass, c6_au)
    as_gf = gribakin_flambaum_length(mass, c6_au, phase_phi=PI/8.0_dp) ! phi = pi/8 时 tan(0)=0 => as = a_bar
    if (abs(as_gf - a_bar) < 1.0e-12_dp .and. a_bar > 0.0_dp) then
        print '(A, F9.5, A, F9.5)', " [PASS] Van der Waals a_bar = ", a_bar, ", GF a_s = ", as_gf
        n_pass = n_pass + 1
    else
        print '(A, F9.5)', " [FAIL] Van der Waals a_bar = ", a_bar
    end if

    ! --------------------------------------------------------------------------
    ! 测试 8: 形状共振 (Shape Resonance) Wigner 时延分析
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    allocate(e_grid(30), delta_scan(30))
    do i = 1, 30
        e_grid(i) = 0.1_dp + real(i - 1, dp) * (0.5_dp - 0.1_dp) / 29.0_dp
        ! 模拟 Breit-Wigner 共振在 E_0 = 0.30 a.u., Gamma = 0.04 a.u.
        ! atan2(Gamma/2, E_0 - E) 随能量单调连续上升 (跨越 0 -> pi/2 -> pi)
        delta_scan(i) = atan2(0.04_dp * 0.5_dp, 0.30_dp - e_grid(i))
    end do
    call analyze_shape_resonance(e_grid, delta_scan, 30, 1.0_dp, res_info, istat)
    if (abs(res_info%e_res - 0.30_dp) < 0.02_dp) then
        print '(A, F8.4, A, F8.4)', " [PASS] Shape resonance detected E_R = ", res_info%e_res, " (Target: 0.30)"
        n_pass = n_pass + 1
    else
        print '(A, F8.4)', " [FAIL] Shape resonance detected E_R = ", res_info%e_res
    end if
    deallocate(e_grid, delta_scan)

    ! --------------------------------------------------------------------------
    ! 测试 9: 双通道密耦定态散射矩阵 (Close-Coupling 2x2 S-Matrix) 幺正性
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    allocate(v11(n_pts), v22(n_pts), v12(n_pts))
    v11 = v_pot
    v22 = v_pot + 0.1_dp
    v12 = 0.05_dp * exp(-(r_grid - 2.0_dp)**2)

    call calc_coupled_channel_smatrix_2x2(r_grid, v11, v22, v12, mass, 0.5_dp, 0.1_dp, &
                                          s_2x2, inel_prob, istat)
    ! 检验通道 1 总通量守恒: |S_11|^2 + |S_12|^2 = 1.0
    s_unitary_norm = abs(s_2x2(1, 1))**2 + abs(s_2x2(1, 2))**2
    if (abs(s_unitary_norm - 1.0_dp) < 0.05_dp .and. inel_prob >= 0.0_dp) then
        print '(A, F9.5, A, F9.5)', " [PASS] 2-Channel S-matrix unitarity |S11|^2+|S12|^2 = ", &
                                     s_unitary_norm, ", P_12 = ", inel_prob
        n_pass = n_pass + 1
    else
        print '(A, F9.5)', " [FAIL] 2-Channel S-matrix unitarity = ", s_unitary_norm
    end if
    deallocate(v11, v22, v12)

    ! --------------------------------------------------------------------------
    ! 测试 10: 全同粒子微分散射截面与量子干涉效应 (Bosons vs Fermions at theta=pi/2)
    ! 在 s-波主导下: Boson 在 pi/2 处截面应为区分粒子的 4 倍; Fermion 在 pi/2 处截面应严格归零
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        real(dp) :: th_pi2(1), ds_dist(1), ds_boson(1), ds_ferm(1)
        th_pi2 = [PI / 2.0_dp]
        call calc_differential_cross_section_identical(0.01_dp, mass, delta_arr, 0, &
                                                      th_pi2, PARTICLE_DISTINGUISHABLE, ds_dist)
        call calc_differential_cross_section_identical(0.01_dp, mass, delta_arr, 0, &
                                                      th_pi2, PARTICLE_IDENTICAL_BOSON, ds_boson)
        call calc_differential_cross_section_identical(0.01_dp, mass, delta_arr, 0, &
                                                      th_pi2, PARTICLE_IDENTICAL_FERMION_POLARIZED, ds_ferm)

        if (abs(ds_boson(1) - 4.0_dp * ds_dist(1)) / ds_boson(1) < 1.0e-5_dp .and. &
            abs(ds_ferm(1)) < 1.0e-12_dp) then
            print '(A, F8.3, A, F8.3, A)', " [PASS] Quantum stats at pi/2: Boson = 4x Dist (", &
                                          ds_boson(1), " vs ", 4.0_dp * ds_dist(1), "), Fermion = 0.0"
            n_pass = n_pass + 1
        else
            print '(A, F8.3, A, F8.3)', " [FAIL] Quantum stats: Boson = ", ds_boson(1), ", Fermion = ", ds_ferm(1)
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 11: 输运截面 (动量传输截面 sigma_m 与 粘滞截面 sigma_v)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        real(dp) :: sig_m, sig_v
        call calc_transport_cross_sections(energy, mass, delta_arr, 3, sig_m, sig_v)
        if (sig_m > 0.0_dp .and. sig_v > 0.0_dp) then
            print '(A, F10.4, A, F10.4)', " [PASS] Transport cross sections: sigma_m = ", sig_m, ", sigma_v = ", sig_v
            n_pass = n_pass + 1
        else
            print '(A, F10.4, A, F10.4)', " [FAIL] Transport cross sections: sigma_m = ", sig_m, ", sigma_v = ", sig_v
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 12: 截面能谱扫描 (Cross Section Spectrum) 与低能极限 4*pi*a_s^2
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        real(dp) :: e_scan(5), sig_tot_scan(5), sig_as_exact
        e_scan = [0.0001_dp, 0.0005_dp, 0.001_dp, 0.005_dp, 0.01_dp]
        call calc_cross_section_spectrum(r_grid, v_pot, mass, e_scan, 5, 1, sig_tot_scan)
        sig_as_exact = 4.0_dp * PI * (as_exact**2)
        ! 最低能量处应当极度接近 4*pi*a_s^2
        if (abs(sig_tot_scan(1) - sig_as_exact) / sig_as_exact < 0.02_dp) then
            print '(A, F9.4, A, F9.4)', " [PASS] Low-energy limit sigma(E->0) = ", sig_tot_scan(1), &
                                         " (4*pi*a_s^2 = ", sig_as_exact, ")"
            n_pass = n_pass + 1
        else
            print '(A, F9.4, A, F9.4)', " [FAIL] Low-energy limit sigma = ", sig_tot_scan(1), &
                                         " (Exact: ", sig_as_exact, ")"
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 13: 微分散射截面 Legendre 展开与各项同性分量 A_0 = sigma_tot / (4*pi)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        real(dp) :: th_grid(181), ds_omega(181), a_coeffs(0:4), fb_asym
        integer  :: j
        do j = 1, 181
            th_grid(j) = real(j - 1, dp) * PI / 180.0_dp
        end do
        call calc_differential_cross_section(energy, mass, delta_arr, 3, th_grid, ds_omega)
        call calc_differential_legendre_expansion(th_grid, ds_omega, 4, a_coeffs, fb_asym)
        if (abs(a_coeffs(0) - sigma_tot / (4.0_dp * PI)) / (sigma_tot / (4.0_dp * PI)) < 0.01_dp) then
            print '(A, F9.4, A, F9.4)', " [PASS] Diff Legendre A_0 = ", a_coeffs(0), &
                                         " (sigma_tot/(4*pi) = ", sigma_tot / (4.0_dp * PI), ")"
            n_pass = n_pass + 1
        else
            print '(A, F9.4)', " [FAIL] Diff Legendre A_0 = ", a_coeffs(0)
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 14: 通用 N 通道定态密耦 (Johnson Log-Derivative) 3 通道全开么正性与非弹性跃迁
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        integer, parameter :: nc = 3, np = 400
        real(dp) :: r_m(np), v_m(nc, nc, np), th(nc), e_tot
        integer  :: l_ch(nc), k_i, k_j, p_idx
        type(multichannel_result_t) :: mc_res
        real(dp) :: r_val, sum_prob, max_unitarity_err

        th = [0.0_dp, 0.05_dp, 0.10_dp]
        l_ch = [0, 0, 0]
        e_tot = 0.25_dp ! 高于所有阈值 -> 3 通道全开 (n_open = 3)

        do p_idx = 1, np
            r_m(p_idx) = 0.5_dp + real(p_idx - 1, dp) * (8.0_dp - 0.5_dp) / real(np - 1, dp)
            r_val = r_m(p_idx)
            ! 对角势能
            v_m(1, 1, p_idx) = -1.0_dp * exp(-(r_val - 2.0_dp)**2)
            v_m(2, 2, p_idx) = -0.8_dp * exp(-(r_val - 2.2_dp)**2)
            v_m(3, 3, p_idx) = -0.6_dp * exp(-(r_val - 2.5_dp)**2)
            ! 通道间非绝热跃迁耦合
            v_m(1, 2, p_idx) = 0.08_dp * exp(-(r_val - 2.1_dp)**2)
            v_m(2, 1, p_idx) = v_m(1, 2, p_idx)
            v_m(2, 3, p_idx) = 0.05_dp * exp(-(r_val - 2.3_dp)**2)
            v_m(3, 2, p_idx) = v_m(2, 3, p_idx)
            v_m(1, 3, p_idx) = 0.02_dp * exp(-(r_val - 2.2_dp)**2)
            v_m(3, 1, p_idx) = v_m(1, 3, p_idx)
        end do

        call calc_multichannel_close_coupling_logder(r_m, v_m, 1.0_dp, e_tot, th, l_ch, mc_res)

        max_unitarity_err = 0.0_dp
        do k_i = 1, mc_res%n_open
            sum_prob = 0.0_dp
            do k_j = 1, mc_res%n_open
                sum_prob = sum_prob + mc_res%prob_matrix(k_j, k_i)
            end do
            max_unitarity_err = max(max_unitarity_err, abs(sum_prob - 1.0_dp))
        end do

        if (mc_res%n_open == 3 .and. max_unitarity_err < 1.0e-3_dp .and. mc_res%prob_matrix(2, 1) > 0.01_dp) then
            print '(A, I1, A, E10.2, A, F8.4)', " [PASS] 3-Channel CC Log-Der: n_open = ", mc_res%n_open, &
                                                  ", Unitarity err = ", max_unitarity_err, &
                                                  ", P(1->2) = ", mc_res%prob_matrix(2, 1)
            n_pass = n_pass + 1
        else
            print '(A, I1, A, E10.2)', " [FAIL] 3-Channel CC Log-Der: n_open = ", mc_res%n_open, &
                                       ", Unitarity err = ", max_unitarity_err
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 15: 开-闭通道共存与 Feshbach 共振散射长度扫描
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        integer, parameter :: nc = 2, np = 400, ne = 10
        real(dp) :: r_m(np), v_m(nc, nc, np), th(nc), e_scan(ne), as_scan(ne), delta_sums(ne)
        integer  :: l_ch(nc), p_idx, ie
        real(dp) :: r_val
        type(multichannel_result_t) :: fb_res

        th = [0.0_dp, 0.40_dp] ! 通道 1 开 (0.05 a.u.), 通道 2 闭 (0.40 a.u.)
        l_ch = [0, 0]

        do p_idx = 1, np
            r_m(p_idx) = 0.6_dp + real(p_idx - 1, dp) * (9.0_dp - 0.6_dp) / real(np - 1, dp)
            r_val = r_m(p_idx)
            v_m(1, 1, p_idx) = -0.5_dp * exp(-(r_val - 2.0_dp)**2)
            v_m(2, 2, p_idx) = -1.2_dp * exp(-(r_val - 2.2_dp)**2) ! 较深闭通道阱，包含 Feshbach 准束缚态
            v_m(1, 2, p_idx) = 0.12_dp * exp(-(r_val - 2.1_dp)**2)
            v_m(2, 1, p_idx) = v_m(1, 2, p_idx)
        end do

        ! 1. 验证单一开通道在闭通道存在下的幺正性 |S_11|^2 = 1.0 (几率不泄漏)
        call calc_multichannel_close_coupling_logder(r_m, v_m, 1.0_dp, 0.08_dp, th, l_ch, fb_res)

        ! 2. 能量扫描
        do ie = 1, ne
            e_scan(ie) = 0.01_dp + real(ie - 1, dp) * (0.25_dp - 0.01_dp) / real(ne - 1, dp)
        end do
        call calc_feshbach_resonance_scan(r_m, v_m, 1.0_dp, e_scan, ne, th, l_ch, as_scan, delta_sums)

        if (fb_res%n_open == 1 .and. fb_res%n_closed == 1 .and. &
            abs(fb_res%prob_matrix(1, 1) - 1.0_dp) < 1.0e-3_dp) then
            print '(A, F9.5, A, F8.4)', " [PASS] Feshbach Open-Closed CC: |S11|^2 = ", fb_res%prob_matrix(1, 1), &
                                         ", a_s(E=0.08) = ", as_scan(4)
            n_pass = n_pass + 1
        else
            print '(A, F9.5)', " [FAIL] Feshbach Open-Closed CC: |S11|^2 = ", fb_res%prob_matrix(1, 1)
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 16: 非含时散射能量本征波函数 u_{l, E}(r) 求解与正交连续态归一化
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        integer, parameter :: nw = 1200
        real(dp) :: r_w(nw), v_w(nw), u_wf_e(nw), u_wf_k(nw), u_wf_1(nw)
        real(dp) :: delta_w, e_test, k_test, amp_num, target_amp_e, target_amp_k
        integer  :: iw

        do iw = 1, nw
            r_w(iw) = 0.01_dp + real(iw - 1, dp) * (25.0_dp - 0.01_dp) / real(nw - 1, dp)
            ! 高斯势阱
            v_w(iw) = -0.8_dp * exp(-(r_w(iw) - 2.0_dp)**2)
        end do
        e_test = 0.20_dp
        k_test = sqrt(2.0_dp * 1.0_dp * e_test)
        target_amp_e = sqrt(2.0_dp * 1.0_dp / (PI * k_test))
        target_amp_k = sqrt(2.0_dp / PI)

        call calc_scattering_wavefunction_ti(r_w, v_w, 1.0_dp, e_test, 0, NORM_ENERGY, u_wf_e, delta_w)
        call calc_scattering_wavefunction_ti(r_w, v_w, 1.0_dp, e_test, 0, NORM_MOMENTUM, u_wf_k, delta_w)
        call calc_scattering_wavefunction_ti(r_w, v_w, 1.0_dp, e_test, 0, NORM_UNIT_AMPLITUDE, u_wf_1, delta_w)

        ! 检验渐近区振幅
        amp_num = sqrt(u_wf_e(nw)**2 + (((u_wf_e(nw) - u_wf_e(nw - 1)) / (r_w(2) - r_w(1))) / k_test)**2)

        if (abs(amp_num - target_amp_e) / target_amp_e < 0.02_dp .and. &
            abs(u_wf_e(1)) < 1.0e-4_dp .and. &
            abs(maxval(abs(u_wf_1(nw-50:nw))) - 1.0_dp) < 0.05_dp) then
            print '(A, F8.4, A, F8.4, A, F8.4)', " [PASS] Continuous eigenstate u_E(r): Amp = ", amp_num, &
                                                 " (Exact: ", target_amp_e, "), delta = ", delta_w
            n_pass = n_pass + 1
        else
            print '(A, F8.4, A, F8.4)', " [FAIL] Continuous eigenstate u_E(r): Amp = ", amp_num, &
                                        " (Exact: ", target_amp_e, ")"
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 17: 多扇区分段网格 (Segmented Grid) 零能散射长度求解与解析解对比
    ! 采用 3 扇区网格: [0.01, 2.5] 细网格, [2.5, 6.0] 中网格, [6.0, 15.0] 粗网格
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        type(segmented_grid_t) :: s_grid
        real(dp) :: bounds(3), steps(3)
        real(dp), allocatable :: v_seg(:)
        real(dp) :: as_seg
        integer  :: isg

        bounds = [2.5_dp, 6.0_dp, 15.0_dp]
        steps  = [0.005_dp, 0.02_dp, 0.05_dp]
        call create_segmented_grid(0.01_dp, bounds, steps, s_grid)

        allocate(v_seg(s_grid%n_total))
        do isg = 1, s_grid%n_total
            if (s_grid%r(isg) <= r_well) then
                v_seg(isg) = -v0_well
            else
                v_seg(isg) = 0.0_dp
            end if
        end do

        call calc_scattering_length_segmented_numerov(s_grid, v_seg, mass, as_seg)
        if (abs(as_seg - as_exact) / abs(as_exact) < 0.005_dp .and. s_grid%n_sectors == 3) then
            print '(A, I2, A, I5, A, F9.5, A, F9.5)', " [PASS] Segmented Numerov (", s_grid%n_sectors, &
                  " sectors, N=", s_grid%n_total, ") a_s = ", as_seg, " (Exact: ", as_exact, ")"
            n_pass = n_pass + 1
        else
            print '(A, F9.5, A, F9.5)', " [FAIL] Segmented Numerov a_s = ", as_seg, " (Exact: ", as_exact, ")"
        end if
        deallocate(v_seg)
    end block

    ! --------------------------------------------------------------------------
    ! 测试 18: 分段扇区网格波函数 u_{l, E}(r) 跨扇区平滑推进与相移验证
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        type(segmented_grid_t) :: s_grid
        real(dp) :: bounds(3), steps(3)
        real(dp), allocatable :: v_seg(:), u_seg(:)
        real(dp) :: delta_seg, k_m, cs_m
        complex(dp) :: s_m, t_m
        real(dp) :: e_t, k_t, amp_t, amp_seg
        integer  :: isg, n_pts_seg

        bounds = [3.0_dp, 8.0_dp, 25.0_dp]
        steps  = [0.005_dp, 0.02_dp, 0.10_dp]
        call create_segmented_grid(0.01_dp, bounds, steps, s_grid)

        n_pts_seg = s_grid%n_total
        allocate(v_seg(n_pts_seg), u_seg(n_pts_seg))

        do isg = 1, n_pts_seg
            v_seg(isg) = -0.8_dp * exp(-(s_grid%r(isg) - 2.0_dp)**2)
        end do

        e_t = 0.20_dp
        k_t = sqrt(2.0_dp * 1.0_dp * e_t)
        amp_t = sqrt(2.0_dp * 1.0_dp / (PI * k_t))

        call calc_scattering_wavefunction_segmented_ti( &
            s_grid, v_seg, 1.0_dp, e_t, 0, NORM_ENERGY, u_seg, delta_seg)

        call calc_phase_shift_segmented( &
            s_grid, v_seg, 1.0_dp, e_t, 0, delta_seg, k_m, s_m, t_m, cs_m)

        amp_seg = sqrt(u_seg(n_pts_seg)**2 + (((u_seg(n_pts_seg) - u_seg(n_pts_seg - 1)) / &
                       s_grid%sector_dr(s_grid%n_sectors)) / k_t)**2)

        if (abs(amp_seg - amp_t) / amp_t < 0.02_dp .and. &
            abs(abs(s_m) - 1.0_dp) < 1.0e-5_dp) then
            print '(A, F8.4, A, F8.4, A, F8.4)', " [PASS] Segmented Wavefunction: Amp = ", amp_seg, &
                  " (Exact: ", amp_t, "), |S| = ", abs(s_m)
            n_pass = n_pass + 1
        else
            print '(A, F8.4, A, F8.4)', " [FAIL] Segmented Wavefunction: Amp = ", amp_seg, &
                  " (Exact: ", amp_t, ")"
        end if
        deallocate(v_seg, u_seg)
    end block

    ! --------------------------------------------------------------------------
    ! 测试 19: 分段扇区网格多通道密耦 Log-Derivative 求解与 S-矩阵幺正性
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        type(segmented_grid_t) :: s_grid
        real(dp) :: bounds(2), steps(2)
        real(dp), allocatable :: v_mat_seg(:, :, :)
        type(multichannel_result_t) :: mc_res
        real(dp) :: thresh(3), e_coll
        integer  :: l_ch(3), isg
        real(dp) :: p_sum, u_err

        bounds = [5.0_dp, 20.0_dp]
        steps  = [0.01_dp, 0.05_dp]
        call create_segmented_grid(0.1_dp, bounds, steps, s_grid)

        allocate(v_mat_seg(3, 3, s_grid%n_total))
        v_mat_seg = 0.0_dp
        do isg = 1, s_grid%n_total
            v_mat_seg(1, 1, isg) = -0.5_dp * exp(-(s_grid%r(isg) - 2.0_dp)**2)
            v_mat_seg(2, 2, isg) = -0.3_dp * exp(-(s_grid%r(isg) - 2.5_dp)**2)
            v_mat_seg(3, 3, isg) = -0.2_dp * exp(-(s_grid%r(isg) - 3.0_dp)**2)
            v_mat_seg(1, 2, isg) =  0.08_dp * exp(-(s_grid%r(isg) - 2.2_dp)**2)
            v_mat_seg(2, 1, isg) =  v_mat_seg(1, 2, isg)
        end do

        thresh = [0.0_dp, 0.05_dp, 0.10_dp]
        l_ch   = [0, 0, 0]
        e_coll = 0.25_dp

        call calc_multichannel_close_coupling_segmented_logder( &
            s_grid, v_mat_seg, 1.0_dp, e_coll, thresh, l_ch, mc_res)

        p_sum = mc_res%prob_matrix(1, 1) + mc_res%prob_matrix(2, 1) + mc_res%prob_matrix(3, 1)
        u_err = abs(p_sum - 1.0_dp)

        if (mc_res%n_open == 3 .and. u_err < 1.0e-5_dp) then
            print '(A, I2, A, ES10.2, A, F8.4)', " [PASS] Segmented CC Log-Der: n_open = ", mc_res%n_open, &
                  ", Unitarity err = ", u_err, ", P(1->2) = ", mc_res%prob_matrix(2, 1)
            n_pass = n_pass + 1
        else
            print '(A, ES10.2)', " [FAIL] Segmented CC Log-Der Unitarity err = ", u_err
        end if
        deallocate(v_mat_seg)
    end block

    deallocate(r_grid, v_pot)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "TI Scattering Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', "SUCCESS: All time-independent scattering tests passed."
    else
        stop 1
    end if
end program test_ti_scattering
