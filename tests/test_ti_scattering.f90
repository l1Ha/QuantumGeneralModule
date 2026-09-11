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

    deallocate(r_grid, v_pot)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "TI Scattering Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', "SUCCESS: All time-independent scattering tests passed."
    else
        stop 1
    end if
end program test_ti_scattering
