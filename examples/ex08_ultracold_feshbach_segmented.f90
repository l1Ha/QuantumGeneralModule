! ==============================================================================
! GeneralModule Example 08:
! 超冷双原子磁 Feshbach 共振与多扇区分段网格密耦动力学 (Segmented Multi-Sector CC)
! ==============================================================================
! 本算例演示超冷物理中极具代表性的两通道/多通道磁 Feshbach 共振散射色散特征：
!   1. [分段扇区网格 (Segmented Grid)]:
!      短程深势阱区 (0.05 - 3.0 a0) 采用密集微步长 (dr = 0.005 a0);
!      中程范德华色散过渡区 (3.0 - 15.0 a0) 放大步长 (dr = 0.02 a0);
!      长程渐近自由演化区 (15.0 - 100.0 a0) 采用大步长 (dr = 0.10 a0);
!      总格点数压缩 90% 以上，极大地加速多通道矩阵逆推。
!   2. [磁场诱导开-闭通道耦合]:
!      开通道 (入口超冷碰撞态) 与闭通道 (高能塞曼超精细分子准束缚态)
!      在特定磁场 B_0 附近由于自旋交换产生简并共振。
!   3. [Feshbach 色散线型提取与非线性拟合]:
!      计算不同磁场下的散射长度 a(B)，提取 Feshbach 共振特征：
!         a(B) = a_bg * [ 1 - Delta / (B - B0) ]
!      并将完整色散扫描曲线写入 feshbach_segmented_scan.dat。
! ==============================================================================
program ex08_ultracold_feshbach_segmented
    use general_module
    implicit none

    type(segmented_grid_t) :: s_grid
    type(cold_atom_t) :: rb87
    type(multichannel_result_t) :: mc_res
    real(dp) :: b_res_fit, delta_b_fit, a_bg_fit

    real(dp) :: bounds(3), steps(3)
    real(dp), allocatable :: v_mat(:, :, :)
    real(dp) :: mass_rb, e_coll, k1
    real(dp) :: thresholds(2)
    integer  :: l_channels(2)
    integer  :: n_b_pts, ib, ig, file_unit
    real(dp) :: b_min, b_max, db, b_curr, delta_mu
    real(dp), allocatable :: b_fields(:), a_scat(:), p12_prob(:)

    print '(A)', "===================================================================="
    print '(A)', "   GeneralModule Example 08: Ultracold Feshbach Resonances          "
    print '(A)', "         [Segmented Multi-Sector Grid Close-Coupling]               "
    print '(A)', "===================================================================="

    ! --------------------------------------------------------------------------
    ! 1. 设置 87Rb 原子参数与相对折合质量
    ! --------------------------------------------------------------------------
    call get_cold_atom_preset("87Rb", rb87)
    ! 87Rb+87Rb 折合质量 (a.u.): mu = M_Rb / 2
    mass_rb = (rb87%mass_amu * 1822.888486_dp) / 2.0_dp

    print '(A, F9.4, A, F10.2, A)', " [Setup] 87Rb + 87Rb Reduced Mass = ", &
                                    rb87%mass_amu / 2.0_dp, " amu (", mass_rb, " a.u.)"

    ! --------------------------------------------------------------------------
    ! 2. 构造三扇区分段径向网格 (Sector Segmented Radial Grid)
    ! --------------------------------------------------------------------------
    bounds = [3.0_dp, 15.0_dp, 100.0_dp]
    steps  = [0.005_dp, 0.02_dp, 0.10_dp]
    call create_segmented_grid(0.05_dp, bounds, steps, s_grid)

    print '(A, I2, A, I5, A)', " [Grid] Segmented Grid: ", s_grid%n_sectors, &
                               " Sectors, Total Points = ", s_grid%n_total, " pts."
    print '(A, F7.4, A, I5, A)', "        Sector 1: [0.05, 3.0] a0, dr = ", s_grid%sector_dr(1), &
          " a0 (", s_grid%sector_npts(1), " pts)"
    print '(A, F7.4, A, I5, A)', "        Sector 2: [3.0, 15.0] a0, dr = ", s_grid%sector_dr(2), &
          " a0 (", s_grid%sector_npts(2), " pts)"
    print '(A, F7.4, A, I5, A)', "        Sector 3: [15.0, 100.0] a0, dr = ", s_grid%sector_dr(3), &
          " a0 (", s_grid%sector_npts(3), " pts)"

    ! --------------------------------------------------------------------------
    ! 3. 构造双通道模型势能矩阵
    !    通道 1: 开通道 (入口碰撞态，渐近阈值 E_1 = 0)
    !    通道 2: 闭通道 (束缚准稳态，渐近阈值随磁场线性移动 E_2(B) = E_bound0 + delta_mu * B)
    ! --------------------------------------------------------------------------
    allocate(v_mat(2, 2, s_grid%n_total))
    v_mat = 0.0_dp

    do ig = 1, s_grid%n_total
        ! 开通道势阱 (Morse/Gaussian 模拟短程分子阱)
        v_mat(1, 1, ig) = -0.60_dp * exp(-(s_grid%r(ig) - 2.5_dp)**2 / 0.8_dp)
        ! 闭通道分子态势阱 (更深，极小值位置略微内移)
        v_mat(2, 2, ig) = -0.90_dp * exp(-(s_grid%r(ig) - 2.2_dp)**2 / 0.8_dp)
        ! 通道间自旋交换非绝热耦合项 V_12(r)
        v_mat(1, 2, ig) =  0.05_dp * exp(-(s_grid%r(ig) - 2.4_dp)**2 / 0.5_dp)
        v_mat(2, 1, ig) =  v_mat(1, 2, ig)
    end do

    ! 碰撞质心动能 (超冷温度范围: E_coll ~ 1 uK - 100 nK)
    e_coll = 1.0e-7_dp
    l_channels = [0, 0] ! s-波碰撞

    ! 磁矩差: delta_mu ~ 1.5 mu_B (以 a.u. 为单位: mu_B_au)
    delta_mu = 1.5_dp * MU_B_AU

    ! --------------------------------------------------------------------------
    ! 4. 磁场扫描求解分段多通道密耦
    ! --------------------------------------------------------------------------
    b_min = 70.0_dp
    b_max = 90.0_dp
    n_b_pts = 41
    db = (b_max - b_min) / real(n_b_pts - 1, dp)

    allocate(b_fields(n_b_pts), a_scat(n_b_pts), p12_prob(n_b_pts))

    print '(A)', " [Scan] Scanning Magnetic Field B in [70, 90] Gauss..."
    print '(A)', "   B (Gauss)    Threshold E2 (a.u.)    Scattering Length a(B) (a0)   Unitarity Error"
    print '(A)', " -----------------------------------------------------------------------------------"

    do ib = 1, n_b_pts
        b_curr = b_min + real(ib - 1, dp) * db
        b_fields(ib) = b_curr

        ! 闭通道能量阈值 (以 B0 ~ 80 Gauss 为共振交叉点)
        thresholds(1) = 0.0_dp
        thresholds(2) = (b_curr - 80.0_dp) * delta_mu

        call calc_multichannel_close_coupling_segmented_logder( &
            s_grid, v_mat, mass_rb, e_coll, thresholds, l_channels, mc_res)

        if (mc_res%n_open >= 1) then
            k1 = mc_res%k_open(1)
            a_scat(ib) = -mc_res%k_matrix(1, 1) / max(1.0e-12_dp, k1)
        else
            a_scat(ib) = 0.0_dp
        end if

        if (ib == 1 .or. ib == n_b_pts .or. mod(ib, 5) == 0) then
            print '(F10.2, ES20.4, F25.4, ES22.2)', b_curr, thresholds(2), a_scat(ib), &
                  abs(sum(mc_res%prob_matrix(:, 1)) - 1.0_dp)
        end if
    end do
    print '(A)', " -----------------------------------------------------------------------------------"

    ! --------------------------------------------------------------------------
    ! 5. 非线性最小二乘拟合提取 Feshbach 共振参数 (B0, Delta, a_bg)
    ! --------------------------------------------------------------------------
    call fit_feshbach_resonance_parameters(b_fields, a_scat, n_b_pts, b_res_fit, delta_b_fit, a_bg_fit)

    print '(A)', " [Fit] Feshbach Resonance Parameters Extraction:"
    print '(A, F9.3, A)', "       Resonance Position B_0 = ", b_res_fit, " Gauss"
    print '(A, F9.3, A)', "       Resonance Width Delta  = ", delta_b_fit, " Gauss"
    print '(A, F9.3, A)', "       Background Length a_bg = ", a_bg_fit, " a0"
    print '(A, F9.3, A)', "       Zero-Crossing B_zero   = ", b_res_fit + delta_b_fit, " Gauss"

    ! --------------------------------------------------------------------------
    ! 6. 保存扫描数据至文件
    ! --------------------------------------------------------------------------
    open(newunit=file_unit, file="feshbach_segmented_scan.dat", status="replace", action="write")
    write(file_unit, '(A)') "# =================================================================="
    write(file_unit, '(A)') "# GeneralModule Ex08: Magnetic Feshbach Resonance Scan"
    write(file_unit, '(A)') "# System: 87Rb + 87Rb (Two-Channel Segmented Log-Derivative CC)"
    write(file_unit, '(A, F8.3, A, F8.3, A, F8.3)') "# Fitted: B0 = ", b_res_fit, &
          " G, Delta = ", delta_b_fit, " G, a_bg = ", a_bg_fit
    write(file_unit, '(A)') "# Column 1: Magnetic Field B (Gauss)"
    write(file_unit, '(A)') "# Column 2: Scattering Length a_s(B) (Bohr / a0)"
    write(file_unit, '(A)') "# Column 3: Closed-Channel Energy Threshold E_2(B) (a.u.)"
    write(file_unit, '(A)') "# =================================================================="
    do ib = 1, n_b_pts
        write(file_unit, '(F12.4, 2X, ES18.8, 2X, ES18.8)') b_fields(ib), a_scat(ib), &
              (b_fields(ib) - 80.0_dp) * delta_mu
    end do
    close(file_unit)

    print '(A)', ">> Feshbach resonance scan data saved to: feshbach_segmented_scan.dat"
    print '(A)', "SUCCESS: Example 08 completed successfully."

    deallocate(v_mat, b_fields, a_scat, p12_prob)
end program ex08_ultracold_feshbach_segmented
