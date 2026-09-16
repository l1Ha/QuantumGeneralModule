! ==============================================================================
! GeneralModule 物理应用算例 ex36
! 主题: 亚稳态稀有气体碰撞潘宁电离与缔合电离动力学全流程仿真
!       (Penning Ionization, Associative Ionization & PIES Spectroscopy)
! 体系: He*(2^3S) + Ar -> He + Ar+ + e- (PI) 及 HeAr+ + e- (AI)
! 物理:
!   1. 复光学势相互作用与自电离跃迁宽度 Gamma(R) 空间分布;
!   2. 碰撞能量依赖的总电离截面、缔合电离分流比与能量阈值效应;
!   3. 潘宁电子能谱 (PIES) 局域核构型发射谱线重构;
!   4. 麦克斯韦热平衡速率常数与低温至室温反应动力学;
!   5. 自旋守恒定则与极化超冷亚稳态原子寿命抑制机理。
! ==============================================================================

program ex36_penning_associative_ionization
    use iso_fortran_env, only: dp => real64
    use general_module
    implicit none

    type(penning_system_t) :: sys, sys_polarized
    type(penning_cross_section_result_t) :: cs_res
    integer :: stat
    integer :: ie, it, ip
    real(dp) :: e_coll, t_kelvin
    real(dp) :: k_tot, k_pi, k_ai
    real(dp) :: e_grid(60), pies(60)
    real(dp) :: a_unpol, beta_unpol, k_unpol
    real(dp) :: a_pol, beta_pol, k_pol
    real(dp) :: v_star, v_plus, gamma_ev

    print *, "======================================================================"
    print *, "  Example 36: Penning & Associative Ionization Dynamics Simulation   "
    print *, "======================================================================"
    print *, "  Reaction Channel 1 (PI): He*(2^3S) + Ar -> He + Ar+ + e-"
    print *, "  Reaction Channel 2 (AI): He*(2^3S) + Ar -> HeAr+ + e-"
    print *, "----------------------------------------------------------------------"

    ! --------------------------------------------------------------------------
    ! 1. 初始化基准体系与相互作用势特性
    ! --------------------------------------------------------------------------
    call init_penning_preset_he_star_ar(sys, stat=stat)
    if (stat /= 0) stop "Failed to initialize Penning system."

    print *, "[1] System Configuration & Interaction Potentials:"
    print '(A, A32)',          "    System Name:       ", sys%name
    print '(A, F8.3, A)',      "    Excitation Energy: ", sys%e_exc_ev, " eV"
    print '(A, F8.3, A)',      "    Target IP (Ar):    ", sys%ip_target_ev, " eV"
    print '(A, F8.3, A)',      "    Excess Energy:     ", sys%excess_energy_ev, " eV"
    print '(A, F8.4, A, F6.2, A)', "    Neutral V*(R):     Well = ", sys%v_star_de_ev, &
                               " eV, Re = ", sys%v_star_re_ang, " A"
    print '(A, F8.4, A, F6.2, A)', "    Ionic V+(R) HeAr+: Well = ", sys%v_plus_de_ev, &
                               " eV, Re = ", sys%v_plus_re_ang, " A"
    print *

    ! 输出势能与自电离跃迁宽度采样
    print *, "    --- Radial Potential and Autoionization Width Profile ---"
    print *, "    R (A)     V*(R) (meV)    V+(R) (eV)    Gamma(R) (meV)"
    print *, "   ---------------------------------------------------------"
    do ip = 1, 5
        e_coll = 2.5_dp + real(ip - 1, dp) * 0.75_dp
        call eval_penning_potentials(sys, e_coll, v_star, v_plus, gamma_ev)
        print '(4X, F5.2, 5X, F10.3, 5X, F9.4, 5X, F10.3)', &
              e_coll, v_star * 1000.0_dp, v_plus, gamma_ev * 1000.0_dp
    end do
    print *

    ! --------------------------------------------------------------------------
    ! 2. 碰撞能量扫描与缔合电离分支比 F_AI(E)
    ! --------------------------------------------------------------------------
    print *, "[2] Collision Energy Scan & Ionization Branching Ratio:"
    print *, "    E_coll (eV)   sigma_tot (A^2)   sigma_PI (A^2)   sigma_AI (A^2)   F_AI (%)"
    print *, "   ------------------------------------------------------------------------"
    do ie = 1, 6
        ! 能量从超冷/低热能 (0.01 eV) 到超热能 (0.50 eV)
        e_coll = 0.01_dp + real(ie - 1, dp)**2 * 0.02_dp
        call calc_penning_cross_sections(sys, e_coll, n_b=128, b_max_ang=8.0_dp, res=cs_res)

        print '(4X, F7.3, 6X, F10.2, 7X, F10.2, 7X, F10.2, 6X, F7.2, A)', &
              cs_res%e_coll_ev, cs_res%sigma_tot_ang2, cs_res%sigma_pi_ang2, &
              cs_res%sigma_ai_ang2, cs_res%ai_branching_ratio * 100.0_dp, " %"
    end do
    print *, "    >> Observation: Associative ionization (AI) dominates at low"
    print *, "       collision energies where particles are trapped in HeAr+ well."
    print *

    ! --------------------------------------------------------------------------
    ! 3. 潘宁电离电子能谱 (PIES)
    ! --------------------------------------------------------------------------
    print *, "[3] Penning Ionization Electron Spectrum (PIES) at E_coll = 0.05 eV:"
    call calc_penning_electron_spectrum(sys, e_coll_ev=0.05_dp, n_e=60, &
                                        e_min_ev=3.70_dp, e_max_ev=4.50_dp, &
                                        e_grid_ev=e_grid, pies_intensity=pies)
    print *, "    E_e (eV)      Intensity (Normalized a.u.)"
    print *, "   --------------------------------------------"
    do ie = 1, 60, 10
        print '(4X, F7.3, 8X, F12.5)', e_grid(ie), pies(ie)
    end do
    print *, "    >> PIES Peak accurately reflects (V* - V+ + E_excess) transition region."
    print *

    ! --------------------------------------------------------------------------
    ! 4. 宽温区麦克斯韦反应速率常数 k(T)
    ! --------------------------------------------------------------------------
    print *, "[4] Thermal Maxwell-Boltzmann Rate Coefficients k(T):"
    print *, "    Temp (K)      k_tot (cm^3/s)        k_PI (cm^3/s)         k_AI (cm^3/s)"
    print *, "   ------------------------------------------------------------------------"
    do it = 1, 5
        t_kelvin = 100.0_dp + real(it - 1, dp) * 150.0_dp
        call calc_penning_thermal_rate(sys, temp_k=t_kelvin, n_e=64, &
                                      k_tot_cm3_s=k_tot, k_pi_cm3_s=k_pi, &
                                      k_ai_cm3_s=k_ai)
        print '(4X, F6.1, 7X, ES14.5, 7X, ES14.5, 7X, ES14.5)', &
              t_kelvin, k_tot, k_pi, k_ai
    end do
    print *

    ! --------------------------------------------------------------------------
    ! 5. 超冷自旋抑制效应 (Spin-Polarization Suppression)
    ! --------------------------------------------------------------------------
    print *, "[5] Ultracold Regime & Spin-Polarized Suppression (BEC Preservation):"
    call calc_ultracold_penning_complex_length(sys, a_unpol, beta_unpol, k_unpol)

    ! 自旋极化态: S=2 (Quintet 态)，因自旋角动量守恒禁阻电离，电离率降低 10^4 倍
    sys_polarized = sys
    sys_polarized%spin_factor = 1.0e-4_dp
    call calc_ultracold_penning_complex_length(sys_polarized, a_pol, beta_pol, k_pol)

    print '(A, ES12.4, A)', "    Unpolarized loss rate K_loss: ", k_unpol, " cm^3/s"
    print '(A, ES12.4, A)', "    Spin-polarized rate K_loss:   ", k_pol, " cm^3/s"
    print '(A, F10.1, A)',  "    Suppression Factor:           ", k_unpol / max(1.0e-30_dp, k_pol), "x"
    print *, "    >> Spin conservation guarantees stability of metastable helium BEC!"
    print *

    print *, "======================================================================"
    print *, "  Example 36 Completed Successfully!                                  "
    print *, "======================================================================"

end program ex36_penning_associative_ionization
