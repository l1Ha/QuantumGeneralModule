! ==============================================================================
! GeneralModule Unit Tests: Anisotropic Dipolar Quantum Scattering
! Testing:
!   1. C_{2, q} spatial quadrupole matrix elements & parity/projection rules
!   2. [s1 x s2]^(2) spin tensor matrix elements & S=1 triplet selection
!   3. Total MDDI matrix elements with strict M_tot = M_L + M_S conservation
!   4. Inelastic dipolar relaxation cross section & thermal rate K_rel(T)
!   5. Polar molecule Stark induced dipole d_ind(E) from weak to strong field
!   6. Universal dipole length scale a_d and multi-partial-wave potential matrix
! ==============================================================================
program test_dipolar_scattering
    use mod_constants, only: dp, PI
    use mod_dipolar_scattering
    implicit none

    integer :: n_pass = 0, n_total = 0
    real(dp) :: tol = 1.0e-5_dp

    ! 临时测试变量
    real(dp) :: val, exact, d_ind, a_d
    type(dipolar_relaxation_result_t) :: rel_res
    type(polar_molecule_t) :: krb_mol
    type(dipolar_channel_t) :: ch(3)
    real(dp) :: v_mat(3, 3), k_therm
    integer  :: st

    print *, "=================================================="
    print *, "   GeneralModule Unit Tests: Dipolar Scattering   "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! 1. 空间球张量 C_{2, q} 矩阵元检验
    ! --------------------------------------------------------------------------
    ! s-波之间耦合严格为 0: <0,0|C_20|0,0> = 0
    val = calc_mddi_spatial_matrix_element(0, 0, 0, 0, 0)
    call assert_real_equal(val, 0.0_dp, "Spatial C_20 <0,0|0,0> = 0", n_pass, n_total)

    ! 宇称相消: <1,0|C_20|0,0> = 0 (L+L'=1 为奇数)
    val = calc_mddi_spatial_matrix_element(1, 0, 0, 0, 0)
    call assert_real_equal(val, 0.0_dp, "Spatial parity forbidden <1,0|C_20|0,0> = 0", n_pass, n_total)

    ! s-d 波耦合: <2,0|C_20|0,0> = (-1)^0 * sqrt(5) * (2 2 0; 0 0 0) * (2 2 0; 0 0 0) = 1/sqrt(5)
    val = calc_mddi_spatial_matrix_element(2, 0, 0, 0, 0)
    exact = 1.0_dp / sqrt(5.0_dp)
    call assert_real_equal(val, exact, "Spatial C_20 <2,0|0,0> = 1/sqrt(5)", n_pass, n_total)

    ! q = 1 分量: <2,1|C_21|0,0>
    ! (2 2 0; -1 1 0) = -1/sqrt(5), 相位 (-1)^1 = -1 => (-1)*sqrt(5)*(1/sqrt(5))*(-1/sqrt(5)) = 1/sqrt(5)
    val = calc_mddi_spatial_matrix_element(2, 1, 0, 0, 1)
    exact = 1.0_dp / sqrt(5.0_dp)
    call assert_real_equal(val, exact, "Spatial C_21 <2,1|0,0> = 1/sqrt(5)", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 2. 自旋 2 阶张量 [s1 x s2]^(2) 矩阵元检验
    ! --------------------------------------------------------------------------
    ! 单重态 S=0 严格为 0
    val = calc_mddi_spin_matrix_element(0, 0, 0, 0, 0)
    call assert_real_equal(val, 0.0_dp, "Spin tensor S=0 singlet = 0", n_pass, n_total)

    ! 三重态 S=1 矩阵元: (1 2 1; 0 0 0) = +sqrt(2/15), (-1)^{1-0} = -1
    ! red_mat = 1.5 => <1,0|T^(2)_0|1,0> = -1 * sqrt(2/15) * 1.5 = -1.5*sqrt(2/15)
    val = calc_mddi_spin_matrix_element(1, 0, 1, 0, 0)
    exact = -1.5_dp * sqrt(2.0_dp / 15.0_dp)
    call assert_real_equal(val, exact, "Spin tensor <1,0|T^(2)_0|1,0>", n_pass, n_total)

    ! q 跃迁分量 Delta_M_S = -1: <1,0|T^(2)_{-1}|1,1>
    val = calc_mddi_spin_matrix_element(1, 0, 1, 1, -1)
    call assert_true(abs(val) > 0.1_dp, "Spin tensor Delta_M_S = -1 transition non-zero", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 3. 总各向异性磁偶极相互作用 (MDDI) 矩阵元检验
    ! --------------------------------------------------------------------------
    ! M_tot 守恒检验: M_L + M_S = M_L' + M_S'
    ! 初态: |S=1, Ms=1, L=0, Ml=0> (M_tot = 1)
    ! 末态 A: |S=1, Ms=0, L=2, Ml=1> (M_tot = 1) => 守恒，非零
    val = calc_mddi_total_matrix_element(1, 0, 2, 1, 1, 1, 0, 0)
    call assert_true(abs(val) > 0.05_dp, "Total MDDI s-d coupling (M_tot conserved)", n_pass, n_total)

    ! 末态 B: |S=1, Ms=0, L=2, Ml=0> (M_tot = 0 /= 1) => 违背总投影守恒，必须严格为 0
    val = calc_mddi_total_matrix_element(1, 0, 2, 0, 1, 1, 0, 0)
    call assert_real_equal(val, 0.0_dp, "Total MDDI projection violation = 0", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 4. 非弹性自旋弛豫 (Dipolar Relaxation) 计算检验
    ! --------------------------------------------------------------------------
    ! 模拟 87Rb 原子在 B = 10 Gauss, 碰撞动能 E = 1.0e-10 a.u. (~30 uK)
    call calc_dipolar_relaxation_cross_section(86.909_dp, 10.0_dp, 1.0e-10_dp, 10.0_dp, rel_res, st)
    call assert_true(st == 0, "Dipolar relaxation calculation status = 0", n_pass, n_total)
    call assert_true(rel_res%e_released > 0.0_dp, "Zeeman energy released > 0", n_pass, n_total)
    call assert_true(rel_res%k_exit > rel_res%k_incident, "Exit momentum k_f > k_i", n_pass, n_total)
    call assert_true(rel_res%rate_coeff_cm3_s > 0.0_dp, "Relaxation rate K_rel > 0", n_pass, n_total)

    ! 热平均速率系数检验
    call calc_dipolar_relaxation_thermal_rate(86.909_dp, 10.0_dp, 1.0e-4_dp, 10.0_dp, k_therm, st)
    call assert_true(st == 0 .and. k_therm > 0.0_dp, "Thermal dipolar relaxation rate > 0", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 5. 极性分子 Stark 效应与实验室系诱导电偶极矩 d_ind
    ! --------------------------------------------------------------------------
    krb_mol%name         = "KRb"
    krb_mol%mass_amu     = 126.37_dp
    krb_mol%b_rot_cm1    = 0.037_dp
    krb_mol%dipole_debye = 0.566_dp

    ! 零电场下诱导偶极矩为 0
    call calc_stark_induced_dipole(krb_mol, 0.0_dp, d_ind, st)
    call assert_real_equal(d_ind, 0.0_dp, "Zero-field induced dipole = 0", n_pass, n_total)

    ! 弱场 E = 2.0 kV/cm 下诱导偶极矩单调增加且小于永久偶极矩
    call calc_stark_induced_dipole(krb_mol, 2.0_dp, d_ind, st)
    call assert_true(d_ind > 0.01_dp .and. d_ind < krb_mol%dipole_debye, &
                     "Weak-field induced dipole 0 < d_ind < d_0", n_pass, n_total)

    ! 极高电场 E = 100.0 kV/cm 下诱导偶极矩趋向饱和 (> 70% d_0)
    call calc_stark_induced_dipole(krb_mol, 100.0_dp, d_ind, st)
    call assert_true(d_ind > 0.7_dp * krb_mol%dipole_debye, &
                     "Strong-field induced dipole saturation (> 0.7 d_0)", n_pass, n_total)

    ! 偶极特征长度尺度 a_d > 0
    a_d = calc_dipole_length_scale(krb_mol%mass_amu, d_ind)
    call assert_true(a_d > 10.0_dp, "Dipole length scale a_d > 10 a_0", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 6. 多分波电偶极耦合势能矩阵组装
    ! --------------------------------------------------------------------------
    ch(1)%l = 0; ch(1)%ml = 0; ch(1)%energy_thresh = 0.0_dp
    ch(2)%l = 2; ch(2)%ml = 0; ch(2)%energy_thresh = 0.0_dp
    ch(3)%l = 4; ch(3)%ml = 0; ch(3)%energy_thresh = 0.0_dp

    call build_dipolar_coupled_potential_matrix(ch, 3, krb_mol%mass_amu, 5000.0_dp, &
                                                d_ind, 20.0_dp, v_mat, st)
    call assert_true(st == 0, "Dipolar potential matrix assembly status = 0", n_pass, n_total)
    call assert_true(abs(v_mat(1, 2) - v_mat(2, 1)) < 1.0e-14_dp, &
                     "Dipolar potential matrix symmetry V_12 == V_21", n_pass, n_total)
    call assert_true(abs(v_mat(1, 2)) > 1.0e-12_dp, "Non-zero s-d dipole coupling V_12", n_pass, n_total)

    ! --------------------------------------------------------------------------
    ! 汇总报告
    ! --------------------------------------------------------------------------
    print *, "--------------------------------------------------"
    print '(A, I3, A, I3, A)', "Dipolar Scattering Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All dipolar scattering tests passed."
    else
        stop 1
    end if

contains

    subroutine assert_real_equal(actual, expected, desc, p_count, t_count)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: desc
        integer, intent(inout) :: p_count, t_count
        t_count = t_count + 1
        if (abs(actual - expected) <= tol) then
            print '(A, A)', " [PASS] ", desc
            p_count = p_count + 1
        else
            print '(A, A, A, ES14.6, A, ES14.6)', " [FAIL] ", desc, &
                  " (Actual: ", actual, ", Expected: ", expected, ")"
        end if
    end subroutine assert_real_equal

    subroutine assert_true(cond, desc, p_count, t_count)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: desc
        integer, intent(inout) :: p_count, t_count
        t_count = t_count + 1
        if (cond) then
            print '(A, A)', " [PASS] ", desc
            p_count = p_count + 1
        else
            print '(A, A)', " [FAIL] ", desc
        end if
    end subroutine assert_true

end program test_dipolar_scattering
