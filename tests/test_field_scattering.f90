!> \brief 外加磁场/电场超冷量子散射与多基组模块自动化单元测试
!> \details 涵盖半整数角动量代数、Breit-Rabi 解析与数值对比、四大基组严格幺正变换、
!>          电子自旋交换势、全耦合哈密顿量与磁 Feshbach 共振色散扫描拟合。
!> \author LiHao
program test_field_scattering
    use mod_constants, only: dp, PI, TWOPI, GAUSS2AU
    use mod_special_functions, only: clebsch_gordan_half, wigner_3j_half, wigner_6j_half, wigner_9j_half
    use mod_field_scattering
    implicit none

    integer :: total_tests = 0
    integer :: passed_tests = 0

    print *, "================================================================"
    print *, "    TEST SUITE: Field Scattering & Multi-Basis Coupled Channels "
    print *, "================================================================"

    call test_half_integer_angular_momentum()
    call test_breit_rabi_energies()
    call test_four_bases_and_unitary_transforms()
    call test_spin_exchange_matrix()
    call test_feshbach_resonance_scan()
    call test_heteronuclear_scattering()

    print *, ""
    print *, "================================================================"
    write(*, '(A, I3, A, I3, A)') " Field Scattering Test Results: ", passed_tests, "/", total_tests, " PASSED"
    print *, "================================================================"

    if (passed_tests /= total_tests) then
        error stop "Some tests failed in test_field_scattering!"
    end if

contains

    subroutine assert_true(cond, name)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        if (cond) then
            passed_tests = passed_tests + 1
            write(*, '(A, A, A)') "  [PASS] ", trim(name), ""
        else
            write(*, '(A, A, A)') "  [FAIL] ", trim(name), ""
        end if
    end subroutine assert_true

    subroutine assert_near(val, target, tol, name)
        real(dp), intent(in) :: val, target, tol
        character(len=*), intent(in) :: name
        real(dp) :: diff
        diff = abs(val - target)
        call assert_true(diff <= tol, name)
        if (diff > tol) then
            write(*, '(A, ES14.6, A, ES14.6, A, ES14.6)') &
                "         Got: ", val, " Expected: ", target, " Diff: ", diff
        end if
    end subroutine assert_near

    ! --------------------------------------------------------------------------
    ! 测试 1: 半整数角动量代数 (Clebsch-Gordan, Wigner 3j, 6j, 9j)
    ! --------------------------------------------------------------------------
    subroutine test_half_integer_angular_momentum()
        real(dp) :: cg, w3j, w6j, w9j
        print *, ""
        print *, "--- 1. Testing Half-Integer Angular Momentum Algebra ---"

        ! <1/2, 1/2; 1/2, -1/2 | 1, 0> = 1/sqrt(2)
        cg = clebsch_gordan_half(1, 1, 1, -1, 2, 0)
        call assert_near(cg, 1.0_dp / sqrt(2.0_dp), 1.0e-14_dp, "CG <1/2, 1/2; 1/2, -1/2 | 1, 0>")

        ! <1/2, 1/2; 1/2, -1/2 | 0, 0> = 1/sqrt(2)
        cg = clebsch_gordan_half(1, 1, 1, -1, 0, 0)
        call assert_near(cg, 1.0_dp / sqrt(2.0_dp), 1.0e-14_dp, "CG <1/2, 1/2; 1/2, -1/2 | 0, 0>")

        ! (1/2 1/2 1; 1/2 -1/2 0) = 1/sqrt(6)
        w3j = wigner_3j_half(1, 1, 2, 1, -1, 0)
        call assert_near(w3j, 1.0_dp / sqrt(6.0_dp), 1.0e-14_dp, "Wigner 3j (1/2, 1/2, 1; 1/2, -1/2, 0)")

        ! Wigner 6j: {1/2 1/2 1; 1/2 1/2 0} = (-1)^(1/2+1/2+1) / sqrt((2*1/2+1)*(2*1/2+1)) = +1 / 2 = 0.5
        w6j = wigner_6j_half(1, 1, 2, 1, 1, 0)
        call assert_near(w6j, 0.5_dp, 1.0e-14_dp, "Wigner 6j {1/2, 1/2, 1; 1/2, 1/2, 0} = 0.5")

        ! Wigner 9j 恒等式测试: 包含 0 时退化为 6j: (-1)^3 / (3) * {1/2 1/2 1; 1/2 1/2 1} = -1/3 * (1/6) = -1/18
        w9j = wigner_9j_half(1, 1, 2, 1, 1, 2, 2, 2, 0)
        call assert_near(w9j, -1.0_dp / 18.0_dp, 1.0e-14_dp, "Wigner 9j with zero element = -1/18")
    end subroutine test_half_integer_angular_momentum

    ! --------------------------------------------------------------------------
    ! 测试 2: 单原子 Breit-Rabi 塞曼-超精细能级与解析公式对照
    ! --------------------------------------------------------------------------
    subroutine test_breit_rabi_energies()
        type(cold_atom_t) :: rb87
        real(dp) :: b_gauss, b_au, delta_hfs_au, x, e_analytic_f2_m1, e_analytic_f1_m1
        real(dp), allocatable :: e_levels(:), psi_vecs(:, :)
        integer  :: n_states
        real(dp) :: split_b0, min_diff_f2, min_diff_f1

        print *, ""
        print *, "--- 2. Testing Breit-Rabi Zeeman-Hyperfine Spectrum (Rb87) ---"

        call get_cold_atom_preset("87Rb", rb87)
        call assert_true(rb87%two_s == 1 .and. rb87%two_i == 3, "Rb87 preset has s=1/2, i=3/2")

        ! 1. 零磁场下的超精细分裂: Delta_E = a_hf * (i + 1/2) = a_hf * 2 = 6.83468261 GHz
        call calc_breit_rabi_energies(rb87, 0.0_dp, e_levels, psi_vecs, n_states)
        call assert_true(n_states == 8, "Rb87 has 8 Zeeman states (s=1/2, i=3/2)")
        split_b0 = (e_levels(4) - e_levels(1)) * AU2GHZ
        call assert_near(split_b0, 6.83468261_dp, 1.0e-5_dp, "Zero-field hyperfine splitting = 6.83468 GHz")

        ! 2. 有限磁场 B = 50.0 Gauss 对比 Breit-Rabi 解析公式
        b_gauss = 50.0_dp
        b_au = b_gauss * GAUSS2AU
        call calc_breit_rabi_energies(rb87, b_gauss, e_levels, psi_vecs, n_states)

        delta_hfs_au = rb87%a_hf_ghz * GHZ2AU * 2.0_dp ! (I + 1/2) = 2
        x = (rb87%g_s * MU_B_AU - rb87%g_i * MU_N_AU) * b_au / delta_hfs_au

        ! 对于 m_F = +1:
        e_analytic_f2_m1 = -delta_hfs_au / 8.0_dp + rb87%g_i * MU_N_AU * b_au * 1.0_dp + &
                           0.5_dp * delta_hfs_au * sqrt(1.0_dp + (4.0_dp * 1.0_dp / 4.0_dp) * x + x**2)
        e_analytic_f1_m1 = -delta_hfs_au / 8.0_dp + rb87%g_i * MU_N_AU * b_au * 1.0_dp - &
                           0.5_dp * delta_hfs_au * sqrt(1.0_dp + (4.0_dp * 1.0_dp / 4.0_dp) * x + x**2)

        min_diff_f2 = minval(abs(e_levels - e_analytic_f2_m1))
        min_diff_f1 = minval(abs(e_levels - e_analytic_f1_m1))
        call assert_near(min_diff_f2, 0.0_dp, 1.0e-12_dp, "Rb87 |F=2, mF=1> Breit-Rabi exact match")
        call assert_near(min_diff_f1, 0.0_dp, 1.0e-12_dp, "Rb87 |F=1, mF=1> Breit-Rabi exact match")

        deallocate(e_levels, psi_vecs)
    end subroutine test_breit_rabi_energies

    ! --------------------------------------------------------------------------
    ! 测试 3: 四大基组构建与幺正变换矩阵性质
    ! --------------------------------------------------------------------------
    subroutine test_four_bases_and_unitary_transforms()
        type(cold_atom_t) :: rb87
        type(field_channel_t), allocatable :: ch_uncoupled(:), ch_f_coupled(:), ch_total_spin(:)
        real(dp), allocatable :: U_unc_f(:, :), U_unc_spin(:, :), U_unc_dress(:, :)
        real(dp), allocatable :: U_prod(:, :), I_ref(:, :)
        integer  :: n_ch, i
        real(dp) :: max_err, b_gauss

        print *, ""
        print *, "--- 3. Testing Four Representation Bases & Unitary Transformations ---"

        call get_cold_atom_preset("87Rb", rb87)
        b_gauss = 100.0_dp

        ! 构建总磁量子数 2*M_tot = 2 (即 M_tot = +1), s-波 (l_max=0) 碰撞子空间
        call build_field_collision_channels(rb87, rb87, BASIS_UNCOUPLED, 2, 0, ch_uncoupled, n_ch)
        call build_field_collision_channels(rb87, rb87, BASIS_F_COUPLED, 2, 0, ch_f_coupled, n_ch)
        call build_field_collision_channels(rb87, rb87, BASIS_TOTAL_SPIN, 2, 0, ch_total_spin, n_ch)

        call assert_true(n_ch > 0, "Channels constructed successfully (n_ch > 0)")

        allocate(U_unc_f(n_ch, n_ch), U_unc_spin(n_ch, n_ch), U_unc_dress(n_ch, n_ch))
        allocate(U_prod(n_ch, n_ch), I_ref(n_ch, n_ch))

        I_ref = 0.0_dp
        do i = 1, n_ch
            I_ref(i, i) = 1.0_dp
        end do

        ! (1) 非耦合基组 -> f-耦合基组变换矩阵 U_{f, unc} 幺正性检验
        call calc_basis_transform_matrix(rb87, rb87, ch_uncoupled, ch_f_coupled, &
                                         n_ch, BASIS_UNCOUPLED, BASIS_F_COUPLED, b_gauss, U_unc_f)
        U_prod = matmul(U_unc_f, transpose(U_unc_f))
        max_err = maxval(abs(U_prod - I_ref))
        call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Unitary: U(f <- unc) * U^T = I")

        ! (2) 非耦合基组 -> 总自旋耦合基组变换矩阵 U_{(S,I)F, unc} 幺正性检验
        call calc_basis_transform_matrix(rb87, rb87, ch_uncoupled, ch_total_spin, &
                                         n_ch, BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, b_gauss, U_unc_spin)
        U_prod = matmul(U_unc_spin, transpose(U_unc_spin))
        max_err = maxval(abs(U_prod - I_ref))
        call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Unitary: U((S,I)F <- unc) * U^T = I")

        ! (3) 非耦合基组 -> 场缀饰本征通道基组变换矩阵 U_{dress, unc} 幺正性检验
        call calc_basis_transform_matrix(rb87, rb87, ch_uncoupled, ch_uncoupled, &
                                         n_ch, BASIS_UNCOUPLED, BASIS_FIELD_DRESSED, b_gauss, U_unc_dress)
        U_prod = matmul(U_unc_dress, transpose(U_unc_dress))
        max_err = maxval(abs(U_prod - I_ref))
        call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Unitary: U(dress <- unc) * U^T = I")

        ! (4) 链式复合变换一致性检验: U_{f, spin} = U_{f, unc} * U_{unc, spin}
        block
            real(dp) :: U_chain(n_ch, n_ch), U_direct(n_ch, n_ch)
            U_chain = matmul(U_unc_f, transpose(U_unc_spin))
            call calc_basis_transform_matrix(rb87, rb87, ch_total_spin, ch_f_coupled, &
                                             n_ch, BASIS_TOTAL_SPIN, BASIS_F_COUPLED, b_gauss, U_direct)
            max_err = maxval(abs(U_chain - U_direct))
            call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Chain rule: U(f <- spin) == U(f<-unc) * U(spin<-unc)^T")
        end block

        deallocate(ch_uncoupled, ch_f_coupled, ch_total_spin)
        deallocate(U_unc_f, U_unc_spin, U_unc_dress, U_prod, I_ref)
    end subroutine test_four_bases_and_unitary_transforms

    ! --------------------------------------------------------------------------
    ! 测试 4: 电子自旋交换矩阵 s1 . s2 算符性质
    ! --------------------------------------------------------------------------
    subroutine test_spin_exchange_matrix()
        type(cold_atom_t) :: rb87
        type(field_channel_t), allocatable :: ch_total_spin(:)
        real(dp), allocatable :: P_exc(:, :)
        integer  :: n_ch, i
        real(dp) :: s_val, expected

        print *, ""
        print *, "--- 4. Testing Spin Exchange Operator s1 . s2 Properties ---"

        call get_cold_atom_preset("87Rb", rb87)
        call build_field_collision_channels(rb87, rb87, BASIS_TOTAL_SPIN, 2, 0, ch_total_spin, n_ch)

        allocate(P_exc(n_ch, n_ch))
        call build_spin_exchange_matrix(ch_total_spin, n_ch, BASIS_TOTAL_SPIN, P_exc)

        ! 在总自旋基组下，s1 . s2 必须严格对角！
        ! 本征值: S=0 -> -3/4; S=1 -> +1/4
        call assert_true(all(abs(P_exc - diag_part(P_exc, n_ch)) < 1.0e-14_dp), &
                         "s1 . s2 is strictly diagonal in Total Spin Basis")

        do i = 1, n_ch
            s_val = real(ch_total_spin(i)%two_S, dp) / 2.0_dp
            expected = 0.5_dp * (s_val * (s_val + 1.0_dp) - 0.75_dp - 0.75_dp)
            call assert_near(P_exc(i, i), expected, 1.0e-14_dp, "s1 . s2 eigenvalue for channel")
        end do

        deallocate(ch_total_spin, P_exc)
    end subroutine test_spin_exchange_matrix

    ! --------------------------------------------------------------------------
    ! 测试 5: 磁 Feshbach 共振色散扫描与拟合
    ! --------------------------------------------------------------------------
    subroutine test_feshbach_resonance_scan()
        real(dp), allocatable :: b_test(:), a_test(:)
        real(dp) :: b0_fit, delta_b_fit, a_bg_fit
        integer  :: n_b, i

        print *, ""
        print *, "--- 5. Testing Magnetic Feshbach Resonance Fit ---"

        ! 1. 检验拟合算法对标准色散公式 a(B) = a_bg * (1 - \Delta B / (B - B0)) 的高保真还原
        n_b = 61
        allocate(b_test(n_b), a_test(n_b))
        do i = 1, n_b
            b_test(i) = 50.0_dp + real(i - 1, dp) * 1.0_dp ! 50 G 到 110 G
            ! 设置真实参数: B0 = 80.0 G, Delta_B = 5.0 G, a_bg = 100.0 a0
            a_test(i) = 100.0_dp * (1.0_dp - 5.0_dp / (b_test(i) - 80.0_dp))
        end do

        call fit_feshbach_resonance_parameters(b_test, a_test, n_b, b0_fit, delta_b_fit, a_bg_fit)

        call assert_near(b0_fit, 80.0_dp, 1.0_dp, "Fitted Resonance Position B0 ~ 80 G")
        call assert_near(delta_b_fit, 5.0_dp, 1.0_dp, "Fitted Resonance Width Delta_B ~ 5 G")
        call assert_near(a_bg_fit, 100.0_dp, 5.0_dp, "Fitted Background Length a_bg ~ 100 a0")

        deallocate(b_test, a_test)
    end subroutine test_feshbach_resonance_scan

    ! --------------------------------------------------------------------------
    ! 测试 6: 异核双原子超冷碰撞体系 (Heteronuclear: 40K + 87Rb, 6Li + 87Rb)
    ! --------------------------------------------------------------------------
    subroutine test_heteronuclear_scattering()
        type(cold_atom_t) :: k40, rb87, li6
        type(field_channel_t), allocatable :: ch_unc(:), ch_f(:), ch_spin(:)
        real(dp), allocatable :: U_f_unc(:, :), U_spin_unc(:, :), U_dress_unc(:, :)
        real(dp), allocatable :: U_prod(:, :), I_ref(:, :)
        real(dp) :: max_err, mu_krb, mu_lirb, b_gauss
        integer  :: n_ch, i

        print *, ""
        print *, "--- 6. Testing Heteronuclear Ultracold Scattering (40K + 87Rb) ---"

        call get_cold_atom_preset("40K", k40)
        call get_cold_atom_preset("87Rb", rb87)
        call get_cold_atom_preset("6Li", li6)

        call assert_true(k40%two_s == 1 .and. k40%two_i == 8, "40K has s=1/2, i=4")
        call assert_true(rb87%two_s == 1 .and. rb87%two_i == 3, "87Rb has s=1/2, i=3/2")
        call assert_true(li6%two_s == 1 .and. li6%two_i == 2, "6Li has s=1/2, i=1")

        ! 验证异核折合质量 \mu = m1 * m2 / (m1 + m2)
        mu_krb = (k40%mass_amu * rb87%mass_amu) / (k40%mass_amu + rb87%mass_amu)
        call assert_near(mu_krb, 27.38006_dp, 0.01_dp, "40K-87Rb reduced mass ~ 27.38 amu")

        mu_lirb = (li6%mass_amu * rb87%mass_amu) / (li6%mass_amu + rb87%mass_amu)
        call assert_near(mu_lirb, 5.62677_dp, 0.01_dp, "6Li-87Rb reduced mass ~ 5.63 amu")

        ! 针对 40K + 87Rb 异核碰撞体系测试：
        ! 40K 具有半整数总角动量 f1, 87Rb 具有整数总角动量 f2, 总磁量子数 M_tot 为半整数
        ! 设置 2*M_tot = -7 (即 M_tot = -7/2), s-波 (l=0)
        call build_field_collision_channels(k40, rb87, BASIS_UNCOUPLED, -7, 0, ch_unc, n_ch)
        call build_field_collision_channels(k40, rb87, BASIS_F_COUPLED, -7, 0, ch_f, n_ch)
        call build_field_collision_channels(k40, rb87, BASIS_TOTAL_SPIN, -7, 0, ch_spin, n_ch)

        call assert_true(n_ch == 12, "40K+87Rb 2*M_tot=-7 subspace dimension = 12 across all bases")

        allocate(U_f_unc(n_ch, n_ch), U_spin_unc(n_ch, n_ch), U_dress_unc(n_ch, n_ch))
        allocate(U_prod(n_ch, n_ch), I_ref(n_ch, n_ch))
        I_ref = 0.0_dp
        do i = 1, n_ch
            I_ref(i, i) = 1.0_dp
        end do

        b_gauss = 540.0_dp ! 40K-87Rb 著名 Feshbach 共振区附近

        ! (1) 异核体系 U(f <- unc) 幺正性
        call calc_basis_transform_matrix(k40, rb87, ch_unc, ch_f, n_ch, &
                                         BASIS_UNCOUPLED, BASIS_F_COUPLED, b_gauss, U_f_unc)
        U_prod = matmul(U_f_unc, transpose(U_f_unc))
        max_err = maxval(abs(U_prod - I_ref))
        call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Heteronuclear 40K-87Rb U(f <- unc) unitary")

        ! (2) 异核体系 U(spin <- unc) 幺正性
        call calc_basis_transform_matrix(k40, rb87, ch_unc, ch_spin, n_ch, &
                                         BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, b_gauss, U_spin_unc)
        U_prod = matmul(U_spin_unc, transpose(U_spin_unc))
        max_err = maxval(abs(U_prod - I_ref))
        call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Heteronuclear 40K-87Rb U(spin <- unc) unitary")

        ! (3) 异核体系 U(dress <- unc) 场缀饰基组幺正性
        call calc_basis_transform_matrix(k40, rb87, ch_unc, ch_unc, n_ch, &
                                         BASIS_UNCOUPLED, BASIS_FIELD_DRESSED, b_gauss, U_dress_unc)
        U_prod = matmul(U_dress_unc, transpose(U_dress_unc))
        max_err = maxval(abs(U_prod - I_ref))
        call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Heteronuclear 40K-87Rb U(dress <- unc) unitary at 540 G")

        ! (4) 异核体系链式变换一致性: U(f <- spin) == U(f <- unc) * U(spin <- unc)^T
        block
            real(dp) :: U_chain(n_ch, n_ch), U_direct(n_ch, n_ch)
            U_chain = matmul(U_f_unc, transpose(U_spin_unc))
            call calc_basis_transform_matrix(k40, rb87, ch_spin, ch_f, n_ch, &
                                             BASIS_TOTAL_SPIN, BASIS_F_COUPLED, b_gauss, U_direct)
            max_err = maxval(abs(U_chain - U_direct))
            call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Heteronuclear chain rule U(f <- spin) matches direct")
        end block

        ! (5) 异核体系总自旋基下的电子交换算符对角性
        block
            real(dp) :: P_exc(n_ch, n_ch)
            call build_spin_exchange_matrix(ch_spin, n_ch, BASIS_TOTAL_SPIN, P_exc)
            max_err = maxval(abs(P_exc - diag_part(P_exc, n_ch)))
            call assert_near(max_err, 0.0_dp, 1.0e-14_dp, "Heteronuclear s1 . s2 strictly diagonal in Total Spin")
        end block

        deallocate(ch_unc, ch_f, ch_spin, U_f_unc, U_spin_unc, U_dress_unc, U_prod, I_ref)
    end subroutine test_heteronuclear_scattering

    pure function diag_part(A, n) result(D)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n, n)
        real(dp) :: D(n, n)
        integer :: i
        D = 0.0_dp
        do i = 1, n
            D(i, i) = A(i, i)
        end do
    end function diag_part

end program test_field_scattering
