!> \brief 库仑三体系统与 Hylleraas-Pekeris 变分法单元测试
program test_coulomb_threebody
    use mod_constants, only: dp
    use mod_coulomb_threebody
    implicit none

    real(dp) :: uvw(3), r_vals(3), vol, vol_manual
    real(dp) :: int_exact, int_quad, e_he, e_ion, alpha, e_best
    real(dp) :: xq(24), wq(24)
    integer  :: n_basis, idx_l(256), idx_m(256), idx_n(256)
    integer  :: i, n_tests, n_pass, stat
    integer, parameter :: NMAX = 7

    n_tests = 0
    n_pass  = 0

    print '(A)', "=========================================================="
    print '(A)', "     GeneralModule Unit Tests: Coulomb Three-Body         "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Pekeris 坐标变换正逆一致性
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    uvw = (/ 3.1_dp, 1.2_dp, 0.9_dp /)
    call perkeris_coordinate_transform(uvw, r_vals)
    ! 逆变换: u = r1 + r2 - r12 等
    if (abs((r_vals(1) + r_vals(2) - r_vals(3)) - uvw(1)) < 1.0e-12_dp .and. &
        abs((r_vals(1) - r_vals(2) + r_vals(3)) - uvw(2)) < 1.0e-12_dp .and. &
        abs((-r_vals(1) + r_vals(2) + r_vals(3)) - uvw(3)) < 1.0e-12_dp) then
        print '(A, 3F8.4, A)', " [PASS] Pekeris transform roundtrip: r = (", r_vals, " )"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Pekeris transform roundtrip mismatch"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: 雅可比体积元显式公式一致性
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    vol = calc_perkeris_volume_element(uvw)
    vol_manual = (uvw(1) + uvw(2)) * (uvw(1) + uvw(3)) * (uvw(2) + uvw(3))
    if (abs(vol - vol_manual) < 1.0e-12_dp .and. vol > 0.0_dp) then
        print '(A, F10.4)', " [PASS] Pekeris volume element consistent: ", vol
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Pekeris volume element mismatch"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: 解析 Gamma 递推积分 vs Gauss-Laguerre 求积交叉验证
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call gauss_laguerre_nodes(24, xq, wq)
    ! 三维张量积独立计算: 每轴代换 x = kappa*coord, du = dx/kappa
    ! int u^a e^{-kappa u} du = (1/kappa) sum_j w_j (x_j/kappa)^a
    block
        real(dp) :: xu(24), wu(24), su, sv, sw
        integer  :: ia
        call gauss_laguerre_nodes(24, xu, wu)
        su = 0.0_dp; sv = 0.0_dp; sw = 0.0_dp
        do ia = 1, 24
            su = su + wu(ia) * (xu(ia) / 2.0_dp)**2   ! u^2, kappa = 2
            sv = sv + wu(ia) * (xu(ia) / 2.0_dp)**1   ! v^1
            sw = sw + wu(ia) * (xu(ia) / 2.0_dp)**3   ! w^3
        end do
        int_quad = (su / 2.0_dp) * (sv / 2.0_dp) * (sw / 2.0_dp)
    end block
    int_exact = calc_hylleraas_analytical_integral(2, 1, 3, 2.0_dp)
    if (abs(int_quad - int_exact) / int_exact < 1.0e-10_dp) then
        print '(A, ES10.3)', " [PASS] Analytical Gamma integral matches GL quadrature: rel_err = ", &
            abs(int_quad - int_exact) / int_exact
        n_pass = n_pass + 1
    else
        print '(A, ES10.3, A, ES10.3)', " [FAIL] Gamma integral mismatch: quad = ", int_quad, &
            " exact = ", int_exact
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: 氦原子基态变分能量收敛到 Pekeris 精确值 (变分原理: E >= E_exact)
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    e_best = 0.0_dp
    do i = 16, 24, 2
        alpha = 0.1_dp * real(i, dp)
        call solve_helium_ground_state_variational(NMAX, 2.0_dp, alpha, e_he, stat)
        if (stat /= 0) cycle
        if (e_best == 0.0_dp .or. e_he < e_best) e_best = e_he
    end do
    ! Pekeris (1958) 精确值: E = -2.9037243770321985 Hartree
    if (abs(e_best + 2.90372437703_dp) < 1.0e-6_dp) then
        print '(A, F14.9, A, ES9.2)', " [PASS] He ground state converged to Pekeris value: E = ", &
            e_best, "  err = ", abs(e_best + 2.90372437703_dp)
        n_pass = n_pass + 1
    else
        print '(A, F14.9)', " [FAIL] He ground state energy not converged: E = ", e_best
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: 高 Z 类氦离子标度律 E(Z) ~ -Z^2 + 5Z/8
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call solve_helium_ground_state_variational(5, 10.0_dp, 9.5_dp, e_ion, stat)
    ! Z=10 类氦离子精确非相对论基态: E = -93.90680628 Hartree
    ! (变分结果含二阶 1/Z 项 -0.15766, 远比一阶公式 -Z^2+5Z/8 = -93.75 精确)
    if (abs(e_ion - (-93.90680628_dp)) < 1.0e-3_dp) then
        print '(A, F12.6, A, ES9.2)', " [PASS] He-like ion Z=10 matches exact: E = ", e_ion, &
            "  err = ", abs(e_ion + 93.90680628_dp)
        n_pass = n_pass + 1
    else
        print '(A, F12.6)', " [FAIL] He-like ion Z=10 energy wrong: E = ", e_ion
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: 基组索引完备性 (组合数 C(n_max+3, 3))
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call hylleraas_basis_indexing(4, n_basis, idx_l, idx_m, idx_n)
    ! n_max = 4: C(7,3) = 35
    if (n_basis == 35 .and. idx_l(1) == 0 .and. idx_m(1) == 0 .and. idx_n(1) == 0) then
        print '(A, I4)', " [PASS] Hylleraas basis indexing complete: n_basis = ", n_basis
        n_pass = n_pass + 1
    else
        print '(A, I4)', " [FAIL] Basis indexing count wrong: ", n_basis
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  Coulomb Three-Body Tests: ", n_pass, " / ", n_tests, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_tests) error stop "Test failures detected in test_coulomb_threebody."
end program test_coulomb_threebody
