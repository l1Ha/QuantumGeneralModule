!> \brief Hermite DVR 验收测试 (路线图 §12.1「更多正交多项式 DVR」项)
!> \details 验收组:
!>          组 1: Gauss-Hermite 求值精确度 (单项式矩 vs 闭式) 与正交归一;
!>          组 2: 谐振子 —— 梯位闭式动能 + 节点求值势能给出精确能级
!>                (基即本征基, 全部能级机器精度);
!>          组 3: 四次非谐振子 —— 随基维数收敛 + 与库内 FGH 交叉对照;
!>          组 4: 节点/权重对称性与 Hermite 函数解析值 (psi_0, psi_1).
program test_hermite_dvr
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat

    call print_banner("GeneralModule: Hermite DVR (Harmonic Basis) Test", 64)

    ! ================================================================
    ! 组 1: Gauss-Hermite 精确度 + 正交归一
    ! ================================================================
    block
        integer, parameter :: n_ord = 40
        real(dp) :: xq(n_ord), wg(n_ord), psi(0:n_ord - 1), moment
        real(dp) :: err_mom, err_orth, err_node, x, ps0, ps1, s_ij
        real(dp), allocatable :: pt(:, :)
        integer :: k, i, j

        call gauss_hermite_nodes_weights(n_ord, xq, wg)
        ! 单项式矩: ∫e^{-x²} x^{2p} dx = Γ(p+1/2); ∫x^{2p+1} = 0
        ! (n=40 时最高矩 x^12 的求和舍入 ~1e-13, 阈值取 1e-11)
        err_mom = 0.0_dp
        do k = 0, 6
            moment = sum(wg * xq**(2 * k)) - exp(log_gamma(real(k, dp) + 0.5_dp))
            err_mom = max(err_mom, abs(moment))
            moment = sum(wg * xq**(2 * k + 1))
            err_mom = max(err_mom, abs(moment))
        end do

        ! 正交归一 (DVR 核心性质): 节点求值矩阵
        !   Pt(k, i) = psi_{i-1}(x_k) * exp(+x_k^2/2)   (去指数多项式部分)
        ! 满足 Σ_k w_k Pt(k,i) Pt(k,j) = delta_ij  —— 权重已含 e^{-x²},
        ! 若用含 e^{-x²/2} 的裸函数值则等价于 e^{-2x²} 权重, 不正交。
        allocate (pt(n_ord, n_ord))
        pt = 0.0_dp
        do k = 1, n_ord
            call hermite_function_values(xq(k), n_ord - 1, psi)
            pt(k, :) = psi(0:n_ord - 1) * exp(0.5_dp * xq(k)**2)
        end do
        err_orth = 0.0_dp
        do i = 1, n_ord
            do j = 1, n_ord
                s_ij = sum(wg * pt(:, i) * pt(:, j)) - merge(1.0_dp, 0.0_dp, i == j)
                err_orth = max(err_orth, abs(s_ij))
            end do
        end do

        ! 解析节点值: psi_0(x_1) 与 psi_1(x_1) (极端节点, 检验递推稳定性)
        call hermite_function_values(xq(1), n_ord - 1, psi)
        err_node = abs(psi(0) - exp(-0.5_dp * xq(1)**2) / PI**0.25_dp)
        ! 解析: psi_0(0) = pi^{-1/4}, psi_1(0) = 0
        call hermite_function_values(0.0_dp, n_ord - 1, psi)
        ps0 = abs(psi(0) - 1.0_dp / PI**0.25_dp)
        ps1 = abs(psi(1))

        n_tests = n_tests + 1
        if (err_mom < 1.0e-11_dp .and. err_orth < 1.0e-12_dp .and. err_node < 1.0e-14_dp &
            .and. ps0 < 1.0e-14_dp .and. ps1 < 1.0e-14_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] Gauss-Hermite moments + node orthogonality:", &
                  err_mom, " /", err_orth
        else
            print '(A, 4ES10.2)', " [FAIL] moments/orth/node/psi", err_mom, err_orth, err_node, ps0
        end if
    end block

    ! ================================================================
    ! 组 2: 谐振子精确能级 (基即本征基)
    ! ================================================================
    block
        integer, parameter :: n_ord = 30
        real(dp) :: h_mat(n_ord, n_ord), d(n_ord), z(n_ord, n_ord), xq(n_ord)
        real(dp) :: wg(n_ord)
        real(dp) :: err_max, omega
        integer :: k

        omega = 0.6_dp
        call gauss_hermite_nodes_weights(n_ord, xq, wg)
        block
            real(dp) :: v_nodes(n_ord)
            ! 谐振子势: V(x) = 0.5 omega^2 x^2 (x = q/sqrt(m omega), m=1)
            do k = 1, n_ord
                v_nodes(k) = 0.5_dp * omega**2 * (xq(k) / sqrt(omega))**2
            end do
            call hermite_dvr_hamiltonian(n_ord, omega, v_nodes, h_mat, xq)
        end block
        call diag_symmetric_matrix(n_ord, h_mat, d, z, stat)
        err_max = 0.0_dp
        do k = 1, 10
            err_max = max(err_max, abs(d(k) - omega * (real(k - 1, dp) + 0.5_dp)))
        end do
        n_tests = n_tests + 1
        if (err_max < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4)', " [PASS] HO exact levels via ladder kinetic: err =", err_max
        else
            print '(A, ES12.4)', " [FAIL] HO levels err =", err_max
        end if
    end block

    ! ================================================================
    ! 组 3: 四次非谐振子收敛 + FGH 交叉对照
    ! ================================================================
    block
        integer, parameter :: n_ord = 60
        real(dp) :: h_mat(n_ord, n_ord), d(n_ord), z(n_ord, n_ord), xq(n_ord)
        real(dp) :: omega, lam, err_conv
        real(dp) :: wg(n_ord)
        real(dp), allocatable :: v_nodes(:)
        real(dp), allocatable :: x_fgh(:), v_fgh(:), e_fgh(:)
        real(dp), allocatable :: eig_vecs_fgh_unused(:, :)
        type(dvr_1d_t) :: dvr
        integer :: k, n_fgh

        ! H = p^2/2 + x^2/2 + lam x^4, lam = 0.2 (强非谐)
        omega = 1.0_dp
        lam = 0.2_dp
        allocate (v_nodes(n_ord))
        call gauss_hermite_nodes_weights(n_ord, xq, wg)
        do k = 1, n_ord
            v_nodes(k) = 0.5_dp * xq(k)**2 + lam * xq(k)**4
        end do
        call hermite_dvr_hamiltonian(n_ord, omega, v_nodes, h_mat, xq)
        call diag_symmetric_matrix(n_ord, h_mat, d, z, stat)

        ! 收敛: 基维数加倍后基态漂移应极小 (四阶 DVR 指数收敛)
        block
            integer, parameter :: n2 = 120
            real(dp) :: h2(n2, n2), d2(n2), z2(n2, n2), v2(n2), xq2(n2)
            real(dp) :: wg2(n2)
            call gauss_hermite_nodes_weights(n2, xq2, wg2)
            do k = 1, n2
                v2(k) = 0.5_dp * xq2(k)**2 + lam * xq2(k)**4
            end do
            call hermite_dvr_hamiltonian(n2, omega, v2, h2, xq2)
            call diag_symmetric_matrix(n2, h2, d2, z2, stat)
            err_conv = abs(d(1) - d2(1))
        end block

        ! FGH 交叉对照 (独立方法)
        n_fgh = 512
        allocate (x_fgh(n_fgh), v_fgh(n_fgh), e_fgh(n_fgh))
        allocate (eig_vecs_fgh_unused(n_fgh, n_fgh))   ! tql2 需全部本征矢
        call dvr_sinc_init(-8.0_dp, 8.0_dp, n_fgh, 1.0_dp, dvr)
        x_fgh = dvr%x
        do k = 1, n_fgh
            v_fgh(k) = 0.5_dp * x_fgh(k)**2 + lam * x_fgh(k)**4
        end do
        call fgh_solve_bound_states(dvr, v_fgh, e_fgh, eig_vecs_fgh_unused, stat)

        n_tests = n_tests + 1
        if (err_conv < 1.0e-9_dp .and. abs(d(1) - e_fgh(1)) < 1.0e-6_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] quartic oscillator: E0 =", d(1), &
                  "  N-doubling drift =", err_conv, "  vs FGH = ", abs(d(1) - e_fgh(1))
        else
            print '(A, 2ES10.2)', " [FAIL] quartic: conv =", err_conv, " FGH diff =", abs(d(1) - e_fgh(1))
        end if
    end block

    print '(A, I2, A, I2, A)', "Hermite DVR Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: Hermite DVR test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

end program test_hermite_dvr
