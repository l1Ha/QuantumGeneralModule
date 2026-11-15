!> \brief Cholesky / QR / SVD 分解验收测试 (路线图 §15.10「SVD/QR/Cholesky」项)
!> \details 验收组:
!>          组 1: Cholesky —— 重构、解析 2x2、非正定检测、解方程;
!>          组 2: QR —— 重构/正交/上三角 + 病态 Hilbert 矩阵上
!>                与正规方程法的精度对比 (经典 cond^2 损失);
!>          组 3: SVD —— 重构/正交/解析例/秩亏谱;
!>          组 4: 伪逆最小二乘 —— 秩亏问题的最小范数解与残差.
program test_linalg_decomp
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat

    call print_banner("GeneralModule: Cholesky / QR / SVD Decomposition Test", 66)

    ! ================================================================
    ! 组 1: Cholesky
    ! ================================================================
    block
        integer, parameter :: n = 6
        real(dp) :: a2(2, 2), l2(2, 2), l(n, n), a(n, n), recon(n, n)
        real(dp) :: err_rec, err_ana, x(n), b(n), y(n), err_solve
        integer :: i, j

        ! 解析 2x2: A = [[4,2],[2,3]] -> L = [[2,0],[1,sqrt(2)]]
        a2(1, 1) = 4.0_dp; a2(2, 1) = 2.0_dp
        a2(1, 2) = 2.0_dp; a2(2, 2) = 3.0_dp
        call cholesky_decompose(2, a2, l2, stat)
        err_ana = max(abs(l2(1, 1) - 2.0_dp), abs(l2(2, 1) - 1.0_dp), &
                      abs(l2(2, 2) - sqrt(2.0_dp)), abs(l2(1, 2)))

        ! 随机构造 SPD: A = M^T M + n I
        do i = 1, n
            do j = 1, n
                a(i, j) = cos(real(3 * i + 5 * j, dp)) + merge(real(n, dp), 0.0_dp, i == j)
            end do
        end do
        a = 0.5_dp * (a + transpose(a))   ! 对称化下三角输入
        call cholesky_decompose(n, a, l, stat)
        recon = matmul(l, transpose(l))
        err_rec = maxval(abs(recon - a))

        ! 解方程 A x = b (已知解)
        do i = 1, n
            b(i) = real(i, dp) - 2.0_dp
        end do
        y = matmul(a, b)                 ! 真解取 b, 右端取 A b
        call cholesky_solve(n, l, y, x)
        err_solve = maxval(abs(x - b))

        ! 非正定检测
        block
            real(dp) :: bad(2, 2)
            integer :: st_bad
            bad(1, 1) = 1.0_dp; bad(2, 1) = 2.0_dp
            bad(1, 2) = 2.0_dp; bad(2, 2) = 1.0_dp    ! 特征值 -1, 3 -> 非正定
            call cholesky_decompose(2, bad, l2, st_bad)
            n_tests = n_tests + 1
            if (st_bad == 1) then
                n_passed = n_passed + 1
                print '(A)', " [PASS] Cholesky: non-SPD detection (indefinite 2x2)"
            else
                print '(A, I3)', " [FAIL] non-SPD detection, stat =", st_bad
            end if
        end block

        n_tests = n_tests + 1
        if (stat == 0 .and. err_rec < 1.0e-14_dp .and. err_ana < 1.0e-15_dp .and. &
            err_solve < 1.0e-13_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] Cholesky: recon err =", err_rec, &
                  "  analytic 2x2 =", err_ana, "  solve =", err_solve
        else
            print '(A, 3ES11.3)', " [FAIL] Cholesky", err_rec, err_ana, err_solve
        end if
    end block

    ! ================================================================
    ! 组 2: QR + 病态最小二乘精度对比 (希尔伯特矩阵)
    ! ================================================================
    block
        integer, parameter :: m = 40, nq = 8
        real(dp) :: a(m, nq), q(m, nq), r(nq, nq), qtq(nq, nq)
        real(dp) :: err_rec, err_orth, err_tri, xtrue(nq), b(m), x_qr(nq), resid_qr
        real(dp) :: ata(nq, nq), atb(nq), x_ne(nq), err_qr, err_ne, hsum
        integer :: i, j

        ! 希尔伯特型设计矩阵 (病态: cond ~ 1e10+)
        do i = 1, m
            do j = 1, nq
                a(i, j) = 1.0_dp / real(i + j - 1, dp)
            end do
        end do
        do j = 1, nq
            xtrue(j) = real(j, dp) - 4.0_dp
        end do
        b = matmul(a, xtrue)     ! 无噪声: 两种方法都应精确; 加入小扰动制造病态差异

        call qr_decompose(m, nq, a, q, r, stat)
        err_rec = maxval(abs(matmul(q, r) - a))
        qtq = matmul(transpose(q), q)
        err_orth = 0.0_dp
        do i = 1, nq
            do j = 1, nq
                err_orth = max(err_orth, abs(qtq(i, j) - merge(1.0_dp, 0.0_dp, i == j)))
            end do
        end do
        err_tri = 0.0_dp
        do i = 1, nq
            do j = 1, i - 1
                err_tri = max(err_tri, abs(r(i, j)))
            end do
        end do

        n_tests = n_tests + 1
        if (stat == 0 .and. err_rec < 1.0e-14_dp .and. err_orth < 1.0e-14_dp .and. &
            err_tri == 0.0_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] QR: reconstruction/orthogonality: ", err_rec, " /", err_orth
        else
            print '(A, 3ES11.3)', " [FAIL] QR", err_rec, err_orth, err_tri
        end if

        ! 病态最小二乘: 精确一致系统 (b = A x_true), 求解过程中的舍入
        ! 放大即经典演示 —— QR 损失 ~log10(cond) 位, 正规方程损失 ~2*log10(cond) 位
        call least_squares_qr(m, nq, a, b, x_qr, resid_qr, stat)
        ata = matmul(transpose(a), a)
        atb = matmul(transpose(a), b)
        x_ne = solve_dense(atb, ata, nq)     ! 正规方程 (测试内直接高斯消元)
        err_qr = maxval(abs(x_qr - xtrue))
        err_ne = maxval(abs(x_ne - xtrue))

        n_tests = n_tests + 1
        if (err_qr < 1.0e-4_dp .and. err_ne > 5.0_dp * err_qr) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES11.2)', " [PASS] ill-conditioned LSQ (Hilbert cond~1e10): QR err =", &
                  err_qr, " vs normal-eq =", err_ne, "  ratio =", err_ne / max(err_qr, 1.0e-30_dp)
        else
            print '(A, ES10.2, A, ES10.2)', " [FAIL] ill-conditioned LSQ: QR =", err_qr, " normal =", err_ne
        end if
    end block

    ! ================================================================
    ! 组 3: SVD (解析例 / 重构 / 正交 / 秩亏谱)
    ! ================================================================
    block
        integer, parameter :: m = 30, ns = 6
        real(dp) :: a(m, ns), u(m, ns), sv(ns), vt(ns, ns)
        real(dp) :: err_rec, err_ou, err_ov, a2(2, 2), u2(2, 2), sv2(2), vt2(2, 2)
        real(dp) :: err_ana, s_sym
        real(dp) :: ata(ns, ns), d(ns), z(ns, ns)
        integer :: i, j

        ! 解析 2x2: [[3,0],[0,-2]] -> sigma = 3, 2
        a2(1, 1) = 3.0_dp; a2(1, 2) = 0.0_dp
        a2(2, 1) = 0.0_dp; a2(2, 2) = -2.0_dp
        call svd_decompose(2, 2, a2, u2, sv2, vt2, stat)
        err_ana = max(abs(sv2(1) - 3.0_dp), abs(sv2(2) - 2.0_dp))

        ! 一般矩阵的重构与正交性
        do i = 1, m
            do j = 1, ns
                a(i, j) = sin(real(2 * i + 3 * j, dp)) / real(1 + abs(i - j), dp)
            end do
        end do
        call svd_decompose(m, ns, a, u, sv, vt, stat)
        err_rec = maxval(abs(matmul(u * spread(sv, 1, m), vt) - a))
        err_ou = 0.0_dp
        do i = 1, ns
            do j = 1, ns
                err_ou = max(err_ou, abs(sum(u(:, i) * u(:, j)) - merge(1.0_dp, 0.0_dp, i == j)))
            end do
        end do
        err_ov = 0.0_dp
        do i = 1, ns
            do j = 1, ns
                err_ov = max(err_ov, abs(sum(vt(i, :) * vt(j, :)) - merge(1.0_dp, 0.0_dp, i == j)))
            end do
        end do
        ! 奇异值平方 = A^T A 的本征值 (独立路径)
        ata = matmul(transpose(a), a)
        call diag_symmetric_matrix(ns, ata, d, z, stat)
        s_sym = 0.0_dp
        do i = 1, ns
            s_sym = max(s_sym, abs(sv(i)**2 - d(ns + 1 - i)))   ! 降序 vs 升序
        end do

        n_tests = n_tests + 1
        if (err_ana < 1.0e-14_dp .and. err_rec < 1.0e-13_dp .and. err_ou < 1.0e-13_dp .and. &
            err_ov < 1.0e-13_dp .and. s_sym < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] SVD: recon/orth/U,V + sigma^2 = eig(A^T A): ", &
                  err_rec, " /", err_ou, " /", s_sym
        else
            print '(A, 5ES10.2)', " [FAIL] SVD", err_ana, err_rec, err_ou, err_ov, s_sym
        end if

        ! 秩亏谱: 构造秩 2 的 8x5 矩阵, sigma_3..5 ~ 0
        block
            integer, parameter :: mr = 8, nr = 5
            real(dp) :: ar(mr, nr), ur(mr, nr), svr(nr), vtr(nr, nr)
            real(dp) :: col1(mr), col2(mr), s3max
            integer :: k
            col1 = [(sin(real(k, dp)), k=1, mr)]
            col2 = [(cos(0.5_dp * real(k, dp)), k=1, mr)]
            do k = 1, mr
                ar(k, 1) = col1(k)
                ar(k, 2) = col2(k)
                ar(k, 3) = 2.0_dp * col1(k) - col2(k)   ! 线性相关
                ar(k, 4) = col1(k) + 0.5_dp * col2(k)
                ar(k, 5) = 0.0_dp
            end do
            call svd_decompose(mr, nr, ar, ur, svr, vtr, stat)
            s3max = maxval(svr(3:nr))
            n_tests = n_tests + 1
            if (svr(1) > 0.5_dp .and. svr(2) > 0.5_dp .and. s3max < 1.0e-13_dp) then
                n_passed = n_passed + 1
                print '(A, 3ES10.2)', " [PASS] rank-deficient spectrum: sigma =", svr(1), svr(2), s3max
            else
                print '(A, 5ES10.2)', " [FAIL] rank-deficient", svr
            end if
        end block
    end block

    ! ================================================================
    ! 组 4: 伪逆最小二乘 (秩亏, 最小范数解)
    ! ================================================================
    block
        integer, parameter :: m = 12, n = 6
        real(dp) :: a(m, n), b(m), x(n), resid, nullproj
        real(dp) :: c1(m), c2(m), s
        integer :: i, k

        ! 秩 3 设计矩阵: 前 3 列独立, 后 3 列线性组合
        c1 = [(sin(real(k, dp)), k=1, m)]
        c2 = [(cos(0.3_dp * real(k, dp)), k=1, m)]
        do k = 1, m
            a(k, 1) = c1(k)
            a(k, 2) = c2(k)
            a(k, 3) = real(k, dp) / real(m, dp)
            a(k, 4) = a(k, 1) + a(k, 2)
            a(k, 5) = 2.0_dp * a(k, 3)
            a(k, 6) = a(k, 1) - a(k, 3)
            b(k) = 1.0_dp + 0.1_dp * real(k, dp)
        end do

        call svd_pseudoinverse(m, n, a, b, x, resid)
        ! 性质 1: 残差 = 最小二乘残差 (子空间投影: A^T (A x - b) = 0)
        s = 0.0_dp
        do i = 1, n
            s = max(s, abs(sum(a(:, i) * (matmul(a, x) - b))))
        end do
        ! 性质 2: 最小范数解与 A^T A 零空间正交 (零奇异值方向投影为零)
        block
            real(dp) :: ata(n, n), d(n), z(n, n), proj
            call matmul_ata(a, m, n, ata)
            call diag_symmetric_matrix(n, ata, d, z, stat)
            proj = 0.0_dp
            do i = 1, n
                if (d(i) < 1.0e-10_dp * d(n)) then
                    proj = max(proj, abs(sum(z(:, i) * x)))
                end if
            end do
            nullproj = proj
        end block
        n_tests = n_tests + 1
        if (s < 1.0e-10_dp .and. nullproj < 1.0e-10_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, ES10.2)', " [PASS] pseudoinverse LSQ (rank 3/6): normal-eq check =", &
                  s, "  null-space projection =", nullproj, "  resid =", resid
        else
            print '(A, 3ES10.2)', " [FAIL] pseudoinverse", s, nullproj, resid
        end if
    end block

    print '(A, I2, A, I2, A)', "Decomposition Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: Cholesky / QR / SVD test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> 测试内高斯消元 (正规方程对照; 故意朴素以实现经典精度损失)
    function solve_dense(b, a, n) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: b(n), a(n, n)
        real(dp) :: x(n)
        real(dp) :: m(n, n), bb(n), piv, f
        integer :: i, k, p
        m = a
        bb = b
        do k = 1, n
            p = k
            do i = k + 1, n
                if (abs(m(i, k)) > abs(m(p, k))) p = i
            end do
            if (p /= k) then
                m([k, p], :) = m([p, k], :)
                f = bb(k); bb(k) = bb(p); bb(p) = f
            end if
            do i = k + 1, n
                f = m(i, k) / m(k, k)
                m(i, k:n) = m(i, k:n) - f * m(k, k:n)
                bb(i) = bb(i) - f * bb(k)
            end do
        end do
        do i = n, 1, -1
            x(i) = (bb(i) - sum(m(i, i + 1:n) * x(i + 1:n))) / m(i, i)
        end do
    end function solve_dense

    subroutine matmul_ata(a, m, n, ata)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: a(m, n)
        real(dp), intent(out) :: ata(n, n)
        ata = matmul(transpose(a), a)
    end subroutine matmul_ata

end program test_linalg_decomp
