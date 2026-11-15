!> \brief 复 Hermitian 矩阵本征求解器验收测试 (§15.10/§15.11 路线图第 1 项)
!> \details 闭式与不变量验收:
!>          组 1: Pauli 矩阵 sigma_x, sigma_y, sigma_z 与一般 2x2 Hermitian
!>                解析本征值 (a+d)/2 +- sqrt(((a-d)/2)^2 + |b|^2);
!>          组 2: 实嵌入一致性 —— 实对称矩阵经复路径应与 diag_symmetric_matrix
!>                逐位一致;
!>          组 3: 随机构造 8x8 Hermitian (A = M + M^dag): 残差 ||Av - lam v||、
!>                正交归一 V^dag V = I、本征值全实、迹/行列式守恒;
!>          组 4: 物理算例 —— 三格点 Peierls 环 (均匀磁通) 闭式谱
!>                2 cos(2 pi (n + phi)/3) 与二能级复耦合 (圆偏振驱动) +-|Omega|/2;
!>          组 5: 简并稳健性 (二重/三重简并且存在) 与相位约定。
program test_hermitian_eig
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat

    call print_banner("GeneralModule: Complex Hermitian Eigen Solver Test", 60)

    ! ================================================================
    ! 组 1: Pauli 矩阵 + 一般 2x2 Hermitian 解析本征值
    ! ================================================================
    block
        complex(dp) :: a2(2, 2), vec(2, 2)
        real(dp) :: d(2), err_max, lam_ana(2)
        real(dp) :: aa, dd, bmod
        integer :: k

        err_max = 0.0_dp
        ! sigma_x
        a2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        a2(1, 2) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        a2(2, 1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        call diag_hermitian_matrix(2, a2, d, vec, stat)
        err_max = max(err_max, abs(d(1) + 1.0_dp), abs(d(2) - 1.0_dp))

        ! sigma_y = [[0, -i],[i, 0]] 本征值 -1, +1
        a2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        a2(1, 2) = cmplx(0.0_dp, -1.0_dp, kind=dp)
        a2(2, 1) = cmplx(0.0_dp, 1.0_dp, kind=dp)
        call diag_hermitian_matrix(2, a2, d, vec, stat)
        err_max = max(err_max, abs(d(1) + 1.0_dp), abs(d(2) - 1.0_dp))

        ! sigma_z
        a2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        a2(1, 1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        a2(2, 2) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
        call diag_hermitian_matrix(2, a2, d, vec, stat)
        err_max = max(err_max, abs(d(1) + 1.0_dp), abs(d(2) - 1.0_dp))

        ! 一般 2x2: A = [[a, b],[b*, d]] 解析本征值
        aa = 0.7_dp; dd = -0.2_dp
        a2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        a2(1, 1) = cmplx(aa, 0.0_dp, kind=dp)
        a2(2, 2) = cmplx(dd, 0.0_dp, kind=dp)
        a2(1, 2) = cmplx(0.3_dp, 0.45_dp, kind=dp)
        a2(2, 1) = conjg(a2(1, 2))
        call diag_hermitian_matrix(2, a2, d, vec, stat)
        bmod = abs(a2(1, 2))
        lam_ana(1) = 0.5_dp * (aa + dd) - sqrt((0.5_dp * (aa - dd))**2 + bmod**2)
        lam_ana(2) = 0.5_dp * (aa + dd) + sqrt((0.5_dp * (aa - dd))**2 + bmod**2)
        err_max = max(err_max, abs(d(1) - lam_ana(1)), abs(d(2) - lam_ana(2)))

        ! 本征矢残差
        do k = 1, 2
            err_max = max(err_max, norm2_complex(matmul(a2, vec(:, k)) - d(k) * vec(:, k)))
        end do

        n_tests = n_tests + 1
        if (stat == 0 .and. err_max < 1.0e-13_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4)', " [PASS] Pauli + analytic 2x2 Hermitian: max err =", err_max
        else
            print '(A, ES12.4)', " [FAIL] 2x2 Hermitian analytic", err_max
        end if
    end block

    ! ================================================================
    ! 组 2: 实嵌入一致性 (复路径对实对称输入应与实路径一致)
    ! ================================================================
    block
        integer, parameter :: n8 = 8
        real(dp) :: ar(n8, n8), dr(n8), zr(n8, n8)
        complex(dp) :: ac(n8, n8), zc(n8, n8)
        real(dp) :: dc(n8), err_ev, err_vec
        integer :: i, j, k

        do i = 1, n8
            do j = 1, n8
                ar(i, j) = cos(real(3 * i + 5 * j, dp)) / real(1 + abs(i - j), dp)
            end do
        end do
        ar = 0.5_dp * (ar + transpose(ar))
        ac = cmplx(ar, 0.0_dp, kind=dp)

        call diag_symmetric_matrix(n8, ar, dr, zr, stat)
        call diag_hermitian_matrix(n8, ac, dc, zc, stat)

        err_ev = maxval(abs(dr - dc))
        ! 本征矢允许整体相位差: 比较 |<v_r|v_c>| = 1
        err_vec = 0.0_dp
        do k = 1, n8
            err_vec = max(err_vec, 1.0_dp - abs(sum(conjg(cmplx(zr(:, k), 0.0_dp, kind=dp)) * zc(:, k))))
        end do

        n_tests = n_tests + 1
        if (err_ev < 1.0e-13_dp .and. err_vec < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4, A, ES12.4)', " [PASS] real-embedding consistency: dE =", err_ev, &
                  " 1-|<v_r|v_c>| =", err_vec
        else
            print '(A, ES12.4, A, ES12.4)', " [FAIL] real-embedding", err_ev, err_vec
        end if
    end block

    ! ================================================================
    ! 组 3: 随机 8x8 Hermitian: 残差 / 正交归一 / 全实 / 迹与行列式守恒
    ! ================================================================
    block
        integer, parameter :: n8 = 8
        complex(dp) :: a8(n8, n8), v8(n8, n8), resid(n8, n8)
        real(dp) :: d8(n8), err_res, err_orth, tr_a, tr_e, det_a, det_e, err_det
        integer :: i, j, k

        ! 构造确定性的"类随机" Hermitian: A = M + M^dag
        do i = 1, n8
            do j = 1, n8
                a8(i, j) = cmplx(sin(real(7 * i + 3 * j, dp)), cos(real(2 * i - 5 * j, dp)), kind=dp)
            end do
        end do
        a8 = a8 + conjg(transpose(a8))

        call diag_hermitian_matrix(n8, a8, d8, v8, stat)

        ! 残差 (逐列)
        err_res = 0.0_dp
        do k = 1, n8
            err_res = max(err_res, norm2_complex(matmul(a8, v8(:, k)) - d8(k) * v8(:, k)))
        end do

        ! 正交归一 V^dag V = I
        resid = matmul(conjg(transpose(v8)), v8)
        err_orth = 0.0_dp
        do i = 1, n8
            do j = 1, n8
                err_orth = max(err_orth, abs(resid(i, j) - merge(1.0_dp, 0.0_dp, i == j)))
            end do
        end do

        ! 迹与行列式守恒
        tr_a = real(sum([(a8(i, i), i=1, n8)]), dp)
        tr_e = sum(d8)
        det_a = abs(det_complex(a8, n8))
        det_e = abs(product(d8))
        err_det = abs(det_a - det_e) / max(1.0_dp, det_a)

        n_tests = n_tests + 1
        if (stat == 0 .and. err_res < 1.0e-12_dp .and. err_orth < 1.0e-13_dp .and. &
            abs(tr_a - tr_e) < 1.0e-12_dp .and. err_det < 1.0e-10_dp) then
            n_passed = n_passed + 1
            print '(A, ES11.3, A, ES11.3, A, ES11.3, A, ES11.3)', &
                  " [PASS] 8x8 Hermitian: resid =", err_res, &
                  " orth =", err_orth, " |dTr| =", abs(tr_a - tr_e), " |dDet| =", err_det
        else
            print '(A, 4ES11.3)', " [FAIL] 8x8 Hermitian", err_res, err_orth, &
                  abs(tr_a - tr_e), err_det
        end if
    end block

    ! ================================================================
    ! 组 4: 物理算例 —— Peierls 环闭式谱 + 复耦合二能级
    ! ================================================================
    block
        integer, parameter :: n3 = 3
        complex(dp) :: ring(n3, n3), v3(n3, n3), b2(2, 2), vv2(2, 2)
        real(dp) :: d3(n3), lam_ex, err_ring, phi_f, d2v(2), err_rabi
        real(dp), parameter :: t_hop = 1.0_dp
        integer :: i, k

        ! 三格点环: H_ij = -t exp(i phi_f) 相邻; 均匀磁通 phi_f = 2 pi Phi/3
        phi_f = TWOPI * 0.17_dp / 3.0_dp
        ring = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n3
            ring(i, mod(i, n3) + 1) = -t_hop * exp(cmplx(0.0_dp, phi_f, kind=dp))
            ring(mod(i, n3) + 1, i) = conjg(ring(i, mod(i, n3) + 1))
        end do
        call diag_hermitian_matrix(n3, ring, d3, v3, stat)

        ! 闭式谱: lam_n = -2t cos(2 pi (n + Phi)/3), n = 0, 1, 2; 排序后逐位比较
        block
            real(dp) :: closed(3), tmp
            integer :: ii, jj
            do ii = 0, 2
                closed(ii + 1) = -2.0_dp * t_hop * cos(TWOPI * (real(ii, dp) + 0.17_dp) / 3.0_dp)
            end do
            do ii = 1, 2
                do jj = 1, 3 - ii
                    if (closed(jj) > closed(jj + 1)) then
                        tmp = closed(jj); closed(jj) = closed(jj + 1); closed(jj + 1) = tmp
                    end if
                end do
            end do
            err_ring = 0.0_dp
            do k = 1, n3
                err_ring = max(err_ring, abs(d3(k) - closed(k)))
            end do
        end block

        ! 复耦合二能级 (圆偏振驱动): H = [[0, Om/2],[Om*/2, 0]] (|Om| 实数时谱 +-|Om|/2)
        b2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        b2(1, 2) = 0.5_dp * cmplx(0.4_dp, 0.3_dp, kind=dp)
        b2(2, 1) = conjg(b2(1, 2))
        call diag_hermitian_matrix(2, b2, d2v, vv2, stat)
        err_rabi = max(abs(d2v(1) + 0.25_dp), abs(d2v(2) - 0.25_dp))

        n_tests = n_tests + 1
        if (err_ring < 1.0e-13_dp .and. err_rabi < 1.0e-14_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4, A, ES12.4)', " [PASS] Peierls ring closed-form spectrum: err =", &
                  err_ring, "  complex-coupling 2-level: ", err_rabi
        else
            print '(A, ES12.4, A, ES12.4)', " [FAIL] Peierls ring / complex coupling", err_ring, err_rabi
        end if

        ! 回归: 全简并谱 (N = 12, Phi = 0, 所有本征值二重简并) —— 早期
        ! "隔一取一"提取因 LAPACK 简并基序取到平行向量而失效; 现方案
        ! (全候选 + 组内 MGS 选取) 必须保持本征矢完整正交
        block
            integer, parameter :: n12 = 12
            complex(dp) :: ring12(n12, n12), v12(n12, n12), gr(n12, n12)
            real(dp) :: d12(n12), err_r12, err_o12, err_e12, closed12(12), tmp12
            integer :: i12, j12, k12

            ring12 = cmplx(0.0_dp, 0.0_dp, kind=dp)   ! Phi = 0: 实对称环
            do i12 = 1, n12
                j12 = mod(i12, n12) + 1
                ring12(i12, j12) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
                ring12(j12, i12) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
            end do
            call diag_hermitian_matrix(n12, ring12, d12, v12, stat)
            err_r12 = 0.0_dp
            do k12 = 1, n12
                err_r12 = max(err_r12, norm2_complex(matmul(ring12, v12(:, k12)) - &
                                                     d12(k12) * v12(:, k12)))
            end do
            gr = matmul(conjg(transpose(v12)), v12)
            err_o12 = 0.0_dp
            do i12 = 1, n12
                do j12 = 1, n12
                    err_o12 = max(err_o12, abs(gr(i12, j12) - merge(1.0_dp, 0.0_dp, i12 == j12)))
                end do
            end do
            ! 闭式谱匹配 (升序)
            do i12 = 0, n12 - 1
                closed12(i12 + 1) = -2.0_dp * cos(TWOPI * real(i12, dp) / real(n12, dp))
            end do
            do i12 = 1, n12 - 1
                do j12 = 1, n12 - i12
                    if (closed12(j12) > closed12(j12 + 1)) then
                        tmp12 = closed12(j12); closed12(j12) = closed12(j12 + 1)
                        closed12(j12 + 1) = tmp12
                    end if
                end do
            end do
            err_e12 = maxval(abs(d12 - closed12))

            n_tests = n_tests + 1
            if (err_r12 < 1.0e-9_dp .and. err_o12 < 1.0e-9_dp .and. err_e12 < 1.0e-12_dp) then
                n_passed = n_passed + 1
                print '(A, ES10.2, A, ES10.2)', " [PASS] fully-degenerate ring (N=12) regression: res/orth =", &
                      err_r12, " /", err_o12
            else
                print '(A, 3ES12.4)', " [FAIL] degenerate ring regression", err_r12, err_o12, err_e12
            end if
        end block
    end block

    ! ================================================================
    ! 组 5: 简并谱稳健性 (单位阵 + 二重简并) 与相位约定
    ! ================================================================
    block
        integer, parameter :: n5 = 5
        complex(dp) :: id5(n5, n5), v5(n5, n5), degen(n5, n5)
        real(dp) :: d5(n5), err_id, err_dg, err_ph
        integer :: i, j, k

        ! 单位阵: 5 重简并 (本征值全为 1)
        id5 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n5
            id5(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do
        call diag_hermitian_matrix(n5, id5, d5, v5, stat)
        err_id = maxval(abs(d5 - 1.0_dp))

        ! 二重简并: diag(2, 2, -1) 加复非对角微扰仅在简并子空间外为零
        degen = cmplx(0.0_dp, 0.0_dp, kind=dp)
        degen(1, 1) = cmplx(2.0_dp, 0.0_dp, kind=dp)
        degen(2, 2) = cmplx(2.0_dp, 0.0_dp, kind=dp)
        degen(3, 3) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
        degen(4, 4) = cmplx(0.5_dp, 0.0_dp, kind=dp)
        degen(5, 5) = cmplx(-2.5_dp, 0.0_dp, kind=dp)
        degen(3, 4) = cmplx(0.2_dp, 0.1_dp, kind=dp)
        degen(4, 3) = conjg(degen(3, 4))
        call diag_hermitian_matrix(n5, degen, d5, v5, stat)
        ! 升序谱: -2.5, lam_-, lam_+, 2, 2
        ! 其中 lam_pm = -0.25 +- sqrt(0.75^2 + |0.2+0.1i|^2) = -0.25 +- sqrt(0.6125)
        block
            real(dp) :: lam_m, lam_p, disc
            disc = sqrt(0.75_dp**2 + (0.2_dp**2 + 0.1_dp**2))
            lam_m = -0.25_dp - disc
            lam_p = -0.25_dp + disc
            err_dg = 0.0_dp
            err_dg = max(err_dg, abs(d5(1) + 2.5_dp), abs(d5(2) - lam_m))
            err_dg = max(err_dg, abs(d5(3) - lam_p), abs(d5(4) - 2.0_dp), abs(d5(5) - 2.0_dp))
        end block

        ! 相位约定: 每列最大模分量为正实数
        err_ph = 0.0_dp
        do k = 1, n5
            block
                integer :: imax
                real(dp) :: amax
                imax = 1
                amax = abs(v5(1, k))
                do i = 2, n5
                    if (abs(v5(i, k)) > amax) then
                        amax = abs(v5(i, k)); imax = i
                    end if
                end do
                if (amax > 1.0e-12_dp) then
                    err_ph = max(err_ph, abs(aimag(v5(imax, k))))
                    err_ph = max(err_ph, max(0.0_dp, -real(v5(imax, k), dp)))
                end if
            end block
        end do

        n_tests = n_tests + 1
        if (err_id < 1.0e-13_dp .and. err_dg < 1.0e-13_dp .and. err_ph < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES12.4, A, ES12.4)', " [PASS] degeneracies (5-fold / 2-fold) + phase convention: ", &
                  err_id, " /", err_dg
        else
            print '(A, 3ES12.4)', " [FAIL] degeneracy/phase", err_id, err_dg, err_ph
        end if

        ! 简并子空间正交性附加检查 (单位阵的全 5 重简并)
        call diag_hermitian_matrix(n5, id5, d5, v5, stat)
        block
            real(dp) :: err_o2
            complex(dp) :: g(n5, n5)
            g = matmul(conjg(transpose(v5)), v5)
            err_o2 = 0.0_dp
            do i = 1, n5
                do j = 1, n5
                    err_o2 = max(err_o2, abs(g(i, j) - merge(1.0_dp, 0.0_dp, i == j)))
                end do
            end do
            n_tests = n_tests + 1
            if (err_o2 < 1.0e-13_dp) then
                n_passed = n_passed + 1
                print '(A, ES12.4)', " [PASS] orthonormality in a 5-fold degenerate subspace: ", err_o2
            else
                print '(A, ES12.4)', " [FAIL] degenerate subspace orthonormality", err_o2
            end if
        end block
    end block

    print '(A, I2, A, I2, A)', "Hermitian Eigen Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: Hermitian eigen solver test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    pure function norm2_complex(v) result(n)
        complex(dp), intent(in) :: v(:)
        real(dp) :: n
        n = sqrt(real(sum(abs(v)**2), dp))
    end function norm2_complex

    !> 复数矩阵行列式 (Gauss 消元, 部分主元), 返回复数值
    function det_complex(a, n) result(det)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a(n, n)
        complex(dp) :: det
        complex(dp) :: m(n, n), piv, tmp
        integer :: i, k, p
        real(dp) :: amax

        m = a
        det = cmplx(1.0_dp, 0.0_dp, kind=dp)
        do k = 1, n
            p = k
            amax = abs(m(k, k))
            do i = k + 1, n
                if (abs(m(i, k)) > amax) then
                    amax = abs(m(i, k)); p = i
                end if
            end do
            if (amax < 1.0e-300_dp) then
                det = cmplx(0.0_dp, 0.0_dp, kind=dp)
                return
            end if
            if (p /= k) then
                do i = 1, n
                    tmp = m(k, i); m(k, i) = m(p, i); m(p, i) = tmp
                end do
                det = -det
            end if
            det = det * m(k, k)
            do i = k + 1, n
                piv = m(i, k) / m(k, k)
                m(i, k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
                m(i, k + 1:n) = m(i, k + 1:n) - piv * m(k, k + 1:n)
            end do
        end do
    end function det_complex

end program test_hermitian_eig
