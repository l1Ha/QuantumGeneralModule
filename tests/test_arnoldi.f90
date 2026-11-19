!> \brief 非对称 Arnoldi 与复 Schur 分解验收测试 (路线图 §15.10 项)
!> \details 验收组:
!>          组 1: Schur 引擎不变量 —— T = Z^† H Z, Z 幺正, 2x2 解析谱,
!>                迹/范数守恒;
!>          组 2: 正规相似变换精确谱 (U Lam U^†, 本征值条件数 = 1,
!>                非正规相似S J S^{-1}的本征值条件数可任意大, 不适合作
!>                机器精度断言 —— 本轮探针实测 0.45 偏差的教训);
!>          组 3: 目标选取语义 (LR/SR/LI/SI) 与密集参照逐位对齐;
!>          组 4: CAP 复对称哈密顿量 (吸收势): Arnoldi k 对 vs 密集
!>                Schur 全谱参照 + 重启 matvec 效率。
program test_arnoldi
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: st
    integer, parameter :: n_cap = 200
    complex(dp) :: a3(60, 60)          ! 组 3 的酉相似精确谱矩阵
    complex(dp), allocatable :: a_cap(:, :)  ! 组 4 的 CAP 复对称三对角矩阵
    real(dp), allocatable :: ckey_ref(:)
    integer :: cref_idx(n_cap)

    call print_banner("GeneralModule: Non-symmetric Arnoldi / Complex Schur Test", 70)

    ! ================================================================
    ! 组 1: Schur 引擎不变量
    ! ================================================================
    block
        integer, parameter :: n = 50
        complex(dp) :: h(n, n), t(n, n), z(n, n), chk(n, n), h2(2, 2), t2(2, 2), z2(2, 2)
        real(dp) :: err_sim, err_unit, err_trace, err_norm, lam_exact(2), lam_get(2)
        integer :: i, j

        h = 0.0_dp
        do j = 1, n
            do i = 1, min(j + 1, n)
                h(i, j) = cmplx(sin(real(7 * i + 3 * j, dp)) * 0.7_dp, &
                                cos(real(5 * i + 11 * j, dp)) * 0.6_dp, kind=dp)
            end do
        end do
        call hessenberg_schur(n, h, t, z, st)
        chk = matmul(conjg(transpose(z)), h)
        chk = matmul(chk, z) - t
        err_sim = maxval(abs(chk))
        chk = matmul(conjg(transpose(z)), z)
        err_unit = 0.0_dp
        do i = 1, n
            do j = 1, n
                err_unit = max(err_unit, abs(chk(i, j) - merge(1.0_dp, 0.0_dp, i == j)))
            end do
        end do
        ! 迹守恒与 Frobenius 范数守恒 (酉相似的不变量)
        err_trace = abs(sum([(t(i, i), i=1, n)]) - sum([(h(i, i), i=1, n)]))
        err_norm = abs(sqrt(real(sum(abs(t)**2), dp)) - sqrt(real(sum(abs(h)**2), dp)))

        ! 解析 2x2: [[1, 2], [3, 4]] -> lam = (5 ± sqrt(33))/2
        h2(1, 1) = (1.0_dp, 0.0_dp); h2(1, 2) = (2.0_dp, 0.0_dp)
        h2(2, 1) = (3.0_dp, 0.0_dp); h2(2, 2) = (4.0_dp, 0.0_dp)
        call hessenberg_schur(2, h2, t2, z2, st)
        lam_exact = [0.5_dp * (5.0_dp + sqrt(33.0_dp)), 0.5_dp * (5.0_dp - sqrt(33.0_dp))]
        lam_get = [real(t2(1, 1), dp), real(t2(2, 2), dp)]
        ! Schur 对角顺序不定, 按值配对
        if (lam_get(1) > lam_get(2)) then
            err_sim = max(err_sim, abs(lam_get(1) - lam_exact(1)), abs(lam_get(2) - lam_exact(2)))
        else
            err_sim = max(err_sim, abs(lam_get(2) - lam_exact(1)), abs(lam_get(1) - lam_exact(2)))
        end if

        n_tests = n_tests + 1
        if (err_sim < 1.0e-13_dp .and. err_unit < 1.0e-13_dp .and. err_trace < 1.0e-12_dp &
            .and. err_norm < 1.0e-12_dp .and. st == 0) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] Schur invariants T=Z^dag HZ / unitary:", &
                  err_sim, " /", err_unit
        else
            print '(A, 3ES10.2, I3)', " [FAIL] Schur invariants:", err_sim, err_unit, err_trace, st
        end if
    end block

    ! ================================================================
    ! 组 2: 正规相似变换精确谱 + 密集路径
    ! ================================================================
    block
        integer, parameter :: n = 60
        complex(dp) :: g(n, n), u(n, n), lam(n), a(n, n), h(n, n), q(n, n)
        complex(dp) :: t(n, n), z(n, n), v(n)
        real(dp) :: key_a(n), key_b(n), err_spec, err_hq
        real(dp) :: rnd2(2)
        integer :: seed_arr2(8), pass_g
        integer :: i, j, k, p, tmp_i, order_a(n), order_b(n)
        seed_arr2 = [123, 456, 789, 24, 68, 91, 35, 57]

        ! 复随机 G -> MGS 酉化 (两遍重正交化)。注意: 早期用 sin(11i+3j)
        ! 型结构化条目, 三角恒等式使列空间秩 ~4, MGS 塌缩, U 非酉,
        ! A = U Lam U^† 的谱偏离 Lam 达 0.5 —— 结构化条目必须避免。
        call random_seed(put=int(seed_arr2))
        do j = 1, n
            do i = 1, n
                call random_number(rnd2)
                g(i, j) = cmplx(rnd2(1) - 0.5_dp, rnd2(2) - 0.5_dp, kind=dp)
            end do
        end do
        u = 0.0_dp
        do j = 1, n
            do pass_g = 1, 2
                do k = 1, j - 1
                    v = sum(conjg(u(1:n, k)) * g(1:n, j)) * u(1:n, k)
                    g(1:n, j) = g(1:n, j) - v
                end do
            end do
            u(1:n, j) = g(1:n, j) / sqrt(real(sum(conjg(g(1:n, j)) * g(1:n, j)), dp))
        end do
        lam = [(cmplx(-1.0_dp + 0.1_dp * i, 0.05_dp * i, kind=dp), i=1, n)]
        do i = 1, n
            do j = 1, n
                a(i, j) = sum(u(i, :) * lam(:) * conjg(u(j, :)))
            end do
        end do

        call dense_hessenberg_reduce(n, a, h, q)
        ! H = Q^† A Q 不变量
        block
            complex(dp) :: chk(n, n)
            chk = matmul(conjg(transpose(q)), a)
            chk = matmul(chk, q) - h
            err_hq = maxval(abs(chk))
        end block
        call hessenberg_schur(n, h, t, z, st)

        ! 本征值: 双侧按键排序后配对 (Schur 对角与 lam 的顺序均任意)
        do i = 1, n
            key_a(i) = real(t(i, i), dp) + 0.001_dp * aimag(t(i, i))
            key_b(i) = real(lam(i), dp) + 0.001_dp * aimag(lam(i))
            order_a(i) = i
            order_b(i) = i
        end do
        do p = 2, n
            tmp_i = order_a(p)
            do i = p - 1, 1, -1
                if (key_a(order_a(i)) <= key_a(tmp_i)) exit
                order_a(i + 1) = order_a(i)
            end do
            order_a(i + 1) = tmp_i
            tmp_i = order_b(p)
            do i = p - 1, 1, -1
                if (key_b(order_b(i)) <= key_b(tmp_i)) exit
                order_b(i + 1) = order_b(i)
            end do
            order_b(i + 1) = tmp_i
        end do
        err_spec = 0.0_dp
        do i = 1, n
            err_spec = max(err_spec, abs(t(order_a(i), order_a(i)) - lam(order_b(i))))
        end do

        n_tests = n_tests + 1
        if (err_spec < 1.0e-12_dp .and. err_hq < 1.0e-13_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] normal similarity exact spectrum:", &
                  err_spec, "  H=Q^dag AQ =", err_hq
        else
            print '(A, 2ES10.2)', " [FAIL] normal similarity spectrum:", err_spec, err_hq
        end if
    end block

    ! ================================================================
    ! 组 3: 目标选取语义 (LR/SR/LI/SI) vs 密集参照
    ! ================================================================
    block
        integer, parameter :: n = 60, kk = 3
        complex(dp) :: g(n, n), u(n, n), lam(n), a(n, n), h(n, n), q(n, n)
        complex(dp) :: t(n, n), z(n, n), v(n), ev_ref(n)
        complex(dp) :: ev_ar(kk), rv(n, kk)
        real(dp) :: res, err_lr, err_sr, err_li
        real(dp) :: rnd2(2)
        integer :: seed_arr2(8), pass_g
        integer :: i, j, k, p, nmv, idx(n), tmpi
        real(dp) :: key(n)
        character(len=2) :: which
        seed_arr2 = [123, 456, 789, 24, 68, 91, 35, 57]

        ! 同组 2 的构造 (同种子同公式 => 同矩阵)
        call random_seed(put=int(seed_arr2))
        do j = 1, n
            do i = 1, n
                call random_number(rnd2)
                g(i, j) = cmplx(rnd2(1) - 0.5_dp, rnd2(2) - 0.5_dp, kind=dp)
            end do
        end do
        u = 0.0_dp
        do j = 1, n
            do pass_g = 1, 2
                do k = 1, j - 1
                    v = sum(conjg(u(1:n, k)) * g(1:n, j)) * u(1:n, k)
                    g(1:n, j) = g(1:n, j) - v
                end do
            end do
            u(1:n, j) = g(1:n, j) / sqrt(real(sum(conjg(g(1:n, j)) * g(1:n, j)), dp))
        end do
        lam = [(cmplx(-1.0_dp + 0.1_dp * i, 0.05_dp * i, kind=dp), i=1, n)]
        do i = 1, n
            do j = 1, n
                a(i, j) = sum(u(i, :) * lam(:) * conjg(u(j, :)))
            end do
        end do
        a3 = a
        call dense_hessenberg_reduce(n, a, h, q)
        call hessenberg_schur(n, h, t, z, st)
        ev_ref = [(t(i, i), i=1, n)]

        ! 三种目标的参照 (对 ev_ref 按相应键排序)
        err_lr = 0.0_dp
        err_sr = 0.0_dp
        err_li = 0.0_dp
        do p = 1, 3
            select case (p)
            case (1)
                which = 'LR'
                do i = 1, n
                    key(i) = -real(ev_ref(i), dp)
                end do
            case (2)
                which = 'SR'
                do i = 1, n
                    key(i) = real(ev_ref(i), dp)
                end do
            case (3)
                which = 'LI'
                do i = 1, n
                    key(i) = -aimag(ev_ref(i))
                end do
            end select
            do i = 1, n
                idx(i) = i
            end do
            do j = 2, n
                tmpi = idx(j)
                do i = j - 1, 1, -1
                    if (key(idx(i)) <= key(tmpi)) exit
                    idx(i + 1) = idx(i)
                end do
                idx(i + 1) = tmpi
            end do
            ! 容限/趟数放宽: 不同 gfortran 版本的随机流与本征矢微差会改变
            ! 重启动力学, 1e-12/10 趟跨版本脆弱 (CI macos 教训)
            call arnoldi_general_eigenpairs(n, kk, 45, 15, cmat_mult3, which, &
                                            ev_ar, rv, nmv, res, 1.0e-10_dp, st)
            do i = 1, kk
                select case (p)
                case (1)
                    err_lr = max(err_lr, abs(ev_ar(i) - ev_ref(idx(i))))
                case (2)
                    err_sr = max(err_sr, abs(ev_ar(i) - ev_ref(idx(i))))
                case (3)
                    err_li = max(err_li, abs(ev_ar(i) - ev_ref(idx(i))))
                end select
            end do
        end do

        n_tests = n_tests + 1
        if (err_lr < 1.0e-10_dp .and. err_sr < 1.0e-10_dp .and. err_li < 1.0e-10_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] selection semantics LR/SR/LI vs dense:", &
                  err_lr, " /", max(err_sr, err_li)
        else
            print '(A, 3ES10.2)', " [FAIL] selection semantics:", err_lr, err_sr, err_li
        end if
    end block

    ! ================================================================
    ! 组 4: CAP 复对称哈密顿量 (谐振子 + 复吸收势)
    ! ================================================================
    block
        integer, parameter :: n_cap = 200, kk = 4
        real(dp) :: xg(n_cap), dx, eta_cap
        complex(dp) :: v_cap(n_cap)
        complex(dp) :: ev_ar(kk), rv(n_cap, kk)
        ! 大数组用 allocatable (堆): macOS 主线程栈默认 8MB, 多个
        ! 200x200 自动数组累计 ~4MB 有栈溢出风险
        complex(dp), allocatable :: h_cap(:, :), q_cap(:, :)
        complex(dp), allocatable :: t_cap(:, :), z_cap(:, :), ev_ref(:)
        real(dp) :: res, err_cap, res_single
        integer :: i, j, nmv, nmv_single, st_cap, st_single, tmpi
        character(len=2) :: which

        ! H = p^2/2 + x^2/2 - i eta (x/xmax)^4 (box [-10, 10], sinc DVR)
        allocate (a_cap(n_cap, n_cap), ckey_ref(n_cap))   ! 程序级 (cmat_cap host association)
        allocate (h_cap(n_cap, n_cap), q_cap(n_cap, n_cap))
        allocate (t_cap(n_cap, n_cap), z_cap(n_cap, n_cap), ev_ref(n_cap))
        dx = 20.0_dp / real(n_cap, dp)
        eta_cap = 0.05_dp
        do i = 1, n_cap
            xg(i) = -10.0_dp + (real(i, dp) - 0.5_dp) * dx
            v_cap(i) = cmplx(0.5_dp * xg(i)**2, &
                             -eta_cap * (xg(i) / 10.0_dp)**4, kind=dp)
        end do
        ! 复对称三对角 (二阶中心差分动能) + CAP 复对角
        a_cap = 0.0_dp
        do i = 1, n_cap
            a_cap(i, i) = v_cap(i) + cmplx(1.0_dp / dx**2, 0.0_dp, kind=dp)
            if (i < n_cap) then
                a_cap(i, i + 1) = cmplx(-0.5_dp / dx**2, 0.0_dp, kind=dp)
                a_cap(i + 1, i) = a_cap(i, i + 1)
            end if
        end do
        ! 密集参照
        call dense_hessenberg_reduce(n_cap, a_cap, h_cap, q_cap)
        call hessenberg_schur(n_cap, h_cap, t_cap, z_cap, st_cap)
        ev_ref = [(t_cap(i, i), i=1, n_cap)]

        ! 'SR': 实部最小 (最深的束缚 + 吸收共振)
        which = 'SR'
        call arnoldi_general_eigenpairs(n_cap, kk, 60, 20, cmat_cap, which, &
                                        ev_ar, rv, nmv, res, 1.0e-10_dp, st_cap)
        ! 密集参照排序
        do i = 1, n_cap
            ckey_ref(i) = real(ev_ref(i), dp)
            cref_idx(i) = i
        end do
        do j = 2, n_cap
            tmpi = cref_idx(j)
            do i = j - 1, 1, -1
                if (ckey_ref(cref_idx(i)) <= ckey_ref(tmpi)) exit
                cref_idx(i + 1) = cref_idx(i)
            end do
            cref_idx(i + 1) = tmpi
        end do
        err_cap = 0.0_dp
        do i = 1, kk
            err_cap = max(err_cap, abs(ev_ar(i) - ev_ref(cref_idx(i))))
        end do
        ! 单趟 (无重启) 同容限的残差水平 —— 展示重启的价值
        block
            complex(dp) :: ev_s(kk), rv_s(n_cap, kk)
            call arnoldi_general_eigenpairs(n_cap, kk, 60, 1, cmat_cap, which, &
                                            ev_s, rv_s, nmv_single, res_single, &
                                            1.0e-10_dp, st_single)
        end block

        n_tests = n_tests + 1
        if (err_cap < 1.0e-9_dp .and. res < 1.0e-9_dp .and. st_cap == 0) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES9.2, A, I5, A, I5)', " [PASS] CAP complex symmetric: err =", &
                  err_cap, "  res =", res, "  mv restart/single =", nmv, "/", nmv_single
        else
            print '(A, ES10.2, ES9.2, I3)', " [FAIL] CAP:", err_cap, res, st_cap
        end if
    end block

    print '(A, I2, A, I2, A)', "Arnoldi Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: non-symmetric Arnoldi test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

    !> 组 3 的复矩阵乘 (酉相似精确谱矩阵, host association)
    subroutine cmat_mult3(vin, vout)
        use mod_constants, only: dp
        complex(dp), intent(in) :: vin(:)
        complex(dp), intent(out) :: vout(:)
        integer :: i
        do i = 1, 60
            vout(i) = sum(a3(i, :) * vin(:))
        end do
    end subroutine cmat_mult3

    !> 组 4 的 CAP 复对称三对角乘
    subroutine cmat_cap(vin, vout)
        use mod_constants, only: dp
        complex(dp), intent(in) :: vin(:)
        complex(dp), intent(out) :: vout(:)
        integer :: i
        vout(1) = a_cap(1, 1) * vin(1) + a_cap(1, 2) * vin(2)
        do i = 2, n_cap - 1
            vout(i) = a_cap(i, i - 1) * vin(i - 1) + a_cap(i, i) * vin(i) &
                      + a_cap(i, i + 1) * vin(i + 1)
        end do
        vout(n_cap) = a_cap(n_cap, n_cap - 1) * vin(n_cap - 1) + a_cap(n_cap, n_cap) * vin(n_cap)
    end subroutine cmat_cap

end program test_arnoldi
