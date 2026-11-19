!> \brief 非对称 Arnoldi 与复 Schur 分解 (路线图 §15.10「非对称 Arnoldi」项)
!> \details 一般复矩阵的 Krylov 本征求解:
!>          1) arnoldi_factorize: 修正 Gram-Schmidt + 二次正交化的 m 步
!>             Arnoldi 分解 A V_m = V_m H_m + h_{m+1,m} v_{m+1} e_m^T;
!>          2) hessenberg_schur: 移位复 QR (Wilkinson 移位 + 异常移位)
!>             将 Hessenberg 矩阵三角化为 Schur 形式 T = Z^† H Z, 旋转
!>             累积给出 Schur 向量 —— 库内此前仅有实对称 (tql2) 与复
!>             Hermitian (加倍归约) 求解器, 一般复矩阵 (如含吸收势
!>             H - iW) 的谱问题无解;
!>          3) hessenberg_eigenvectors: 三角回代 + 近简并组内 MGS;
!>          4) arnoldi_general_eigenpairs: 重启式 k 对选取 (ARPACK
!>             记法 LM/SM/LR/SR/LI/SI)。残差后验估计
!>             r_i = |h_{m+1,m}| |(Z y_tri)_m|: 回代矢量 y_tri 是上三角
!>             本征矢, 对 i < m 恒有 y_tri(m) = 0, 直接用 y_tri(m) 是
!>             常见错误; 必须经 Schur 向量 Z 变换后取末分量。
!>          5) dense_hessenberg_reduce: 复 Householder 稠密约简, 供
!>             全谱参照路径 (验收测试)。
!> \author LiHao
module mod_arnoldi
    use mod_constants, only: dp
    implicit none
    private

    public :: arnoldi_factorize
    public :: hessenberg_schur
    public :: hessenberg_eigenvectors
    public :: dense_hessenberg_reduce
    public :: arnoldi_general_eigenpairs

    real(dp), parameter :: EPS_MACH = 2.220446049250313e-16_dp

contains

    !> \brief 复 Givens 旋转: G = [[c, s], [-conj(s), c]], c 实数,
    !>        使 G [a; b] = [r; 0]。推导: -conj(s) a + c b = 0
    !>        => s = (c/a) conj(b), c = |a|/sqrt(|a|^2+|b|^2)。
    pure subroutine complex_givens(a, b, c, s)
        complex(dp), intent(in) :: a, b
        real(dp), intent(out) :: c
        complex(dp), intent(out) :: s
        real(dp) :: r

        if (abs(b) == 0.0_dp) then
            c = 1.0_dp
            s = cmplx(0.0_dp, 0.0_dp, kind=dp)
        else if (abs(a) == 0.0_dp) then
            c = 0.0_dp
            s = cmplx(1.0_dp, 0.0_dp, kind=dp)
        else
            r = sqrt(abs(a)**2 + abs(b)**2)
            c = abs(a) / r
            ! -conj(s)*a + c*b = 0 => s = c*conj(b)/conj(a)
            ! (早期版本漏掉 conj(a): -conj(s)*a + c*b = c*b*(1 - a/conj(a))
            !  对复 a 不为零, Schur 向量失去幺正性 —— 不变量检查捕获)
            s = (c / conjg(a)) * conjg(b)
        end if
    end subroutine complex_givens

    !> \brief m 步 Arnoldi 分解 (复 MGS + 二次正交化)
    !> \param[out] v_basis (n, m_eff) 正交基 (m_eff <= m_max)
    !> \param[out] h_hess 上 Hessenberg 矩阵 (1..m_eff 有效)
    !> \param[out] h_sub h_{m_eff+1, m_eff} 残差幅度 (breakdown 时为 0)
    !> \param[out] n_mv 实际 matvec 次数; stat: 0 成功, 2 起始向量零
    subroutine arnoldi_factorize(n, m_max, c_mult, v_start, v_basis, &
                                 h_hess, h_sub, m_eff, n_mv, stat)
        integer, intent(in) :: n, m_max
        interface
            subroutine c_mult(v_in, v_out)
                use mod_constants, only: dp
                complex(dp), intent(in) :: v_in(:)
                complex(dp), intent(out) :: v_out(:)
            end subroutine c_mult
        end interface
        complex(dp), intent(in) :: v_start(:)
        complex(dp), intent(out) :: v_basis(:, :), h_hess(:, :)
        complex(dp), intent(out) :: h_sub
        integer, intent(out) :: m_eff, n_mv, stat

        complex(dp) :: w(n), proj
        real(dp) :: tol_break, nrm0, reorth_max
        integer :: j, i, pass2

        stat = 0
        n_mv = 0
        m_eff = 0
        h_hess = cmplx(0.0_dp, 0.0_dp, kind=dp)
        tol_break = 1.0e-13_dp * sqrt(real(max(n, 1), dp))

        nrm0 = sqrt(real(sum(conjg(v_start) * v_start), dp))
        if (nrm0 <= 1.0e-300_dp) then
            stat = 2
            h_sub = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if
        v_basis(:, 1) = v_start / nrm0

        do j = 1, m_max
            call c_mult(v_basis(:, j), w)
            n_mv = n_mv + 1
            ! 修正 Gram-Schmidt; 第一遍最大投影超阈值时做第二遍 (二次正交化)
            do pass2 = 1, 2
                reorth_max = 0.0_dp
                do i = 1, j
                    proj = sum(conjg(v_basis(1:n, i)) * w(1:n))
                    reorth_max = max(reorth_max, abs(proj))
                    h_hess(i, j) = h_hess(i, j) + proj
                    w(1:n) = w(1:n) - proj * v_basis(1:n, i)
                end do
                if (pass2 == 1 .and. reorth_max > 1.0e-9_dp) cycle
                exit
            end do
            h_sub = cmplx(sqrt(real(sum(conjg(w) * w), dp)), 0.0_dp, kind=dp)
            if (abs(h_sub) <= tol_break) then
                ! 不变子空间: Krylov 空间饱和, 合法提前终止
                h_sub = cmplx(0.0_dp, 0.0_dp, kind=dp)
                m_eff = j
                return
            end if
            if (j == m_max) then
                m_eff = m_max
                return
            end if
            v_basis(:, j + 1) = w / abs(h_sub)
            h_hess(j + 1, j) = h_sub
        end do
    end subroutine arnoldi_factorize

    !> \brief 移位复 QR: Hessenberg -> 上三角 Schur 形式 T = Z^† H Z
    !> \details Wilkinson 移位 (活动窗尾部 2x2 本征值中靠近 H(ihi,ihi)
    !>          者); deflation 判据 |H(l,l-1)| <= eps(|H(l-1,l-1)|+|H(l,l)|)
    !>          加绝对地板; 每 10 次停滞注入异常移位 mu = H(ihi,ihi)
    !>          + 0.75|H(ihi,ihi-1)| 破除对称滞环 (EISPACK COMQR 惯例)。
    !>          行旋转须先保存旧行 (原地覆盖会毁掉旋转所需的旧值);
    !>          T/Z 仅 1..m 子块有意义, 其余置零 (intent(out) 契约)。
    subroutine hessenberg_schur(m, h_in, t, z, stat)
        integer, intent(in) :: m
        complex(dp), intent(in) :: h_in(:, :)
        complex(dp), intent(out) :: t(:, :), z(:, :)
        integer, intent(out) :: stat

        complex(dp) :: sn(m), tmp_k(m), tmp_k1(m), lam1, lam2, disc, mu
        real(dp) :: cr_arr(m)
        real(dp) :: cr, defl_floor
        integer :: ihi, l, k, iters_this, total_iters, i, nrow

        stat = 0
        t = cmplx(0.0_dp, 0.0_dp, kind=dp)
        z = cmplx(0.0_dp, 0.0_dp, kind=dp)
        t(1:m, 1:m) = h_in(1:m, 1:m)
        do i = 1, m
            z(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do
        defl_floor = max(sqrt(real(sum(abs(t(1:m, 1:m))**2), dp)), 1.0e-300_dp) * EPS_MACH
        ihi = m
        iters_this = 0
        total_iters = 0

        do while (ihi > 1)
            ! 找活动窗起点 l (副对角可忽略处 deflation)
            l = ihi
            do while (l > 1)
                if (abs(t(l, l - 1)) <= EPS_MACH * (abs(t(l - 1, l - 1)) + abs(t(l, l))) &
                    .or. abs(t(l, l - 1)) <= defl_floor) then
                    t(l, l - 1) = cmplx(0.0_dp, 0.0_dp, kind=dp)
                    exit
                end if
                l = l - 1
            end do
            if (l == ihi) then
                ihi = ihi - 1
                iters_this = 0
                cycle
            end if

            iters_this = iters_this + 1
            total_iters = total_iters + 1
            if (total_iters > 100 * m) then
                stat = 1
                return
            end if
            if (mod(iters_this, 10) == 0) then
                mu = t(ihi, ihi) + cmplx(0.75_dp, 0.0_dp, kind=dp) * abs(t(ihi, ihi - 1))
            else
                disc = sqrt((t(ihi - 1, ihi - 1) - t(ihi, ihi))**2 &
                            + 4.0_dp * t(ihi - 1, ihi) * t(ihi, ihi - 1))
                lam1 = 0.5_dp * (t(ihi - 1, ihi - 1) + t(ihi, ihi)) + 0.5_dp * disc
                lam2 = 0.5_dp * (t(ihi - 1, ihi - 1) + t(ihi, ihi)) - 0.5_dp * disc
                if (abs(lam1 - t(ihi, ihi)) <= abs(lam2 - t(ihi, ihi))) then
                    mu = lam1
                else
                    mu = lam2
                end if
            end if

            ! 显式移位 QR 步: 减移位; 行旋转消副对角 (R); 列旋转还原
            ! Hessenberg (R Q); 旋转同时累计 Z <- Z Q, 终态 T = Z^† H_orig Z
            do i = l, ihi
                t(i, i) = t(i, i) - mu
            end do
            do k = l, ihi - 1
                call complex_givens(t(k, k), t(k + 1, k), cr, sn(k))
                cr_arr(k) = cr
                ! 行旋转必须覆盖 k..m 全部非零列: 活动窗的行在已收敛列
                ! ihi+1..m 上有耦合元, 只转 k..ihi 会留下与 Z 不一致的
                ! 旧值, 破坏 T = Z^† H Z (不变量检查捕获)
                tmp_k(1:m - k + 1) = t(k, k:m)
                tmp_k1(1:m - k + 1) = t(k + 1, k:m)
                t(k, k:m) = cr * tmp_k(1:m - k + 1) + sn(k) * tmp_k1(1:m - k + 1)
                t(k + 1, k:m) = -conjg(sn(k)) * tmp_k(1:m - k + 1) &
                                + cr * tmp_k1(1:m - k + 1)
                t(k + 1, k) = cmplx(0.0_dp, 0.0_dp, kind=dp)   ! 结构性精确为零
            end do
            do k = l, ihi - 1
                ! c_k 必须按 k 取 (cr_arr(k)): 标量 cr 在此已只剩行循环
                ! 末次的残留值, 曾致全部列旋转用错 c, Schur 向量失去幺正
                nrow = min(k + 1, ihi)
                tmp_k(1:nrow) = t(1:nrow, k)
                tmp_k1(1:nrow) = t(1:nrow, k + 1)
                t(1:nrow, k) = cr_arr(k) * tmp_k(1:nrow) + conjg(sn(k)) * tmp_k1(1:nrow)
                t(1:nrow, k + 1) = -sn(k) * tmp_k(1:nrow) + cr_arr(k) * tmp_k1(1:nrow)
                ! Schur 向量累计 (全行)
                tmp_k(1:m) = z(1:m, k)
                tmp_k1(1:m) = z(1:m, k + 1)
                z(1:m, k) = cr_arr(k) * tmp_k(1:m) + conjg(sn(k)) * tmp_k1(1:m)
                z(1:m, k + 1) = -sn(k) * tmp_k(1:m) + cr_arr(k) * tmp_k1(1:m)
            end do
            do i = l, ihi
                t(i, i) = t(i, i) + mu
            end do
        end do
    end subroutine hessenberg_schur

    !> \brief Schur 形式的本征矢: 三角回代 + 近简并组内 MGS 正交化
    !> \details 第 i 个本征值 lam_i = T(i,i) 的本征矢 y 支持 1..i
    !>          (相异谱时), 回代 y(j) = -(sum_k T(j,k) y(k)) / (T(j,j)-lam_i);
    !>          T(j,j) 与 lam_i 近重合时以 eps*||T|| 地板防零除, 近简并
    !>          组内对已得向量 MGS, 塌缩时确定性扰动重试一次。
    subroutine hessenberg_eigenvectors(m, t, z, vecs, stat)
        integer, intent(in) :: m
        complex(dp), intent(in) :: t(:, :), z(:, :)
        complex(dp), intent(out) :: vecs(:, :)
        integer, intent(out) :: stat

        complex(dp) :: y(m), x(m), den, ssum, noise
        real(dp) :: tnorm, clus_tol, floor_den, nrm, proj
        integer :: i, j, k, jj, n_retry

        stat = 0
        vecs = cmplx(0.0_dp, 0.0_dp, kind=dp)
        tnorm = max(sqrt(real(sum(abs(t(1:m, 1:m))**2), dp)), 1.0e-300_dp)
        clus_tol = sqrt(EPS_MACH) * tnorm
        floor_den = EPS_MACH * tnorm

        do i = 1, m
            n_retry = 0
            do
                y = cmplx(0.0_dp, 0.0_dp, kind=dp)
                y(i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
                do j = i - 1, 1, -1
                    ssum = cmplx(0.0_dp, 0.0_dp, kind=dp)
                    do k = j + 1, i
                        ssum = ssum + t(j, k) * y(k)
                    end do
                    den = t(j, j) - t(i, i)
                    if (abs(den) < floor_den) then
                        den = cmplx(floor_den, 0.0_dp, kind=dp)
                        if (n_retry > 0) y(j) = y(j) + noise * cmplx(0.5_dp, 0.5_dp, kind=dp)
                    end if
                    y(j) = -ssum / den
                end do
                x = matmul(z(1:m, 1:m), y)
                nrm = sqrt(real(sum(conjg(x) * x), dp))
                if (nrm <= 1.0e-300_dp) then
                    stat = 1
                    return
                end if
                x = x / nrm
                ! 近简并组内正交化 (与已求出的近重合本征矢)
                do jj = 1, i - 1
                    if (abs(t(jj, jj) - t(i, i)) > clus_tol) cycle
                    proj = sum(conjg(vecs(1:m, jj)) * x(1:m))
                    x(1:m) = x(1:m) - proj * vecs(1:m, jj)
                end do
                nrm = sqrt(real(sum(conjg(x) * x), dp))
                if (nrm > 0.05_dp .or. n_retry >= 1) then
                    if (nrm > 1.0e-300_dp) x = x / nrm
                    exit
                end if
                ! 塌缩: 确定性扰动后重试一次
                n_retry = n_retry + 1
                noise = cmplx(sin(real(i, dp) * 7.13_dp + 1.7_dp), &
                              cos(real(i, dp) * 3.71_dp + 0.9_dp), kind=dp)
            end do
            vecs(1:m, i) = x
        end do
    end subroutine hessenberg_eigenvectors

    !> \brief 稠密复矩阵 -> 上 Hessenberg (复 Householder), H = Q^† A Q
    subroutine dense_hessenberg_reduce(n, a_in, h, q)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a_in(:, :)
        complex(dp), intent(out) :: h(:, :), q(:, :)
        integer :: k, i
        complex(dp) :: a(n, n), u(n), alpha, wrow
        real(dp) :: xn, un2, two_over

        a = a_in(1:n, 1:n)
        q = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n
            q(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do
        do k = 1, n - 2
            xn = sqrt(real(sum(conjg(a(k + 1:n, k)) * a(k + 1:n, k)), dp))
            if (xn <= 1.0e-300_dp) cycle
            u = cmplx(0.0_dp, 0.0_dp, kind=dp)
            u(1:n - k) = a(k + 1:n, k)
            if (abs(u(1)) <= 1.0e-300_dp) then
                alpha = cmplx(-xn, 0.0_dp, kind=dp)
            else
                alpha = -xn * u(1) / abs(u(1))    ! -exp(i arg x1) * |x|
            end if
            u(1) = u(1) - alpha
            un2 = real(sum(conjg(u) * u), dp)
            if (un2 <= 1.0e-300_dp) cycle
            two_over = 2.0_dp / un2
            ! 左乘: 行 k+1..n, 列 k..n
            do i = k, n
                wrow = sum(conjg(u(1:n - k)) * a(k + 1:n, i))
                a(k + 1:n, i) = a(k + 1:n, i) - two_over * wrow * u(1:n - k)
            end do
            ! 右乘: 行 1..n, 列 k+1..n
            do i = 1, n
                wrow = sum(a(i, k + 1:n) * u(1:n - k))
                a(i, k + 1:n) = a(i, k + 1:n) - two_over * wrow * conjg(u(1:n - k))
            end do
            ! Q <- Q P
            do i = 1, n
                wrow = sum(q(i, k + 1:n) * u(1:n - k))
                q(i, k + 1:n) = q(i, k + 1:n) - two_over * wrow * conjg(u(1:n - k))
            end do
            a(k + 2:n, k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end do
        h = a
    end subroutine dense_hessenberg_reduce

    !> \brief 重启式 Arnoldi 求 k 个目标本征对 (ARPACK 记法)
    !> \details 设计: 锁定 + 重启正交化 (保持算子 A 不动, 残差恒相对 A):
    !>          每趟 Arnoldi 后, 残差 < tol 的 Ritz 对进入锁定表 (按方向
    !>          平行去重); 重启向量 = 未收敛最优 Ritz 矢量 + 几何衰减噪声,
    !>          并对锁定矢量正交化 (隐式紧缩, 使后续趟探索新方向)。
    !>          重启空间的 Krylov 会对同一 lambda 产生渐近重复 Ritz 对,
    !>          方向平行判据 |<x_i, x_kept>| > 0.999 尺度无关地滤除
    !>          (残差带判据 |dlam| <= 10*max(res) 在单趟大残差阶段会误吞
    !>          gap~1 的真谱, 已弃用); 平行候选若残差更小则原位升级。
    !> \param[in] which 'LM'/'SM' (模最大/最小), 'LR'/'SR' (实部最大/
    !>              最小), 'LI'/'SI' (虚部最大/最小; CAP 应用常用 'SR')
    !> \param[out] evals(k) rvecs(n,k) 锁定优先、键序排列的 k 个本征对
    !> \param[out] res_max 最坏 Ritz 对残差 |h_{m+1,m}| |(Z y)_m|
    !> \param[out] n_matvec 累计 matvec
    !> \param[out] stat 0 全部收敛 (res < tol), 1 尽力而为, 2/3/4 内部失败
    subroutine arnoldi_general_eigenpairs(n_dim, k, m_max, n_restart, c_mult, &
                                          which, evals, rvecs, n_matvec, res_max, &
                                          tol, stat)
        integer, intent(in) :: n_dim, k, m_max, n_restart
        interface
            subroutine c_mult(v_in, v_out)
                use mod_constants, only: dp
                complex(dp), intent(in) :: v_in(:)
                complex(dp), intent(out) :: v_out(:)
            end subroutine c_mult
        end interface
        character(len=*), intent(in) :: which
        complex(dp), intent(out) :: evals(k), rvecs(:, :)
        integer, intent(out) :: n_matvec
        real(dp), intent(out) :: res_max
        real(dp), intent(in), optional :: tol
        integer, intent(out) :: stat

        complex(dp), allocatable :: v_basis(:, :), h_hess(:, :), t(:, :), z(:, :)
        complex(dp), allocatable :: svecs(:, :), v_start(:)
        complex(dp) :: h_sub, noise_v, proj
        real(dp) :: proj_dup
        real(dp) :: rnd(2)
        complex(dp) :: cand_lam(m_max), cand_vec(n_dim)
        real(dp) :: cand_key(m_max)
        logical :: is_dup
        complex(dp) :: lam_tmp, vec_tmp(n_dim)
        real(dp) :: tmp_r
        complex(dp) :: lock_lam(k), lock_vec(n_dim, k)
        complex(dp) :: open_lam(k), open_vec(n_dim, k)
        real(dp) :: cand_res, lock_res(k), open_res(k)
        real(dp) :: tol_v, noise_amp, nrm, nrm_o
        integer :: order(m_max)
        integer :: seed_arr(8)
        integer :: m_eff, st, n_mv, pass, i, j, q, p
        integer :: n_lock, n_open, cand_n

        stat = 1
        tol_v = 1.0e-10_dp
        if (present(tol)) tol_v = tol
        n_matvec = 0
        res_max = 1.0e30_dp
        seed_arr = [42, 17, 93, 5, 61, 28, 74, 9]
        call random_seed(put=seed_arr)

        allocate (v_basis(n_dim, m_max), h_hess(m_max, m_max))
        allocate (t(m_max, m_max), z(m_max, m_max), svecs(m_max, m_max))
        allocate (v_start(n_dim))

        ! 确定性起始向量 (可复现, 与厚重启 Lanczos 轮同惯例)
        do j = 1, n_dim
            call random_number(rnd)
            v_start(j) = cmplx(rnd(1) - 0.5_dp, rnd(2) - 0.5_dp, kind=dp) &
                         + cmplx(0.13_dp * sin(real(j, dp)), 0.11_dp * cos(real(j, dp)), kind=dp)
        end do
        v_start(1:n_dim) = v_start(1:n_dim) &
                           / sqrt(real(sum(conjg(v_start(1:n_dim)) * v_start(1:n_dim)), dp))
        n_lock = 0
        n_open = 0
        open_res = 1.0e30_dp

        do pass = 1, max(n_restart, 1)
            call arnoldi_factorize(n_dim, m_max, c_mult, v_start, v_basis, &
                                   h_hess, h_sub, m_eff, n_mv, st)
            n_matvec = n_matvec + n_mv
            if (st /= 0 .or. m_eff < 1) then
                stat = 2
                return
            end if
            call hessenberg_schur(m_eff, h_hess, t, z, st)
            if (st /= 0) then
                stat = 3
                return
            end if
            call hessenberg_eigenvectors(m_eff, t, z, svecs, st)
            if (st /= 0) then
                stat = 4
                return
            end if

            ! 候选按键排序
            do i = 1, m_eff
                cand_key(i) = target_key(t(i, i))
            end do
            call arg_sort(m_eff, cand_key, order)

            ! 逐候选处理: 与锁定表平行者跳过 (平行且残差更小则升级锁定对);
            ! 残差 < tol 收敛者锁定; 其余并入未收敛最优集 (键序贪心)
            do q = 1, m_eff
                i = order(q)
                if (q > k) cycle   ! 仅键序前 k 名是目标 (否则谱顶的快收敛对
                                   !  会以 'SR' 等目标误入锁定表)
                cand_lam(i) = t(i, i)
                cand_res = abs(h_sub) * abs(svecs(m_eff, i))
                cand_vec(1:n_dim) = matmul(v_basis(1:n_dim, 1:m_eff), svecs(1:m_eff, i))
                ! 与锁定表比对 (方向平行)
                is_dup = .false.
                do p = 1, n_lock
                    proj_dup = abs(sum(conjg(lock_vec(1:n_dim, p)) * cand_vec(1:n_dim)))
                    if (proj_dup > 0.999_dp) then
                        if (cand_res < lock_res(p)) then
                            lock_lam(p) = cand_lam(i)
                            lock_vec(1:n_dim, p) = cand_vec(1:n_dim)
                            lock_res(p) = cand_res
                        end if
                        is_dup = .true.
                        exit
                    end if
                end do
                if (is_dup) cycle
                if (cand_res < tol_v .and. n_lock < k) then
                    n_lock = n_lock + 1
                    lock_lam(n_lock) = cand_lam(i)
                    lock_vec(1:n_dim, n_lock) = cand_vec(1:n_dim)
                    lock_res(n_lock) = cand_res
                    cycle
                end if
                ! 未收敛: 并入 open 集 (按键贪心, 平行升级)
                call merge_open(cand_lam(i), cand_vec, cand_res)
            end do

            ! 锁定发生后, 从 open 集清除与其平行的陈旧条目 (该本征对已
            ! 收敛, open 集里的旧残差只会在最终 res_max 里自相矛盾)
            q = 1
            do while (q <= n_open)
                is_dup = .false.
                do p = 1, n_lock
                    proj_dup = abs(sum(conjg(lock_vec(1:n_dim, p)) * open_vec(1:n_dim, q)))
                    if (proj_dup > 0.999_dp) then
                        is_dup = .true.
                        exit
                    end if
                end do
                if (is_dup) then
                    do p = q, n_open - 1
                        open_lam(p) = open_lam(p + 1)
                        open_vec(1:n_dim, p) = open_vec(1:n_dim, p + 1)
                        open_res(p) = open_res(p + 1)
                    end do
                    n_open = n_open - 1
                else
                    q = q + 1
                end if
            end do

            if (n_lock >= k) then
                stat = 0
                exit
            end if

            ! 重启向量 = 未收敛最优 + 几何衰减噪声, 对锁定矢量正交化
            v_start = cmplx(0.0_dp, 0.0_dp, kind=dp)
            do i = 1, n_open
                if (sqrt(real(sum(conjg(open_vec(1:n_dim, i)) * open_vec(1:n_dim, i)), dp)) &
                    > 1.0e-300_dp) then
                    v_start(1:n_dim) = v_start(1:n_dim) + open_vec(1:n_dim, i)
                end if
            end do
            noise_amp = 1.0e-2_dp * (0.1_dp**real(pass - 1, dp))
            do j = 1, n_dim
                call random_number(rnd)
                v_start(j) = v_start(j) &
                             + noise_amp * cmplx(rnd(1) - 0.5_dp, rnd(2) - 0.5_dp, kind=dp)
            end do
            do p = 1, n_lock
                proj = sum(conjg(lock_vec(1:n_dim, p)) * v_start(1:n_dim))
                v_start(1:n_dim) = v_start(1:n_dim) - proj * lock_vec(1:n_dim, p)
            end do
            nrm = sqrt(real(sum(conjg(v_start(1:n_dim)) * v_start(1:n_dim)), dp))
            if (nrm <= 1.0e-300_dp) then
                do j = 1, n_dim
                    call random_number(rnd)
                    v_start(j) = cmplx(rnd(1) - 0.5_dp, rnd(2) - 0.5_dp, kind=dp)
                end do
                do p = 1, n_lock
                    proj = sum(conjg(lock_vec(1:n_dim, p)) * v_start(1:n_dim))
                    v_start(1:n_dim) = v_start(1:n_dim) - proj * lock_vec(1:n_dim, p)
                end do
                nrm = sqrt(real(sum(conjg(v_start(1:n_dim)) * v_start(1:n_dim)), dp))
            end if
            v_start(1:n_dim) = v_start(1:n_dim) / nrm
        end do

        ! 输出: 锁定对优先 (键序), 未收敛最优补足
        cand_n = 0
        do p = 1, n_lock
            cand_n = cand_n + 1
            evals(cand_n) = lock_lam(p)
            rvecs(1:n_dim, cand_n) = lock_vec(1:n_dim, p)
        end do
        do p = 1, n_open
            if (cand_n >= k) exit
            cand_n = cand_n + 1
            evals(cand_n) = open_lam(p)
            rvecs(1:n_dim, cand_n) = open_vec(1:n_dim, p)
        end do
        do p = cand_n + 1, k
            evals(p) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            rvecs(1:n_dim, p) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end do
        ! res_max 只统计实际返回的对 (锁定 + 填充用的 open 前段):
        ! 未被返回的 open 条目不参与, 否则 stat=0 与 res_max > tol 自相矛盾
        res_max = 0.0_dp
        do p = 1, n_lock
            res_max = max(res_max, lock_res(p))
        end do
        do p = 1, min(n_open, max(k - n_lock, 0))
            res_max = max(res_max, open_res(p))
        end do
        if (res_max <= tol_v .and. cand_n >= k) stat = 0

    contains

        !> 目标键 (ARPACK 记法): 键小者排前
        pure function target_key(lam_v) result(kv)
            complex(dp), intent(in) :: lam_v
            real(dp) :: kv
            select case (which)
            case ('LM')
                kv = -abs(lam_v)
            case ('SM')
                kv = abs(lam_v)
            case ('SR')
                kv = real(lam_v, dp)
            case ('LI')
                kv = -aimag(lam_v)
            case ('SI')
                kv = aimag(lam_v)
            case default   ! 'LR'
                kv = -real(lam_v, dp)
            end select
        end function target_key

        !> 插入排序产生键升序的索引排列 (不搬数据)
        pure subroutine arg_sort(cnt, keys, order_out)
            integer, intent(in) :: cnt
            real(dp), intent(in) :: keys(:)
            integer, intent(out) :: order_out(:)
            integer :: p, q, tmp

            do p = 1, cnt
                order_out(p) = p
            end do
            do p = 2, cnt
                tmp = order_out(p)
                do q = p - 1, 1, -1
                    if (keys(order_out(q)) <= keys(tmp)) exit
                    order_out(q + 1) = order_out(q)
                end do
                order_out(q + 1) = tmp
            end do
        end subroutine arg_sort

        !> 未收敛最优集: 最多 k 个, 按键升序; 平行候选残差更小则升级
        subroutine merge_open(lam_v, vec_v, res_v)
            complex(dp), intent(in) :: lam_v
            complex(dp), intent(in) :: vec_v(:)
            real(dp), intent(in) :: res_v
            integer :: p, q2, pos
            real(dp) :: kv

            kv = target_key(lam_v)
            do p = 1, n_open
                proj_dup = abs(sum(conjg(open_vec(1:n_dim, p)) * vec_v(1:n_dim)))
                if (proj_dup > 0.999_dp) then
                    if (res_v < open_res(p)) then
                        open_lam(p) = lam_v
                        open_vec(1:n_dim, p) = vec_v(1:n_dim)
                        open_res(p) = res_v
                    end if
                    return
                end if
            end do
            if (n_open < k) then
                n_open = n_open + 1
                pos = n_open
            else
                ! 已满: 键最大者被替换 (若新键更小)
                pos = 1
                do p = 2, k
                    if (target_key(open_lam(p)) > target_key(open_lam(pos))) pos = p
                end do
                if (kv >= target_key(open_lam(pos))) return
            end if
            open_lam(pos) = lam_v
            open_vec(1:n_dim, pos) = vec_v(1:n_dim)
            open_res(pos) = res_v
            ! 重排 (键升序, 简单插入)
            do p = 2, n_open
                do q2 = p - 1, 1, -1
                    if (target_key(open_lam(q2)) <= target_key(open_lam(q2 + 1))) exit
                    ! swap q2, q2+1
                    lam_tmp = open_lam(q2)
                    open_lam(q2) = open_lam(q2 + 1)
                    open_lam(q2 + 1) = lam_tmp
                    vec_tmp(1:n_dim) = open_vec(1:n_dim, q2)
                    open_vec(1:n_dim, q2) = open_vec(1:n_dim, q2 + 1)
                    open_vec(1:n_dim, q2 + 1) = vec_tmp(1:n_dim)
                    tmp_r = open_res(q2)
                    open_res(q2) = open_res(q2 + 1)
                    open_res(q2 + 1) = tmp_r
                end do
            end do
        end subroutine merge_open

    end subroutine arnoldi_general_eigenpairs

end module mod_arnoldi
