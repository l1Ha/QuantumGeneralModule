!> \brief 线性代数与快速傅里叶变换（FFT）核心算法模块
!> \details 提供实对称矩阵三对角化及本征值/特征向量求解（基于现代化 Householder QL 算法），
!>          以及 1D/2D 复数快速傅里叶变换（FFT/IFFT），全模块零外部库依赖。
!> \author LiHao
module mod_linear_algebra
    use mod_constants, only: dp, PI, TWOPI, EYE
    implicit none
    private

    public :: diag_symmetric_matrix
    public :: diag_hermitian_matrix
    public :: cholesky_decompose
    public :: cholesky_solve
    public :: qr_decompose
    public :: least_squares_qr
    public :: svd_decompose
    public :: svd_pseudoinverse
    public :: fft_1d
    public :: fft_2d
    public :: inv_real_matrix
    public :: inv_complex_matrix

contains

    !> \brief 求解 n x n 实对称矩阵的全部本征值与本征向量（按升序排列）
    !> \param[in] n 矩阵阶数
    !> \param[in] a_in 输入实对称矩阵
    !> \param[out] d 本征值向量 (n)
    !> \param[out] z 本征向量矩阵 (n, n)，第 i 列对应 d(i) 的特征向量
    !> \param[out] stat 状态码 (0 表示正常收敛)
    subroutine diag_symmetric_matrix(n, a_in, d, z, stat)
        integer, intent(in) :: n
        real(dp), intent(in) :: a_in(n, n)
        real(dp), intent(out) :: d(n)
        real(dp), intent(out) :: z(n, n)
        integer, intent(out) :: stat

        real(dp) :: e(n)

        stat = 0
        z = a_in

        ! 1. Householder 三对角化并累积正交变换
        call tred2(n, z, d, e)

        ! 2. QL 隐式位移迭代求解本征值与本征向量（自动升序排列）
        call tql2(n, d, e, z, stat)
    end subroutine diag_symmetric_matrix

    !> \brief 求解 n x n 复 Hermitian 矩阵的全部本征值与本征向量（本征值升序）
    !> \details 加倍归约: 写 A = B + iC (B 实对称, C 实反对称), 构造 2n x 2n
    !>          实对称矩阵
    !>            M = [ B  -C ]
    !>                [ C   B ] ,
    !>          其本征值恰为 A 的每个本征值二重出现。对 M 用已验证的
    !>          Householder-QL 路径对角化, 取 n 个互异本征值; 由实本征矢
    !>          [u; v] 重构复本征矢 w = u + iv (A w = lam w 严格成立)。
    !>          归一化 ||w|| = 1, 并固定最大模分量为正实数的相位约定。
    !>          代价 O((2n)^3) —— 对小中型矩阵稳健且零外部依赖;
    !>          更大规模请走 Lanczos/Arnoldi 迭代路径 (§14.1, §15.10)。
    !> \param[in] n 矩阵阶数
    !> \param[in] a_in 输入复 Hermitian 矩阵
    !> \param[out] d 本征值向量 (n, 升序)
    !> \param[out] z 复本征向量矩阵 (n, n)，第 i 列对应 d(i)
    !> \param[out] stat 状态码 (0 正常)
    subroutine diag_hermitian_matrix(n, a_in, d, z, stat)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a_in(n, n)
        real(dp), intent(out) :: d(n)
        complex(dp), intent(out) :: z(n, n)
        integer, intent(out) :: stat

        real(dp), allocatable :: m2(:, :), d2(:), z2(:, :)
        real(dp) :: nrm, b_val, c_val
        integer :: i, j, k

        allocate(m2(2 * n, 2 * n), d2(2 * n), z2(2 * n, 2 * n))

        ! 严格对称的加倍矩阵 M = [[B, -C],[C, B]]:
        !   B_ij = (Re A_ij + Re A_ji)/2 (实对称化), C_ij = Im A_ij (反对称化)
        !   对称性要求 M_{n+j, i} = M_{i, n+j} = -C_ij (= C_ji, 因 C 反对称)
        m2 = 0.0_dp
        do j = 1, n
            do i = 1, n
                b_val = 0.5_dp * (real(a_in(i, j), dp) + real(a_in(j, i), dp))
                c_val = 0.5_dp * (aimag(a_in(i, j)) - aimag(a_in(j, i)))
                m2(i, j) = b_val
                m2(n + i, n + j) = b_val
                m2(i, n + j) = -c_val
                m2(n + j, i) = -c_val
            end do
        end do

        call diag_symmetric_matrix(2 * n, m2, d2, z2, stat)
        if (stat /= 0) then
            d = 0.0_dp
            z = cmplx(0.0_dp, 0.0_dp, kind=dp)
            deallocate(m2, d2, z2)
            return
        end if

        ! 稳健选取 n 个独立本征矢:
        ! 加倍矩阵 M 的**任意**本征矢 [u; v] 都给出 A 的合法本征矢
        ! w = u + i v (A w = lam w 对 eigenspace 中任意向量成立), 故:
        !   1) 对 M 的全部 2n 列构造候选 w_j (j = 1..2n);
        !   2) 按简并组划分 (组大小 g = 2m, m 为该本征值在 A 中的重数);
        !   3) 组内 MGS 依次正交化, 塌缩 (范数远小于 1) 的候选丢弃,
        !      保留前 m 个独立向量 —— 不依赖 LAPACK 简并基序的"配对",
        !      对任意简并度稳健 (早期"隔一取一"版本在 4 重简并处因基序
        !      取到平行向量而失效, 已被本方案替代)。
        block
            integer :: kend, ii, jj, m_keep, n_keep, n_surv
            real(dp) :: nrm2, tol_deg
            complex(dp) :: proj_c
            complex(dp), allocatable :: cand(:, :)
            allocate(cand(n, 2 * n))
            do jj = 1, 2 * n
                do ii = 1, n
                    cand(ii, jj) = cmplx(z2(ii, jj), z2(n + ii, jj), kind=dp)
                end do
            end do
            n_keep = 0
            jj = 1
            do while (jj <= 2 * n)
                tol_deg = 1.0e-10_dp * max(1.0_dp, abs(d2(jj)))
                kend = jj
                do while (kend < 2 * n)
                    if (abs(d2(kend + 1) - d2(jj)) > tol_deg) exit
                    kend = kend + 1
                end do
                m_keep = (kend - jj + 1) / 2      ! A 中该本征值的重数
                n_surv = 0
                do ii = jj, kend
                    if (n_surv >= m_keep) exit
                    ! 与已保留向量正交化
                    do k = 1, n_keep
                        proj_c = sum(conjg(z(:, k)) * cand(:, ii))
                        cand(:, ii) = cand(:, ii) - proj_c * z(:, k)
                    end do
                    nrm2 = sqrt(real(sum(abs(cand(:, ii))**2), dp))
                    if (nrm2 > 1.0e-8_dp) then
                        n_keep = n_keep + 1
                        z(:, n_keep) = cand(:, ii) / nrm2
                        d(n_keep) = d2(jj)
                        n_surv = n_surv + 1
                    end if
                end do
                jj = kend + 1
            end do
            deallocate(cand)
        end block

        ! 归一化 + 相位约定 (最大模分量为正实数)
        do k = 1, n
            nrm = sqrt(real(sum(abs(z(:, k))**2), dp))
            if (nrm > 0.0_dp) z(:, k) = z(:, k) / nrm
            block
                integer :: imax
                real(dp) :: amax
                complex(dp) :: ph
                imax = 1
                amax = abs(z(1, k))
                do i = 2, n
                    if (abs(z(i, k)) > amax) then
                        amax = abs(z(i, k)); imax = i
                    end if
                end do
                if (amax > 0.0_dp) then
                    ph = z(imax, k) / amax
                    z(:, k) = z(:, k) * conjg(ph)
                end if
            end block
        end do

        deallocate(m2, d2, z2)
    end subroutine diag_hermitian_matrix

    !> \brief 安全计算 sqrt(a^2 + b^2) 避免数值溢出
    pure function pythag(a, b) result(p)
        real(dp), intent(in) :: a, b
        real(dp) :: p, absa, absb
        absa = abs(a)
        absb = abs(b)
        if (absa > absb) then
            p = absa * sqrt(1.0_dp + (absb / absa)**2)
        else if (absb > 0.0_dp) then
            p = absb * sqrt(1.0_dp + (absa / absb)**2)
        else
            p = 0.0_dp
        end if
    end function pythag

    !> \brief Householder 约化实对称矩阵为三对角形式并累积正交变换矩阵 (EISPACK TRED2)
    subroutine tred2(n, z, d, e)
        integer, intent(in) :: n
        real(dp), intent(inout) :: z(n, n)
        real(dp), intent(out) :: d(n)
        real(dp), intent(out) :: e(n)

        integer :: i, j, k, l, ii, jp1
        real(dp) :: f, g, h, hh, scale

        do i = 1, n
            d(i) = z(n, i)
        end do

        if (n == 1) then
            d(1) = z(1, 1)
            z(1, 1) = 1.0_dp
            e(1) = 0.0_dp
            return
        end if

        do ii = 2, n
            i = n + 2 - ii
            l = i - 1
            h = 0.0_dp
            scale = 0.0_dp
            if (l < 2) then
                e(i) = d(l)
                do j = 1, l
                    d(j) = z(l, j)
                    z(i, j) = 0.0_dp
                    z(j, i) = 0.0_dp
                end do
                d(i) = h
                cycle
            end if

            do k = 1, l
                scale = scale + abs(d(k))
            end do

            if (scale <= 1.0e-35_dp) then
                e(i) = d(l)
                do j = 1, l
                    d(j) = z(l, j)
                    z(i, j) = 0.0_dp
                    z(j, i) = 0.0_dp
                end do
            else
                do k = 1, l
                    d(k) = d(k) / scale
                    h = h + d(k) * d(k)
                end do
                f = d(l)
                g = -sign(sqrt(h), f)
                e(i) = scale * g
                h = h - f * g
                d(l) = f - g

                do j = 1, l
                    e(j) = 0.0_dp
                end do

                do j = 1, l
                    f = d(j)
                    z(j, i) = f
                    g = e(j) + z(j, j) * f
                    jp1 = j + 1
                    if (l >= jp1) then
                        do k = jp1, l
                            g = g + z(k, j) * d(k)
                            e(k) = e(k) + z(k, j) * f
                        end do
                    end if
                    e(j) = g
                end do

                f = 0.0_dp
                do j = 1, l
                    e(j) = e(j) / h
                    f = f + e(j) * d(j)
                end do
                hh = f / (h + h)

                do j = 1, l
                    e(j) = e(j) - hh * d(j)
                end do

                do j = 1, l
                    f = d(j)
                    g = e(j)
                    do k = j, l
                        z(k, j) = z(k, j) - f * e(k) - g * d(k)
                    end do
                    d(j) = z(l, j)
                    z(i, j) = 0.0_dp
                end do
            end if
            d(i) = h
        end do

        ! 累积正交变换矩阵
        do i = 2, n
            l = i - 1
            z(n, l) = z(l, l)
            z(l, l) = 1.0_dp
            h = d(i)
            if (abs(h) > 1.0e-35_dp) then
                do k = 1, l
                    d(k) = z(k, i) / h
                end do
                do j = 1, l
                    g = 0.0_dp
                    do k = 1, l
                        g = g + z(k, i) * z(k, j)
                    end do
                    do k = 1, l
                        z(k, j) = z(k, j) - g * d(k)
                    end do
                end do
            end if
            do k = 1, l
                z(k, i) = 0.0_dp
            end do
        end do

        do i = 1, n
            d(i) = z(n, i)
            z(n, i) = 0.0_dp
        end do
        z(n, n) = 1.0_dp
        e(1) = 0.0_dp
    end subroutine tred2

    !> \brief QL 隐式位移算法求实对称三对角矩阵本征值与本征向量 (EISPACK TQL2)
    subroutine tql2(n, d, e, z, ierr)
        integer, intent(in) :: n
        real(dp), intent(inout) :: d(n), e(n), z(n, n)
        integer, intent(out) :: ierr

        integer :: i, j, k, l, m, ii, l1, l2, mml
        real(dp) :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1
        real(dp) :: eps_prec

        ierr = 0
        if (n == 1) return
        eps_prec = epsilon(1.0_dp)

        do i = 2, n
            e(i - 1) = e(i)
        end do
        f = 0.0_dp
        tst1 = 0.0_dp
        e(n) = 0.0_dp

        do l = 1, n
            j = 0
            h = abs(d(l)) + abs(e(l))
            if (tst1 < h) tst1 = h

            ! 查找极小次对角元素
            do m = l, n
                if (abs(e(m)) <= eps_prec * tst1) exit
            end do

            if (m /= l) then
                do
                    if (j == 60) then
                        ierr = l
                        return
                    end if
                    j = j + 1

                    l1 = l + 1
                    l2 = l1 + 1
                    g = d(l)
                    p = (d(l1) - g) / (2.0_dp * e(l))
                    r = pythag(p, 1.0_dp)
                    d(l) = e(l) / (p + sign(r, p))
                    d(l1) = e(l) * (p + sign(r, p))
                    dl1 = d(l1)
                    h = g - d(l)
                    if (l2 <= n) then
                        do i = l2, n
                            d(i) = d(i) - h
                        end do
                    end if
                    f = f + h

                    ! QL 变换
                    p = d(m)
                    c = 1.0_dp
                    c2 = c
                    c3 = c
                    el1 = e(l1)
                    s = 0.0_dp
                    s2 = s
                    mml = m - l

                    do ii = 1, mml
                        c3 = c2
                        c2 = c
                        s2 = s
                        i = m - ii
                        g = c * e(i)
                        h = c * p
                        r = pythag(p, e(i))
                        e(i + 1) = s * r
                        s = e(i) / r
                        c = p / r
                        p = c * d(i) - s * g
                        d(i + 1) = h + s * (c * g + s * d(i))
                        do k = 1, n
                            h = z(k, i + 1)
                            z(k, i + 1) = s * z(k, i) + c * h
                            z(k, i) = c * z(k, i) - s * h
                        end do
                    end do

                    p = -s * s2 * c3 * el1 * e(l) / dl1
                    e(l) = s * p
                    d(l) = c * p
                    if (abs(e(l)) <= eps_prec * tst1) exit
                end do
            end if
            d(l) = d(l) + f
        end do

        ! 本征值与本征向量升序排列
        do ii = 2, n
            i = ii - 1
            k = i
            p = d(i)
            do j = ii, n
                if (d(j) < p) then
                    k = j
                    p = d(j)
                end if
            end do
            if (k /= i) then
                d(k) = d(i)
                d(i) = p
                do j = 1, n
                    p = z(j, i)
                    z(j, i) = z(j, k)
                    z(j, k) = p
                end do
            end if
        end do
    end subroutine tql2

    !> \brief 一维快速傅里叶变换 (Cooley-Tukey 算法，输入长度须为 2 的幂)
    !> \param[inout] data_vec 复数数据向量
    !> \param[in] isign -1: 正变换 (t -> omega), +1: 逆变换 (omega -> t)
    subroutine fft_1d(data_vec, isign)
        complex(dp), intent(inout) :: data_vec(:)
        integer, intent(in) :: isign

        integer :: n, i, j, m, mmax, istep
        real(dp) :: theta, wtemp, wpr, wpi, wr, wi
        complex(dp) :: temp

        n = size(data_vec)
        ! 校验是否为 2 的幂次
        if (iand(n, n - 1) /= 0) then
            ! 非 2 的幂次调用通用 DFT
            call dft_slow(data_vec, isign)
            return
        end if

        ! 位反转重排
        j = 1
        do i = 1, n
            if (j > i) then
                temp = data_vec(j)
                data_vec(j) = data_vec(i)
                data_vec(i) = temp
            end if
            m = n / 2
            do while (m >= 2 .and. j > m)
                j = j - m
                m = m / 2
            end do
            j = j + m
        end do

        ! Danielson-Lanczos 蝶形运算
        mmax = 1
        do while (n > mmax)
            istep = 2 * mmax
            theta = -real(isign, dp) * (TWOPI / real(istep, dp))
            wtemp = sin(0.5_dp * theta)
            wpr = -2.0_dp * wtemp * wtemp
            wpi = sin(theta)
            wr = 1.0_dp
            wi = 0.0_dp
            do m = 1, mmax
                do i = m, n, istep
                    j = i + mmax
                    temp = cmplx(wr, wi, kind=dp) * data_vec(j)
                    data_vec(j) = data_vec(i) - temp
                    data_vec(i) = data_vec(i) + temp
                end do
                wtemp = wr
                wr = wr * wpr - wi * wpi + wr
                wi = wi * wpr + wtemp * wpi + wi
            end do
            mmax = istep
        end do

        ! 逆变换归一化
        if (isign > 0) then
            data_vec = data_vec / real(n, dp)
        end if
    end subroutine fft_1d

    !> \brief 通用非 2 的幂次离散傅里叶变换 (保底回退例程)
    subroutine dft_slow(data_vec, isign)
        complex(dp), intent(inout) :: data_vec(:)
        integer, intent(in) :: isign
        integer :: n, j, k
        complex(dp), allocatable :: temp(:)
        real(dp) :: angle

        n = size(data_vec)
        allocate(temp(n))
        temp = data_vec

        do k = 1, n
            data_vec(k) = (0.0_dp, 0.0_dp)
            do j = 1, n
                angle = -real(isign, dp) * TWOPI * real((j - 1) * (k - 1), dp) / real(n, dp)
                data_vec(k) = data_vec(k) + temp(j) * exp(EYE * angle)
            end do
        end do

        if (isign > 0) then
            data_vec = data_vec / real(n, dp)
        end if
        deallocate(temp)
    end subroutine dft_slow

    !> \brief 二维快速傅里叶变换（转置+连续列访问高性能算法，零内层动态内存分配，缓存友好）
    subroutine fft_2d(data_mat, isign)
        complex(dp), intent(inout) :: data_mat(:, :)
        integer, intent(in) :: isign
        integer :: nx, ny, i
        complex(dp), allocatable :: temp(:, :)

        nx = size(data_mat, 1)
        ny = size(data_mat, 2)
        allocate(temp(ny, nx))

        ! 1. 对所有列做 1D FFT（沿 Fortran 列主序连续内存访问）
        !$omp parallel do schedule(static)
        do i = 1, ny
            call fft_1d(data_mat(:, i), isign)
        end do
        !$omp end parallel do

        ! 2. 矩阵转置：行转化为列 (分块 + OpenMP 块并行)
        call transpose_blocked(data_mat, temp)

        ! 3. 对转置后的矩阵所有列做 1D FFT（相当于对原始矩阵的行做变换，再次保证纯连续内存步长）
        !$omp parallel do schedule(static)
        do i = 1, nx
            call fft_1d(temp(:, i), isign)
        end do
        !$omp end parallel do

        ! 4. 转置还原 (分块 + OpenMP 块并行)
        call transpose_blocked(temp, data_mat)
        deallocate(temp)
    end subroutine fft_2d

    !> \brief 分块转置 b = a^T (缓存友好; 编译 OpenMP 时块循环并行)
    !> \details fft_2d 的转置步骤原为串行 `transpose()` 内在函数, 在
    !>          OpenMP 扩展测量中成为 Amdahl 瓶颈 (1->20 线程仅 2.7 倍);
    !>          分块 (block size 32) 转置同时改善缓存局部性, 并暴露
    !>          块级并行维度 (collapse(2))。
    subroutine transpose_blocked(a, b)
        complex(dp), intent(in) :: a(:, :)
        complex(dp), intent(out) :: b(:, :)
        integer, parameter :: bs = 32
        integer :: i0, j0, i, j, nx, ny, ni, nj

        nx = size(a, 1)
        ny = size(a, 2)
        !$omp parallel do collapse(2) schedule(static) private(i0, j0, i, j, ni, nj)
        do j0 = 1, nx, bs
            do i0 = 1, ny, bs
                nj = min(bs, nx - j0 + 1)
                ni = min(bs, ny - i0 + 1)
                do j = j0, j0 + nj - 1
                    do i = i0, i0 + ni - 1
                        b(i, j) = a(j, i)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine transpose_blocked

    !> \brief 通用实方阵求逆（Gauss-Jordan 全主元消去法，零外部依赖）
    subroutine inv_real_matrix(n, a, a_inv, stat)
        integer, intent(in)   :: n
        real(dp), intent(in)  :: a(n, n)
        real(dp), intent(out) :: a_inv(n, n)
        integer, intent(out)  :: stat

        real(dp) :: aug(n, 2 * n), temp, pivot, factor
        integer  :: i, j, k, max_row

        stat = 0
        aug(:, 1:n) = a
        aug(:, n + 1:2 * n) = 0.0_dp
        do i = 1, n
            aug(i, n + i) = 1.0_dp
        end do

        do i = 1, n
            max_row = i
            pivot = abs(aug(i, i))
            do k = i + 1, n
                if (abs(aug(k, i)) > pivot) then
                    pivot = abs(aug(k, i))
                    max_row = k
                end if
            end do

            if (pivot < 1.0e-30_dp) then
                stat = -1
                a_inv = 0.0_dp
                return
            end if

            if (max_row /= i) then
                do j = 1, 2 * n
                    temp = aug(i, j)
                    aug(i, j) = aug(max_row, j)
                    aug(max_row, j) = temp
                end do
            end if

            factor = aug(i, i)
            aug(i, :) = aug(i, :) / factor

            do k = 1, n
                if (k /= i) then
                    factor = aug(k, i)
                    aug(k, :) = aug(k, :) - factor * aug(i, :)
                end if
            end do
        end do

        a_inv = aug(:, n + 1:2 * n)
    end subroutine inv_real_matrix

    !> \brief 通用复方阵求逆（Gauss-Jordan 全主元消去法，零外部依赖）
    subroutine inv_complex_matrix(n, a, a_inv, stat)
        integer, intent(in)      :: n
        complex(dp), intent(in)  :: a(n, n)
        complex(dp), intent(out) :: a_inv(n, n)
        integer, intent(out)     :: stat

        complex(dp) :: aug(n, 2 * n), temp, factor
        real(dp)    :: pivot
        integer     :: i, j, k, max_row

        stat = 0
        aug(:, 1:n) = a
        aug(:, n + 1:2 * n) = (0.0_dp, 0.0_dp)
        do i = 1, n
            aug(i, n + i) = (1.0_dp, 0.0_dp)
        end do

        do i = 1, n
            max_row = i
            pivot = abs(aug(i, i))
            do k = i + 1, n
                if (abs(aug(k, i)) > pivot) then
                    pivot = abs(aug(k, i))
                    max_row = k
                end if
            end do

            if (pivot < 1.0e-30_dp) then
                stat = -1
                a_inv = (0.0_dp, 0.0_dp)
                return
            end if

            if (max_row /= i) then
                do j = 1, 2 * n
                    temp = aug(i, j)
                    aug(i, j) = aug(max_row, j)
                    aug(max_row, j) = temp
                end do
            end if

            factor = aug(i, i)
            aug(i, :) = aug(i, :) / factor

            do k = 1, n
                if (k /= i) then
                    factor = aug(k, i)
                    aug(k, :) = aug(k, :) - factor * aug(i, :)
                end if
            end do
        end do

        a_inv = aug(:, n + 1:2 * n)
    end subroutine inv_complex_matrix

    !> \brief Cholesky 分解 A = L L^T (A 对称正定, L 下三角)
    !> \details 经典平方根算法逐列递推; 对角主元非正时判定非正定 (stat = 1)。
    !> \param[in] n 阶数
    !> \param[in] a 输入对称矩阵 (仅用下三角)
    !> \param[out] l 下三角因子
    !> \param[out] stat 0 正常; 1 非正定
    subroutine cholesky_decompose(n, a, l, stat)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n)
        real(dp), intent(out) :: l(n, n)
        integer, intent(out) :: stat

        real(dp) :: ssum
        integer :: i, j

        stat = 0
        l = 0.0_dp
        do j = 1, n
            ssum = a(j, j) - sum(l(j, 1:j - 1)**2)
            if (ssum <= 0.0_dp) then
                stat = 1
                return
            end if
            l(j, j) = sqrt(ssum)
            do i = j + 1, n
                ssum = a(i, j) - sum(l(i, 1:j - 1) * l(j, 1:j - 1))
                l(i, j) = ssum / l(j, j)
            end do
        end do
    end subroutine cholesky_decompose

    !> \brief 用 Cholesky 因子解 A x = b (L y = b 前代; L^T x = y 回代)
    subroutine cholesky_solve(n, l, b, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: l(n, n), b(n)
        real(dp), intent(out) :: x(n)

        real(dp) :: y(n), ssum
        integer :: i, k

        do i = 1, n
            ssum = b(i)
            do k = 1, i - 1
                ssum = ssum - l(i, k) * y(k)
            end do
            y(i) = ssum / l(i, i)
        end do
        do i = n, 1, -1
            ssum = y(i)
            do k = i + 1, n
                ssum = ssum - l(k, i) * x(k)
            end do
            x(i) = ssum / l(i, i)
        end do
    end subroutine cholesky_solve

    !> \brief 薄 QR 分解 A = Q R (Householder 反射, m >= n)
    !> \details A (m x n) -> Q (m x n, 列正交), R (n x n, 上三角)。
    !>          反射向量在约化过程中保存, Q 由 H_1 H_2 ... H_k 作用单位阵
    !>          逆序累计; 正交性对病态矩阵仍保持机器精度 (对比正规方程
    !>          法失去约 2*cond 位数字, 见测试组 2)。
    !> \param[in] m,n 维数 (m >= n)
    !> \param[in] a 输入矩阵
    !> \param[out] q 正交因子 (m x n)
    !> \param[out] r 上三角因子 (n x n)
    !> \param[out] stat 0 正常; 1 m < n
    subroutine qr_decompose(m, n, a, q, r, stat)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: a(m, n)
        real(dp), intent(out) :: q(m, n), r(n, n)
        integer, intent(out) :: stat

        real(dp), allocatable :: w(:, :), vs(:, :), v(:)
        real(dp) :: alpha, vTv
        integer :: j, k, i

        stat = 0
        if (m < n) then
            stat = 1
            q = 0.0_dp; r = 0.0_dp
            return
        end if

        allocate (w(m, n), vs(m, n), v(m))
        w = a
        vs = 0.0_dp
        do j = 1, n
            v = 0.0_dp
            v(j:m) = w(j:m, j)
            alpha = sqrt(sum(v(j:m)**2))
            if (alpha > 0.0_dp) then
                if (v(j) >= 0.0_dp) then
                    v(j) = v(j) + alpha
                else
                    v(j) = v(j) - alpha
                end if
                vTv = sum(v(j:m)**2)
                if (vTv > 0.0_dp) then
                    vs(:, j) = v
                    do k = j, n
                        alpha = 2.0_dp * sum(v(j:m) * w(j:m, k)) / vTv
                        w(j:m, k) = w(j:m, k) - alpha * v(j:m)
                    end do
                end if
            end if
        end do

        r = 0.0_dp
        do j = 1, n
            do k = j, n
                r(j, k) = w(j, k)
            end do
        end do

        ! Q = H_1 H_2 ... H_n 作用于单位阵: 逆序应用反射
        q = 0.0_dp
        do i = 1, n
            q(i, i) = 1.0_dp
        end do
        do j = n, 1, -1
            v = vs(:, j)
            vTv = sum(v(j:m)**2)
            if (vTv > 0.0_dp) then
                do k = 1, n
                    alpha = 2.0_dp * sum(v(j:m) * q(j:m, k)) / vTv
                    q(j:m, k) = q(j:m, k) - alpha * v(j:m)
                end do
            end if
        end do

        deallocate (w, vs, v)
    end subroutine qr_decompose

    !> \brief 用薄 QR 解最小二乘 min ||A x - b|| (m >= n)
    !> \param[in] m,n 维数
    !> \param[in] a 设计矩阵 (m x n)
    !> \param[in] b 右端 (m)
    !> \param[out] x 解 (n)
    !> \param[out] resid 残差范数 ||A x - b||
    !> \param[out] stat 0 正常; 其它来自 qr_decompose
    subroutine least_squares_qr(m, n, a, b, x, resid, stat)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: a(m, n), b(m)
        real(dp), intent(out) :: x(n), resid
        integer, intent(out) :: stat

        real(dp), allocatable :: q(:, :), r(:, :), qtb(:), y(:)
        real(dp) :: ssum
        integer :: i, k
        logical :: singular

        allocate (q(m, n), r(n, n), qtb(n), y(n))
        call qr_decompose(m, n, a, q, r, stat)
        if (stat /= 0) then
            x = 0.0_dp; resid = 0.0_dp
            deallocate (q, r, qtb, y)
            return
        end if
        qtb = matmul(transpose(q), b)

        ! R x = Q^T b 回代 (检测零主元 -> 奇异, 用最小范数补充: 置零)
        singular = .false.
        y = 0.0_dp
        do i = n, 1, -1
            ssum = qtb(i)
            do k = i + 1, n
                ssum = ssum - r(i, k) * y(k)
            end do
            if (abs(r(i, i)) > 1.0e-300_dp) then
                y(i) = ssum / r(i, i)
            else
                y(i) = 0.0_dp
                singular = .true.
            end if
        end do
        x = y
        resid = sqrt(sum((matmul(a, x) - b)**2))
        if (singular) stat = 2
        deallocate (q, r, qtb, y)
    end subroutine least_squares_qr

    !> \brief 单边 Jacobi SVD: A = U diag(s) V^T (m >= n)
    !> \details 对 A 的列做 Jacobi 旋转正交化, 收敛后列范数即奇异值;
    !>          零奇异值对应 U 列由标准基 Gram-Schmidt 补全, 保证 U 列正交。
    !> \param[in] m,n 维数 (m >= n)
    !> \param[in] a 输入 (m x n)
    !> \param[out] u 左奇异向量 (m x n)
    !> \param[out] sv 奇异值 (n, 降序, 非负)
    !> \param[out] vt 右奇异向量转置 (n x n)
    !> \param[out] stat 0 正常; 1 m < n; 2 未在迭代上限内收敛
    subroutine svd_decompose(m, n, a, u, sv, vt, stat)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: a(m, n)
        real(dp), intent(out) :: u(m, n), sv(n), vt(n, n)
        integer, intent(out) :: stat

        real(dp), allocatable :: b(:, :), v(:, :)
        real(dp) :: alpha, beta, gamma, zeta, c, s_rot, tau, t
        real(dp) :: eps_conv, nrm
        integer :: p, q, i, j, sweep, kmax, k
        real(dp), allocatable :: tmpcol(:)
        logical :: changed

        stat = 0
        if (m < n) then
            stat = 1
            u = 0.0_dp; sv = 0.0_dp; vt = 0.0_dp
            return
        end if

        allocate (b(m, n), v(n, n), tmpcol(m))
        b = a
        v = 0.0_dp
        do i = 1, n
            v(i, i) = 1.0_dp
        end do

        eps_conv = 1.0e-15_dp
        kmax = 60
        do sweep = 1, kmax
            changed = .false.
            do p = 1, n - 1
                do q = p + 1, n
                    alpha = sum(b(:, p)**2)
                    beta = sum(b(:, q)**2)
                    gamma = sum(b(:, p) * b(:, q))
                    if (abs(gamma) <= eps_conv * sqrt(alpha * beta)) cycle
                    changed = .true.
                    zeta = (beta - alpha) / (2.0_dp * gamma)
                    t = sign(1.0_dp, zeta) / (abs(zeta) + sqrt(1.0_dp + zeta**2))
                    c = 1.0_dp / sqrt(1.0_dp + t**2)
                    s_rot = c * t
                    ! 旋转列 p,q: b' = b J; 同时累计 V
                    do i = 1, m
                        tau = b(i, p)
                        b(i, p) = c * tau - s_rot * b(i, q)
                        b(i, q) = s_rot * tau + c * b(i, q)
                    end do
                    do i = 1, n
                        tau = v(i, p)
                        v(i, p) = c * tau - s_rot * v(i, q)
                        v(i, q) = s_rot * tau + c * v(i, q)
                    end do
                end do
            end do
            if (.not. changed) exit
            if (sweep == kmax) stat = 2
        end do

        ! 奇异值 = 列范数; 按降序排序 (U, V 列同步)
        do i = 1, n
            sv(i) = sqrt(sum(b(:, i)**2))
        end do
        do p = 1, n - 1
            do q = p + 1, n
                if (sv(q) > sv(p)) then
                    t = sv(p); sv(p) = sv(q); sv(q) = t
                    tmpcol = b(:, p); b(:, p) = b(:, q); b(:, q) = tmpcol
                    tmpcol(1:n) = v(1:n, p); v(1:n, p) = v(1:n, q); v(1:n, q) = tmpcol(1:n)
                end if
            end do
        end do

        ! U 列 = b 列归一化; 零奇异值 -> 标准基 Gram-Schmidt 补全
        u = 0.0_dp
        do j = 1, n
            if (sv(j) > 1.0e-14_dp * max(1.0_dp, sv(1))) then
                u(:, j) = b(:, j) / sv(j)
            else
                ! 补全: 从标准基中找与已定列正交的向量
                do k = 1, m
                    tmpcol = 0.0_dp
                    tmpcol(k) = 1.0_dp
                    do i = 1, j - 1
                        tmpcol = tmpcol - sum(u(:, i) * tmpcol) * u(:, i)
                    end do
                    nrm = sqrt(sum(tmpcol**2))
                    if (nrm > 1.0e-8_dp) then
                        u(:, j) = tmpcol / nrm
                        exit
                    end if
                end do
            end if
        end do
        vt = transpose(v)
        deallocate (b, v, tmpcol)
    end subroutine svd_decompose

    !> \brief SVD 伪逆最小范数最小二乘解 x = V Sigma^+ U^T b
    !> \details 秩亏问题的标准稳健解法 (正规方程/QR 在此失效):
    !>          零奇异值方向截断, 给出最小二乘中的最小范数解。
    !> \param[in] m,n 维数
    !> \param[in] a,b 设计矩阵与右端
    !> \param[out] x 解 (n)
    !> \param[out] resid 残差范数
    !> \param[in] rtol (可选) 奇异值截断相对阈值, 默认 1e-12
    subroutine svd_pseudoinverse(m, n, a, b, x, resid, rtol)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: a(m, n), b(m)
        real(dp), intent(out) :: x(n), resid
        real(dp), intent(in), optional :: rtol

        real(dp), allocatable :: u(:, :), vt(:, :), sv(:), c(:)
        real(dp) :: tol
        integer :: i, stat_s

        allocate (u(m, n), vt(n, n), sv(n), c(n))
        call svd_decompose(m, n, a, u, sv, vt, stat_s)
        tol = 1.0e-12_dp
        if (present(rtol)) tol = rtol
        c = matmul(transpose(u), b)
        x = 0.0_dp
        do i = 1, n
            if (sv(i) > tol * max(1.0_dp, sv(1))) x = x + (c(i) / sv(i)) * vt(i, :)
        end do
        resid = sqrt(sum((matmul(a, x) - b)**2))
        deallocate (u, vt, sv, c)
    end subroutine svd_pseudoinverse

end module mod_linear_algebra
