!> \brief 线性代数与快速傅里叶变换（FFT）核心算法模块
!> \details 提供实对称矩阵三对角化及本征值/特征向量求解（基于现代化 Householder QL 算法），
!>          以及 1D/2D 复数快速傅里叶变换（FFT/IFFT），全模块零外部库依赖。
!> \author LiHao
module mod_linear_algebra
    use mod_constants, only: dp, PI, TWOPI, EYE
    implicit none
    private

    public :: diag_symmetric_matrix
    public :: fft_1d
    public :: fft_2d

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

    !> \brief 安全计算 sqrt(a^2 + b^2) 避免数值溢出
    pure function pythag(a, b) result(p)
        real(dp), intent(in) :: a, b
        real(dp) :: p, absa, absb
        absa = abs(a)
        absb = abs(b)
        if (absa > absb) then
            p = absa * sqrt(1.0_dp + (absb / absa)**2)
        else if (absb /= 0.0_dp) then
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

            if (scale == 0.0_dp) then
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
            if (h /= 0.0_dp) then
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
        real(dp) :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2

        ierr = 0
        if (n == 1) return

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
                tst2 = tst1 + abs(e(m))
                if (tst2 == tst1) exit
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
                    el1 = e(l1)
                    s = 0.0_dp
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
                    tst2 = tst1 + abs(e(l))
                    if (tst2 <= tst1) exit
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

    !> \brief 二维快速傅里叶变换
    subroutine fft_2d(data_mat, isign)
        complex(dp), intent(inout) :: data_mat(:, :)
        integer, intent(in) :: isign
        integer :: nx, ny, ix, iy
        complex(dp), allocatable :: row(:), col(:)

        nx = size(data_mat, 1)
        ny = size(data_mat, 2)

        allocate(row(ny))
        allocate(col(nx))

        ! 对每行做 1D FFT
        do ix = 1, nx
            row = data_mat(ix, :)
            call fft_1d(row, isign)
            data_mat(ix, :) = row
        end do

        ! 对每列做 1D FFT
        do iy = 1, ny
            col = data_mat(:, iy)
            call fft_1d(col, isign)
            data_mat(:, iy) = col
        end do

        deallocate(row, col)
    end subroutine fft_2d

end module mod_linear_algebra
