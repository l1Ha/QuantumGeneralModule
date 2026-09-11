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
        integer :: i, j, k
        real(dp) :: p

        stat = 0
        z = a_in

        ! 1. Householder 三对角化
        call tred2(n, z, d, e)

        ! 2. QL 隐式位移迭代求特征值与特征向量
        call tql2(n, d, e, z, stat)
        if (stat /= 0) return

        ! 3. 特征值与特征向量升序排序 (Bubble / Insertion sort)
        do i = 1, n - 1
            k = i
            p = d(i)
            do j = i + 1, n
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
    end subroutine diag_symmetric_matrix

    !> \brief Householder 约化实对称矩阵为三对角形式
    subroutine tred2(n, z, d, e)
        integer, intent(in) :: n
        real(dp), intent(inout) :: z(n, n)
        real(dp), intent(out) :: d(n)
        real(dp), intent(out) :: e(n)

        integer :: i, j, k, l
        real(dp) :: f, g, h, hh, scale

        do i = 1, n
            d(i) = z(n, i)
        end do

        do i = n, 2, -1
            l = i - 1
            h = 0.0_dp
            scale = 0.0_dp
            if (l > 1) then
                do k = 1, l
                    scale = scale + abs(d(k))
                end do
                if (scale == 0.0_dp) then
                    e(i) = d(l)
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
                        do k = j + 1, l
                            g = g + z(k, j) * d(k)
                            e(k) = e(k) + z(k, j) * f
                        end do
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
                            z(k, j) = z(k, j) - (f * e(k) + g * d(k))
                        end do
                        d(j) = z(l, j)
                        z(i, j) = 0.0_dp
                    end do
                end if
            else
                e(i) = d(l)
            end if
            d(i) = h
        end do

        d(1) = 0.0_dp
        e(1) = 0.0_dp

        do i = 1, n
            l = i - 1
            if (d(i) /= 0.0_dp) then
                do j = 1, l
                    g = 0.0_dp
                    do k = 1, l
                        g = g + z(k, i) * z(j, k)
                    end do
                    do k = 1, l
                        z(j, k) = z(j, k) - g * z(k, i)
                    end do
                end do
            end if
            d(i) = z(i, i)
            z(i, i) = 1.0_dp
            do j = 1, l
                z(j, i) = 0.0_dp
                z(i, j) = 0.0_dp
            end do
        end do
    end subroutine tred2

    !> \brief QL 隐式位移算法求三对角矩阵的特征值与特征向量
    subroutine tql2(n, d, e, z, ierr)
        integer, intent(in) :: n
        real(dp), intent(inout) :: d(n)
        real(dp), intent(inout) :: e(n)
        real(dp), intent(inout) :: z(n, n)
        integer, intent(out) :: ierr

        integer :: i, j, k, l, m, iter
        real(dp) :: b, c, f, g, h, p, r, s

        ierr = 0
        if (n == 1) return

        do i = 2, n
            e(i - 1) = e(i)
        end do
        e(n) = 0.0_dp

        do l = 1, n
            iter = 0
            iterate_m: do
                do m = l, n - 1
                    b = abs(d(m)) + abs(d(m + 1))
                    if (abs(e(m)) + b == b) exit
                end do
                if (m == l) exit iterate_m

                if (iter >= 60) then
                    ierr = l
                    return
                end if
                iter = iter + 1

                g = (d(l + 1) - d(l)) / (2.0_dp * e(l))
                r = sqrt(g * g + 1.0_dp)
                g = d(m) - d(l) + e(l) / (g + sign(r, g))
                s = 1.0_dp
                c = 1.0_dp
                p = 0.0_dp

                do i = m - 1, l, -1
                    f = s * e(i)
                    b = c * e(i)
                    r = sqrt(f * f + g * g)
                    e(i + 1) = r
                    if (r == 0.0_dp) then
                        d(i + 1) = d(i + 1) - p
                        e(m) = 0.0_dp
                        cycle iterate_m
                    end if
                    s = f / r
                    c = g / r
                    g = d(i + 1) - p
                    r = (d(i) - g) * s + 2.0_dp * c * b
                    p = s * r
                    d(i + 1) = g + p
                    g = c * r - b
                    do k = 1, n
                        f = z(k, i + 1)
                        z(k, i + 1) = s * z(k, i) + c * f
                        z(k, i) = c * z(k, i) - s * f
                    end do
                end do
                d(l) = d(l) - p
                e(l) = g
                e(m) = 0.0_dp
            end do iterate_m
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
