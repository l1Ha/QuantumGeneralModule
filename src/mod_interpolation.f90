!> \brief 高精度科学计算数据样条插值与势能面外推模块
!> \details 提供一维自然/固定边界三次样条插值（Cubic Spline）、一阶导数（受力）、二阶导数连续计算，
!>          以及针对分子从头算势能曲面的短程指数排斥与长程范德华渐近外推例程。
!> \author LiHao
module mod_interpolation
    use mod_constants, only: dp
    implicit none
    private

    public :: spline_1d_t
    public :: spline_init
    public :: spline_eval
    public :: spline_eval_deriv
    public :: spline_eval_deriv2
    public :: spline_clean
    public :: interpolate_pes_to_grid

    ! 边界条件类型
    integer, parameter, public :: BC_NATURAL = 0  !< 自然边界条件: y''(1) = y''(n) = 0
    integer, parameter, public :: BC_CLAMPED = 1  !< 固定一阶导数边界: 指定 y'(1) 和 y'(n)

    !> \brief 一维三次样条结构体
    type :: spline_1d_t
        integer :: n = 0
        real(dp), allocatable :: x(:)   !< 原始离散横坐标 (须单调递增)
        real(dp), allocatable :: y(:)   !< 原始离散纵坐标
        real(dp), allocatable :: y2(:)  !< 各节点计算所得的二阶导数数组 M_i = y''_i
    end type spline_1d_t

contains

    !> \brief 初始化三次样条插值对象（构建三对角方程并以 Thomas 算法求解二阶导）
    !> \param[in] x 离散自变量数组 (n)，必须严格单调递增
    !> \param[in] y 离散函数值数组 (n)
    !> \param[out] spl 初始化的样条结构体
    !> \param[in] bc_type 边界条件 (可选, 默认为 BC_NATURAL)
    !> \param[in] yp1 左端点一阶导数值 (仅在 BC_CLAMPED 时生效)
    !> \param[in] ypn 右端点一阶导数值 (仅在 BC_CLAMPED 时生效)
    subroutine spline_init(x, y, spl, bc_type, yp1, ypn)
        real(dp), intent(in) :: x(:), y(:)
        type(spline_1d_t), intent(out) :: spl
        integer, intent(in), optional :: bc_type
        real(dp), intent(in), optional :: yp1, ypn

        integer :: n, i, bct
        real(dp), allocatable :: h(:), alpha(:), l(:), mu(:), z(:)
        real(dp) :: dy_left, dy_right

        n = min(size(x), size(y))
        if (n < 2) return

        bct = BC_NATURAL
        if (present(bc_type)) bct = bc_type

        spl%n = n
        allocate(spl%x(n), spl%y(n), spl%y2(n))
        spl%x = x(1:n)
        spl%y = y(1:n)

        allocate(h(n-1), alpha(n), l(n), mu(n), z(n))

        do i = 1, n - 1
            h(i) = spl%x(i+1) - spl%x(i)
            if (h(i) <= 0.0_dp) h(i) = 1.0e-14_dp  ! 保护非增序列
        end do

        ! 计算差商向量 alpha: alpha_i = 6 * ( (y_{i+1}-y_i)/h_i - (y_i-y_{i-1})/h_{i-1} )
        do i = 2, n - 1
            alpha(i) = (6.0_dp / h(i)) * (spl%y(i+1) - spl%y(i)) - &
                       (6.0_dp / h(i-1)) * (spl%y(i) - spl%y(i-1))
        end do

        ! 左边界设置
        if (bct == BC_CLAMPED .and. present(yp1)) then
            dy_left = yp1
            l(1) = 2.0_dp * h(1)
            mu(1) = 0.5_dp
            z(1) = (6.0_dp / h(1)) * (spl%y(2) - spl%y(1)) - 6.0_dp * dy_left
            z(1) = z(1) / l(1)
        else
            l(1) = 1.0_dp
            mu(1) = 0.0_dp
            z(1) = 0.0_dp
        end if

        ! Thomas 前向消元
        do i = 2, n - 1
            l(i) = 2.0_dp * (spl%x(i+1) - spl%x(i-1)) - h(i-1) * mu(i-1)
            mu(i) = h(i) / l(i)
            z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
        end do

        ! 右边界设置
        if (bct == BC_CLAMPED .and. present(ypn)) then
            dy_right = ypn
            l(n) = h(n-1) * (2.0_dp - mu(n-1))
            z(n) = (6.0_dp * dy_right - (6.0_dp / h(n-1)) * (spl%y(n) - spl%y(n-1)) - h(n-1) * z(n-1)) / l(n)
            spl%y2(n) = z(n)
        else
            l(n) = 1.0_dp
            z(n) = 0.0_dp
            spl%y2(n) = 0.0_dp
        end if

        ! Thomas 回代求解 y2 (即 M_i)
        do i = n - 1, 1, -1
            spl%y2(i) = z(i) - mu(i) * spl%y2(i+1)
        end do

        deallocate(h, alpha, l, mu, z)
    end subroutine spline_init

    !> \brief 估算样条插值函数值 y(x_val)
    pure function spline_eval(spl, x_val) result(y_val)
        type(spline_1d_t), intent(in) :: spl
        real(dp), intent(in) :: x_val
        real(dp) :: y_val

        integer :: klo, khi, k
        real(dp) :: h, a, b

        if (spl%n < 2) then
            y_val = 0.0_dp
            return
        end if

        ! 二分查找定位所在区间 [x(klo), x(khi)]
        klo = 1
        khi = spl%n
        do while (khi - klo > 1)
            k = (khi + klo) / 2
            if (spl%x(k) > x_val) then
                khi = k
            else
                klo = k
            end if
        end do

        h = spl%x(khi) - spl%x(klo)
        if (abs(h) < 1.0e-15_dp) then
            y_val = spl%y(klo)
            return
        end if

        a = (spl%x(khi) - x_val) / h
        b = (x_val - spl%x(klo)) / h

        y_val = a * spl%y(klo) + b * spl%y(khi) + &
                ((a**3 - a) * spl%y2(klo) + (b**3 - b) * spl%y2(khi)) * (h**2) / 6.0_dp
    end function spline_eval

    !> \brief 估算样条插值一阶导数值 y'(x_val) (用于计算原子力 F = -dV/dx)
    pure function spline_eval_deriv(spl, x_val) result(dy_val)
        type(spline_1d_t), intent(in) :: spl
        real(dp), intent(in) :: x_val
        real(dp) :: dy_val

        integer :: klo, khi, k
        real(dp) :: h, a, b

        if (spl%n < 2) then
            dy_val = 0.0_dp
            return
        end if

        klo = 1
        khi = spl%n
        do while (khi - klo > 1)
            k = (khi + klo) / 2
            if (spl%x(k) > x_val) then
                khi = k
            else
                klo = k
            end if
        end do

        h = spl%x(khi) - spl%x(klo)
        if (abs(h) < 1.0e-15_dp) then
            dy_val = 0.0_dp
            return
        end if

        a = (spl%x(khi) - x_val) / h
        b = (x_val - spl%x(klo)) / h

        dy_val = (spl%y(khi) - spl%y(klo)) / h - &
                 (3.0_dp * a**2 - 1.0_dp) * h * spl%y2(klo) / 6.0_dp + &
                 (3.0_dp * b**2 - 1.0_dp) * h * spl%y2(khi) / 6.0_dp
    end function spline_eval_deriv

    !> \brief 估算样条插值二阶导数值 y''(x_val)
    pure function spline_eval_deriv2(spl, x_val) result(d2y_val)
        type(spline_1d_t), intent(in) :: spl
        real(dp), intent(in) :: x_val
        real(dp) :: d2y_val

        integer :: klo, khi, k
        real(dp) :: h, a, b

        if (spl%n < 2) then
            d2y_val = 0.0_dp
            return
        end if

        klo = 1
        khi = spl%n
        do while (khi - klo > 1)
            k = (khi + klo) / 2
            if (spl%x(k) > x_val) then
                khi = k
            else
                klo = k
            end if
        end do

        h = spl%x(khi) - spl%x(klo)
        if (abs(h) < 1.0e-15_dp) then
            d2y_val = 0.0_dp
            return
        end if

        a = (spl%x(khi) - x_val) / h
        b = (x_val - spl%x(klo)) / h

        d2y_val = a * spl%y2(klo) + b * spl%y2(khi)
    end function spline_eval_deriv2

    !> \brief 清理样条动态内存
    subroutine spline_clean(spl)
        type(spline_1d_t), intent(inout) :: spl
        if (allocated(spl%x)) deallocate(spl%x)
        if (allocated(spl%y)) deallocate(spl%y)
        if (allocated(spl%y2)) deallocate(spl%y2)
        spl%n = 0
    end subroutine spline_clean

    !> \brief 从离散点集投影到连续格点，含短程指数排斥与长程范德华衰减外推
    !> \param[in] r_raw 原始从头算核间距点集 (n_raw)
    !> \param[in] v_raw 原始从头算势能点集 (n_raw)
    !> \param[in] r_grid 目标规则格点 (n_grid)
    !> \param[out] v_grid 目标格点上的平滑势能 (n_grid)
    !> \param[in] v_inf 渐近解离极限能量 (可选, 默认为 v_raw(n_raw))
    subroutine interpolate_pes_to_grid(r_raw, v_raw, r_grid, v_grid, v_inf)
        real(dp), intent(in) :: r_raw(:), v_raw(:)
        real(dp), intent(in) :: r_grid(:)
        real(dp), intent(out) :: v_grid(:)
        real(dp), intent(in), optional :: v_inf

        type(spline_1d_t) :: spl
        integer :: n_raw, n_grid, i
        real(dp) :: r1, rn, v1, vn, slope1, c6, v_asymp, r_val

        n_raw = size(r_raw)
        n_grid = size(r_grid)

        call spline_init(r_raw, v_raw, spl, BC_NATURAL)

        r1 = r_raw(1)
        rn = r_raw(n_raw)
        v1 = v_raw(1)
        vn = v_raw(n_raw)
        slope1 = spline_eval_deriv(spl, r1)

        v_asymp = vn
        if (present(v_inf)) v_asymp = v_inf

        ! 拟合 C6 范德华系数: V(rn) ~ v_asymp - C6 / rn^6
        c6 = max(0.0_dp, (v_asymp - vn) * (rn**6))

        do i = 1, n_grid
            r_val = r_grid(i)
            if (r_val < r1) then
                ! 短程排斥外推
                v_grid(i) = v1 + abs(slope1) * (r1 - r_val) + 2.0_dp * (r1 - r_val)**2
            else if (r_val > rn) then
                ! 长程范德华外推: V(R) = v_asymp - C6 / R^6
                if (c6 > 0.0_dp) then
                    v_grid(i) = v_asymp - c6 / (r_val**6)
                else
                    v_grid(i) = v_asymp
                end if
            else
                v_grid(i) = spline_eval(spl, r_val)
            end if
        end do

        call spline_clean(spl)
    end subroutine interpolate_pes_to_grid

end module mod_interpolation
