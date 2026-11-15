!> \brief 高精度科学计算数据样条插值与势能面外推模块
!> \details 提供一维自然/固定边界三次样条插值（Cubic Spline）、一阶导数（受力）、二阶导数连续计算，
!>          二维双三次张量积样条（PES 插值/梯度，§15.10 路线图项），
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
    public :: spline_2d_t
    public :: spline_2d_init
    public :: spline_2d_eval
    public :: spline_2d_grad
    public :: spline_2d_clean
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

    !> \brief 二维双三次张量积样条结构体 (PES 插值)
    !> \details 在规则网格 (x_i, y_j) 上存储 z 值与三个导数数组:
    !>          zxx = d2z/dx2, zyy = d2z/dy2, zxxyy = d4z/(dx2 dy2),
    !>          由沿两个方向的一维三次样条反复构造 (张量积):
    !>            (1) 逐行 (固定 x_i) 沿 y 做 1D 样条 -> zyy(i, :)
    !>            (2) 逐列 (固定 y_j) 沿 x 做 1D 样条 -> zxx(:, j)
    !>            (3) 逐列沿 x 对 zyy 做 1D 样条 -> zxxyy(:, j)
    !>          胞元内用 Hermite 型基函数 (T 线性帽, U 三次样条基):
    !>            S(x,y) = sum_{p,q} [ z T_p T_q + zxx U_p T_q + zyy T_p U_q + zxxyy U_p U_q ]
    !>          保证节点处精确插值 (S = z, d2S/dx2 = zxx, d2S/dy2 = zyy,
    !>          d4S/dx2dy2 = zxxyy) 与 C1 连续性; 自然边界条件。
    type :: spline_2d_t
        integer :: nx = 0, ny = 0
        real(dp), allocatable :: x(:)              !< 网格 x (nx, 严格递增)
        real(dp), allocatable :: y(:)              !< 网格 y (ny, 严格递增)
        real(dp), allocatable :: z(:, :)           !< 网格值 (nx, ny)
        real(dp), allocatable :: zxx(:, :)         !< d2z/dx2 节点值
        real(dp), allocatable :: zyy(:, :)         !< d2z/dy2 节点值
        real(dp), allocatable :: zxxyy(:, :)       !< d4z/(dx2 dy2) 节点值
    end type spline_2d_t

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

    !> \brief 初始化二维双三次张量积样条 (PES 插值)
    !> \details 张量积构造: (1) 逐行沿 y 做 1D 样条得 zyy;
    !>          (2) 逐列沿 x 做 1D 样条得 zxx; (3) 逐列沿 x 对 zyy
    !>          做 1D 样条得混合导数 zxxyy。全部复用 1D 自然样条核。
    !> \param[in] x 网格 x (nx, 严格递增)
    !> \param[in] y 网格 y (ny, 严格递增)
    !> \param[in] z 网格值 (nx, ny)
    !> \param[out] spl 二维样条结构体
    !> \param[out] stat 状态码 (0 正常)
    subroutine spline_2d_init(x, y, z, spl, stat)
        real(dp), intent(in) :: x(:), y(:), z(:, :)
        type(spline_2d_t), intent(out) :: spl
        integer, intent(out) :: stat

        type(spline_1d_t) :: s1
        integer :: nx, ny, i, j
        real(dp), allocatable :: row(:), col(:)

        stat = 0
        nx = size(x)
        ny = size(y)
        if (size(z, 1) /= nx .or. size(z, 2) /= ny .or. nx < 2 .or. ny < 2) then
            stat = 1
            return
        end if

        spl%nx = nx
        spl%ny = ny
        allocate (spl%x(nx), spl%y(ny), spl%z(nx, ny))
        allocate (spl%zxx(nx, ny), spl%zyy(nx, ny), spl%zxxyy(nx, ny))
        spl%x = x
        spl%y = y
        spl%z = z

        allocate (row(ny), col(nx))

        ! (1) 逐行 (固定 x_i) 沿 y: zyy
        do i = 1, nx
            call spline_init(y, z(i, :), s1, BC_NATURAL)
            do j = 1, ny
                spl%zyy(i, j) = spline_eval_deriv2(s1, y(j))
            end do
            call spline_clean(s1)
        end do

        ! (2) 逐列 (固定 y_j) 沿 x: zxx
        do j = 1, ny
            call spline_init(x, z(:, j), s1, BC_NATURAL)
            do i = 1, nx
                spl%zxx(i, j) = spline_eval_deriv2(s1, x(i))
            end do
            call spline_clean(s1)
        end do

        ! (3) 逐列沿 x 对 zyy 做样条: zxxyy (混合导数)
        do j = 1, ny
            row = spl%zyy(:, j)
            call spline_init(x, row, s1, BC_NATURAL)
            do i = 1, nx
                spl%zxxyy(i, j) = spline_eval_deriv2(s1, x(i))
            end do
            call spline_clean(s1)
        end do

        deallocate (row, col)
    end subroutine spline_2d_init

    !> \brief 胞元定位辅助: 返回包含 xv 的胞元左端索引 (钳制到边界胞元)
    pure function locate_cell(g, n, xv) result(idx)
        real(dp), intent(in) :: g(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: xv
        integer :: idx
        integer :: lo, hi, mid
        if (xv <= g(1)) then
            idx = 1
            return
        end if
        if (xv >= g(n)) then
            idx = n - 1
            return
        end if
        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (xv >= g(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        idx = lo
    end function locate_cell

    !> \brief Hermite 型基函数 (T 线性帽, U 三次样条基) 及一阶导
    !> \details 胞元 [g_i, g_{i+1}], 参数 xv:
    !>          T_i = (g_{i+1}-xv)/h, T_{i+1} = (xv-g_i)/h;
    !>          U_i = ((g_{i+1}-xv)^3/h - h(g_{i+1}-xv))/6,
    !>          U_{i+1} = ((xv-g_i)^3/h - h(xv-g_i))/6;
    !>          导数 dT_i = -1/h, dT_{i+1} = 1/h,
    !>          dU_i = (h - 3(g_{i+1}-xv)^2/h)/6,
    !>          dU_{i+1} = (3(xv-g_i)^2/h - h)/6。
    pure subroutine hermite_basis(g_i, g_ip1, xv, t0, t1, u0, u1, dt0, dt1, du0, du1)
        real(dp), intent(in) :: g_i, g_ip1, xv
        real(dp), intent(out) :: t0, t1, u0, u1, dt0, dt1, du0, du1
        real(dp) :: h, d0, d1
        h = g_ip1 - g_i
        d0 = g_ip1 - xv
        d1 = xv - g_i
        t0 = d0 / h
        t1 = d1 / h
        u0 = (d0**3 / h - h * d0) / 6.0_dp
        u1 = (d1**3 / h - h * d1) / 6.0_dp
        dt0 = -1.0_dp / h
        dt1 = 1.0_dp / h
        du0 = (h - 3.0_dp * d0**2 / h) / 6.0_dp
        du1 = (3.0_dp * d1**2 / h - h) / 6.0_dp
    end subroutine hermite_basis

    !> \brief 二维双三次样条求值 (自然边界; 网格外钳制到边界胞元外推)
    pure function spline_2d_eval(spl, xv, yv) result(val)
        type(spline_2d_t), intent(in) :: spl
        real(dp), intent(in) :: xv, yv
        real(dp) :: val

        integer :: i, j
        real(dp) :: t0x, t1x, u0x, u1x, dt0x, dt1x, du0x, du1x
        real(dp) :: t0y, t1y, u0y, u1y, dt0y, dt1y, du0y, du1y
        real(dp) :: z00, z10, z01, z11
        real(dp) :: zxx00, zxx10, zxx01, zxx11
        real(dp) :: zyy00, zyy10, zyy01, zyy11
        real(dp) :: zxy00, zxy10, zxy01, zxy11

        i = locate_cell(spl%x, spl%nx, xv)
        j = locate_cell(spl%y, spl%ny, yv)
        call hermite_basis(spl%x(i), spl%x(i + 1), xv, t0x, t1x, u0x, u1x, &
                           dt0x, dt1x, du0x, du1x)
        call hermite_basis(spl%y(j), spl%y(j + 1), yv, t0y, t1y, u0y, u1y, &
                           dt0y, dt1y, du0y, du1y)

        z00 = spl%z(i, j); z10 = spl%z(i + 1, j)
        z01 = spl%z(i, j + 1); z11 = spl%z(i + 1, j + 1)
        zxx00 = spl%zxx(i, j); zxx10 = spl%zxx(i + 1, j)
        zxx01 = spl%zxx(i, j + 1); zxx11 = spl%zxx(i + 1, j + 1)
        zyy00 = spl%zyy(i, j); zyy10 = spl%zyy(i + 1, j)
        zyy01 = spl%zyy(i, j + 1); zyy11 = spl%zyy(i + 1, j + 1)
        zxy00 = spl%zxxyy(i, j); zxy10 = spl%zxxyy(i + 1, j)
        zxy01 = spl%zxxyy(i, j + 1); zxy11 = spl%zxxyy(i + 1, j + 1)

        val = z00 * t0x * t0y + z10 * t1x * t0y + z01 * t0x * t1y + z11 * t1x * t1y &
              + zxx00 * u0x * t0y + zxx10 * u1x * t0y + zxx01 * u0x * t1y + zxx11 * u1x * t1y &
              + zyy00 * t0x * u0y + zyy10 * t1x * u0y + zyy01 * t0x * u1y + zyy11 * t1x * u1y &
              + zxy00 * u0x * u0y + zxy10 * u1x * u0y + zxy01 * u0x * u1y + zxy11 * u1x * u1y
    end function spline_2d_eval

    !> \brief 二维双三次样条解析梯度 (dS/dx, dS/dy)
    pure subroutine spline_2d_grad(spl, xv, yv, dzdx, dzdy)
        type(spline_2d_t), intent(in) :: spl
        real(dp), intent(in) :: xv, yv
        real(dp), intent(out) :: dzdx, dzdy

        integer :: i, j
        real(dp) :: t0x, t1x, u0x, u1x, dt0x, dt1x, du0x, du1x
        real(dp) :: t0y, t1y, u0y, u1y, dt0y, dt1y, du0y, du1y

        i = locate_cell(spl%x, spl%nx, xv)
        j = locate_cell(spl%y, spl%ny, yv)
        call hermite_basis(spl%x(i), spl%x(i + 1), xv, t0x, t1x, u0x, u1x, &
                           dt0x, dt1x, du0x, du1x)
        call hermite_basis(spl%y(j), spl%y(j + 1), yv, t0y, t1y, u0y, u1y, &
                           dt0y, dt1y, du0y, du1y)

        dzdx = spl%z(i, j) * dt0x * t0y + spl%z(i + 1, j) * dt1x * t0y &
               + spl%z(i, j + 1) * dt0x * t1y + spl%z(i + 1, j + 1) * dt1x * t1y &
               + spl%zxx(i, j) * du0x * t0y + spl%zxx(i + 1, j) * du1x * t0y &
               + spl%zxx(i, j + 1) * du0x * t1y + spl%zxx(i + 1, j + 1) * du1x * t1y &
               + spl%zyy(i, j) * dt0x * u0y + spl%zyy(i + 1, j) * dt1x * u0y &
               + spl%zyy(i, j + 1) * dt0x * u1y + spl%zyy(i + 1, j + 1) * dt1x * u1y &
               + spl%zxxyy(i, j) * du0x * u0y + spl%zxxyy(i + 1, j) * du1x * u0y &
               + spl%zxxyy(i, j + 1) * du0x * u1y + spl%zxxyy(i + 1, j + 1) * du1x * u1y

        dzdy = spl%z(i, j) * t0x * dt0y + spl%z(i + 1, j) * t1x * dt0y &
               + spl%z(i, j + 1) * t0x * dt1y + spl%z(i + 1, j + 1) * t1x * dt1y &
               + spl%zxx(i, j) * u0x * dt0y + spl%zxx(i + 1, j) * u1x * dt0y &
               + spl%zxx(i, j + 1) * u0x * dt1y + spl%zxx(i + 1, j + 1) * u1x * dt1y &
               + spl%zyy(i, j) * t0x * du0y + spl%zyy(i + 1, j) * t1x * du0y &
               + spl%zyy(i, j + 1) * t0x * du1y + spl%zyy(i + 1, j + 1) * t1x * du1y &
               + spl%zxxyy(i, j) * u0x * du0y + spl%zxxyy(i + 1, j) * u1x * du0y &
               + spl%zxxyy(i, j + 1) * u0x * du1y + spl%zxxyy(i + 1, j + 1) * u1x * du1y
    end subroutine spline_2d_grad

    !> \brief 释放二维样条内存
    subroutine spline_2d_clean(spl)
        type(spline_2d_t), intent(inout) :: spl
        if (allocated(spl%x)) deallocate (spl%x)
        if (allocated(spl%y)) deallocate (spl%y)
        if (allocated(spl%z)) deallocate (spl%z)
        if (allocated(spl%zxx)) deallocate (spl%zxx)
        if (allocated(spl%zyy)) deallocate (spl%zyy)
        if (allocated(spl%zxxyy)) deallocate (spl%zxxyy)
        spl%nx = 0
        spl%ny = 0
    end subroutine spline_2d_clean

end module mod_interpolation
