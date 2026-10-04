!> \brief Smolyak 稀疏网格求积模块 (Sparse Grid Quadrature)
!> \details 基于书中第18.4节与第18.8节规范实现高维稀疏网格求积：
!>          - 嵌套 Clenshaw-Curtis 一维积分规则 (点数 m = 2^(l-1) + 1)
!>          - 二维 Smolyak 稀疏组合网格构造 U_N^D f
!>          - 稀疏积分权重与节点生成
!> \author LiHao
module mod_sparse_grid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_constants, only: PI
    implicit none
    private

    public :: clenshaw_curtis_1d
    public :: smolyak_grid_2d_t
    public :: build_smolyak_grid_2d
    public :: destroy_smolyak_grid_2d
    public :: smolyak_integrate_2d

    !> \brief 二维 Smolyak 稀疏网格数据类型
    type :: smolyak_grid_2d_t
        integer :: n_points = 0               !< 稀疏网格总点数
        integer :: level = 0                  !< Smolyak 截断级别 N
        real(dp), allocatable :: pts(:, :)    !< 坐标矩阵 (2, n_points)，[-1, 1]^2
        real(dp), allocatable :: wts(:)       !< 积分权重向量 (n_points)
    end type smolyak_grid_2d_t

contains

    !> \brief 计算一维嵌套 Clenshaw-Curtis 求积规则点与权重
    !> \param[in] n 点数 (对于级别 l>=2, n = 2^(l-1) + 1; 对于 l=1, n = 1)
    subroutine clenshaw_curtis_1d(n, x, w)
        integer, intent(in)   :: n
        real(dp), intent(out) :: x(n), w(n)

        integer  :: i, k
        real(dp) :: theta, sum_term, b_coeff

        if (n <= 1) then
            x(1) = 0.0_dp
            w(1) = 2.0_dp
            return
        end if

        ! 节点分布: x_k = cos(pi * (k - 1) / (n - 1))
        do i = 1, n
            theta = PI * real(i - 1, dp) / real(n - 1, dp)
            x(i) = cos(theta)
        end do

        ! 权重解析式 (Waldvogel 紧凑算法/标准求积匹配)
        do i = 1, n
            theta = PI * real(i - 1, dp) / real(n - 1, dp)
            sum_term = 0.0_dp
            do k = 1, (n - 1) / 2
                b_coeff = 2.0_dp
                if (2 * k == n - 1) b_coeff = 1.0_dp
                sum_term = sum_term + b_coeff * cos(2.0_dp * real(k, dp) * theta) / (4.0_dp * real(k, dp)**2 - 1.0_dp)
            end do

            b_coeff = 1.0_dp
            if (i == 1 .or. i == n) b_coeff = 0.5_dp
            w(i) = (b_coeff * 2.0_dp / real(n - 1, dp)) * (1.0_dp - sum_term)
        end do
    end subroutine clenshaw_curtis_1d

    !> \brief 返回级别 level 对应的 Clenshaw-Curtis 点数
    pure integer function cc_num_points(level)
        integer, intent(in) :: level
        if (level <= 1) then
            cc_num_points = 1
        else
            cc_num_points = 2**(level - 1) + 1
        end if
    end function cc_num_points

    !> \brief 构造二维 Smolyak 稀疏网格
    !> \details 公式: U_N^2 = sum_{N-1 <= l1+l2 <= N} (-1)^{N - (l1+l2)} C(1, N - (l1+l2)) (U_{l1} x U_{l2})
    subroutine build_smolyak_grid_2d(n_level, grid)
        integer, intent(in)                 :: n_level
        type(smolyak_grid_2d_t), intent(out):: grid

        integer :: l1, l2, m1, m2, i1, i2, n_comb
        real(dp) :: c_comb, px, py, pw
        real(dp), allocatable :: x1(:), w1(:), x2(:), w2(:)
        real(dp), allocatable :: raw_pts(:, :), raw_wts(:)
        integer  :: total_alloc, cur_pt, p_idx
        logical  :: found

        grid%level = n_level

        ! 预统计网格容量
        total_alloc = 0
        do l1 = 1, n_level
            do l2 = 1, n_level
                if (l1 + l2 >= n_level - 1 .and. l1 + l2 <= n_level) then
                    total_alloc = total_alloc + cc_num_points(l1) * cc_num_points(l2)
                end if
            end do
        end do

        allocate(raw_pts(2, total_alloc), raw_wts(total_alloc))
        cur_pt = 0

        do l1 = 1, n_level
            m1 = cc_num_points(l1)
            allocate(x1(m1), w1(m1))
            call clenshaw_curtis_1d(m1, x1, w1)

            do l2 = 1, n_level
                if (l1 + l2 >= n_level - 1 .and. l1 + l2 <= n_level) then
                    m2 = cc_num_points(l2)
                    allocate(x2(m2), w2(m2))
                    call clenshaw_curtis_1d(m2, x2, w2)

                    ! Smolyak 组合系数: (-1)^(N - |l|) * C(D-1, N - |l|), 其中 D=2, C(1, 0)=1, C(1, 1)=1
                    n_comb = n_level - (l1 + l2)
                    if (n_comb == 0) then
                        c_comb = 1.0_dp
                    else
                        c_comb = -1.0_dp
                    end if

                    do i2 = 1, m2
                        do i1 = 1, m1
                            px = x1(i1)
                            py = x2(i2)
                            pw = c_comb * w1(i1) * w2(i2)

                            ! 合并重复节点
                            found = .false.
                            do p_idx = 1, cur_pt
                                if (abs(raw_pts(1, p_idx) - px) < 1.0e-12_dp .and. &
                                    abs(raw_pts(2, p_idx) - py) < 1.0e-12_dp) then
                                    raw_wts(p_idx) = raw_wts(p_idx) + pw
                                    found = .true.
                                    exit
                                end if
                            end do

                            if (.not. found) then
                                cur_pt = cur_pt + 1
                                raw_pts(1, cur_pt) = px
                                raw_pts(2, cur_pt) = py
                                raw_wts(cur_pt)    = pw
                            end if
                        end do
                    end do
                    deallocate(x2, w2)
                end if
            end do
            deallocate(x1, w1)
        end do

        grid%n_points = cur_pt
        allocate(grid%pts(2, cur_pt), grid%wts(cur_pt))
        grid%pts = raw_pts(:, 1:cur_pt)
        grid%wts = raw_wts(1:cur_pt)

        deallocate(raw_pts, raw_wts)
    end subroutine build_smolyak_grid_2d

    !> \brief 销毁 Smolyak 稀疏网格
    subroutine destroy_smolyak_grid_2d(grid)
        type(smolyak_grid_2d_t), intent(inout) :: grid
        if (allocated(grid%pts)) deallocate(grid%pts)
        if (allocated(grid%wts)) deallocate(grid%wts)
        grid%n_points = 0
        grid%level = 0
    end subroutine destroy_smolyak_grid_2d

    !> \brief 执行二维稀疏网格数值积分
    pure function smolyak_integrate_2d(grid, f_vals) result(integral_val)
        type(smolyak_grid_2d_t), intent(in) :: grid
        real(dp), intent(in)                :: f_vals(grid%n_points)
        real(dp) :: integral_val
        integral_val = sum(grid%wts * f_vals)
    end function smolyak_integrate_2d

end module mod_sparse_grid
