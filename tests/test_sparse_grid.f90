!> \brief Smolyak 稀疏网格求积单元测试 (Sparse Grid Quadrature Unit Tests)
program test_sparse_grid
    use mod_constants, only: dp
    use mod_sparse_grid
    implicit none

    integer, parameter :: m = 9
    real(dp) :: x1d(m), w1d(m)
    type(smolyak_grid_2d_t) :: grid
    real(dp), allocatable :: f_vals(:)
    real(dp) :: sum_w1d, sum_w2d, int_calc, int_exact, err_int
    integer  :: i, n_tests, n_pass

    n_tests = 0
    n_pass  = 0

    print '(A)', "=========================================================="
    print '(A)', "       GeneralModule Unit Tests: Sparse Grid (Smolyak)    "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: 一维 Clenshaw-Curtis 权重和检验 (常数 1 在 [-1, 1] 积分恒为 2)
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call clenshaw_curtis_1d(m, x1d, w1d)
    sum_w1d = sum(w1d)

    if (abs(sum_w1d - 2.0_dp) < 1.0e-12_dp) then
        print '(A, F8.5)', " [PASS] 1D Clenshaw-Curtis weights sum to 2.0: ", sum_w1d
        n_pass = n_pass + 1
    else
        print '(A, F8.5)', " [FAIL] 1D CC weights sum mismatch: ", sum_w1d
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: 一维 Clenshaw-Curtis 单项式积分精度检验 (int x^2 dx = 2/3)
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    int_calc = sum(w1d * (x1d**2))
    int_exact = 2.0_dp / 3.0_dp

    if (abs(int_calc - int_exact) < 1.0e-12_dp) then
        print '(A, ES10.3)', " [PASS] 1D CC monomial x^2 quadrature error: ", abs(int_calc - int_exact)
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] 1D CC monomial error: ", abs(int_calc - int_exact)
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: 二维 Smolyak 稀疏网格节点坐标范围与生成检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call build_smolyak_grid_2d(n_level=4, grid=grid)

    if (grid%n_points > 0 .and. maxval(abs(grid%pts)) <= 1.0_dp + 1.0e-12_dp) then
        print '(A, I4, A)', " [PASS] 2D Smolyak grid level 4 built with ", grid%n_points, " sparse nodes"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] 2D Smolyak grid nodes out of range [-1, 1]^2"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: 二维 Smolyak 稀疏权重和检验 (常数 1 在 [-1, 1]^2 积分恒等于 4)
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    sum_w2d = sum(grid%wts)

    if (abs(sum_w2d - 4.0_dp) < 1.0e-12_dp) then
        print '(A, F8.5)', " [PASS] 2D Smolyak weights sum to 4.0: ", sum_w2d
        n_pass = n_pass + 1
    else
        print '(A, F8.5)', " [FAIL] 2D Smolyak weights sum mismatch: ", sum_w2d
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: 二维可分离光滑函数积分检验: f(x, y) = cos(x) * cos(y)
    !         解析值 int_{-1}^1 int_{-1}^1 cos(x)cos(y) dx dy = 4 * sin(1)^2
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    allocate(f_vals(grid%n_points))
    do i = 1, grid%n_points
        f_vals(i) = cos(grid%pts(1, i)) * cos(grid%pts(2, i))
    end do

    int_calc = smolyak_integrate_2d(grid, f_vals)
    int_exact = 4.0_dp * (sin(1.0_dp)**2)
    err_int = abs(int_calc - int_exact)
    deallocate(f_vals)

    if (err_int < 1.0e-2_dp) then
        print '(A, ES10.3)', " [PASS] 2D Smolyak integration of cos(x)*cos(y) matches analytical: err = ", err_int
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] 2D Smolyak integration error too large: err = ", err_int
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: 释放网格资源检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call destroy_smolyak_grid_2d(grid)
    if (grid%n_points == 0 .and. .not. allocated(grid%pts)) then
        print '(A)', " [PASS] 2D Smolyak grid successfully deallocated"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] 2D Smolyak grid deallocation failed"
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  Sparse Grid Tests: ", n_pass, " / ", n_tests, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_tests) error stop "Test failures detected in test_sparse_grid."
end program test_sparse_grid
