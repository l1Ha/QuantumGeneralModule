!> \brief 张量列车与矩阵乘积态单元测试 (Tensor Train Unit Tests)
program test_tensor_train
    use mod_constants, only: dp
    use mod_tensor_train
    implicit none

    type(tt_tensor_3d_t) :: tt_a, tt_b
    integer, parameter :: n1 = 4, n2 = 5, n3 = 3
    integer, parameter :: r1 = 2, r2 = 2
    real(dp) :: dense_a(n1, n2, n3), dense_b(n1, n2, n3)
    real(dp) :: dot_tt, dot_dense, norm_tt, norm_dense
    real(dp) :: diff_dot, diff_norm
    integer  :: i1, i2, i3, a1, a2
    integer  :: n_tests, n_pass

    n_tests = 0
    n_pass  = 0

    print '(A)', "=========================================================="
    print '(A)', "       GeneralModule Unit Tests: Tensor Train (TT/MPS)    "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: TT 核初始化与结构检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call init_tt_tensor_3d(tt_a, n1, n2, n3, r1, r2)
    call init_tt_tensor_3d(tt_b, n1, n2, n3, r1, r2)

    if (tt_a%n1 == n1 .and. tt_a%r1 == r1 .and. tt_a%r2 == r2) then
        print '(A)', " [PASS] 3D Tensor Train cores successfully allocated"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] TT allocation failed"
    end if

    ! 填充确定数值
    do a1 = 1, r1; do i1 = 1, n1
        tt_a%g1(1, i1, a1) = sin(real(i1 + 2 * a1, dp))
        tt_b%g1(1, i1, a1) = cos(real(2 * i1 - a1, dp))
    end do; end do

    do a2 = 1, r2; do i2 = 1, n2; do a1 = 1, r1
        tt_a%g2(a1, i2, a2) = cos(real(a1 * i2 + a2, dp))
        tt_b%g2(a1, i2, a2) = sin(real(2 * a1 + i2 - a2, dp))
    end do; end do; end do

    do i3 = 1, n3; do a2 = 1, r2
        tt_a%g3(a2, i3, 1) = sin(real(3 * a2 + i3, dp))
        tt_b%g3(a2, i3, 1) = cos(real(a2 - 2 * i3, dp))
    end do; end do

    ! --------------------------------------------------------------------------
    ! Test 2: 稠密张量还原检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call tt_3d_to_dense(tt_a, dense_a)
    call tt_3d_to_dense(tt_b, dense_b)

    if (maxval(abs(dense_a)) > 0.0_dp .and. maxval(abs(dense_b)) > 0.0_dp) then
        print '(A)', " [PASS] TT cores contracted into dense tensor"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Dense tensor reconstruction is zero"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: 快速 TT 内积 vs 稠密全张量点乘
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    dot_tt = tt_3d_inner_product(tt_a, tt_b)
    dot_dense = sum(dense_a * dense_b)
    diff_dot = abs(dot_tt - dot_dense) / max(abs(dot_dense), 1.0e-30_dp)

    if (diff_dot < 1.0e-12_dp) then
        print '(A, ES10.3)', " [PASS] Fast TT inner product matches dense Frobenius dot: diff = ", diff_dot
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] TT inner product mismatch: diff = ", diff_dot
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: TT 范数与全张量 Frobenius 范数一致性
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    norm_tt = tt_3d_norm(tt_a)
    norm_dense = sqrt(sum(dense_a**2))
    diff_norm = abs(norm_tt - norm_dense) / max(norm_dense, 1.0e-30_dp)

    if (diff_norm < 1.0e-12_dp) then
        print '(A, ES10.3)', " [PASS] TT Frobenius norm matches dense tensor norm: diff = ", diff_norm
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] TT norm mismatch: diff = ", diff_norm
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: 释放资源后状态检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call destroy_tt_tensor_3d(tt_a)
    call destroy_tt_tensor_3d(tt_b)
    if (tt_a%n1 == 0 .and. .not. allocated(tt_a%g1)) then
        print '(A)', " [PASS] TT tensors successfully deallocated"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] TT deallocation failed"
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  Tensor Train Tests: ", n_pass, " / ", n_tests, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_tests) error stop "Test failures detected in test_tensor_train."
end program test_tensor_train
