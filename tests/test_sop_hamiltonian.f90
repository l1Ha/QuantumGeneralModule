!> \brief 和积哈密顿量单元测试 (Sum-of-Products Hamiltonian Unit Tests)
program test_sop_hamiltonian
    use mod_constants, only: dp
    use mod_sop_hamiltonian
    implicit none

    type(sop_hamiltonian_3d_t) :: sop
    integer, parameter :: n1 = 3, n2 = 4, n3 = 2
    integer, parameter :: n_tot = n1 * n2 * n3
    real(dp) :: in_vec(n1, n2, n3), out_vec(n1, n2, n3)
    real(dp) :: in_vec2(n1, n2, n3), out_vec2(n1, n2, n3), out_linear(n1, n2, n3)
    real(dp) :: h_kron(n_tot, n_tot), dense_ref(n_tot)
    real(dp) :: max_diff, rel_err
    integer  :: i1, i2, i3, j1, j2, j3, row, col, m
    integer  :: n_tests, n_pass

    n_tests = 0
    n_pass  = 0

    print '(A)', "=========================================================="
    print '(A)', "       GeneralModule Unit Tests: SOP Hamiltonian          "
    print '(A)', "=========================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: 初始化与结构体分配检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call init_sop_hamiltonian_3d(sop, n1, n2, n3, n_terms=2)

    if (sop%n1 == n1 .and. sop%n2 == n2 .and. sop%n3 == n3 .and. sop%n_terms == 2) then
        print '(A)', " [PASS] 3D SOP Hamiltonian successfully allocated"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] 3D SOP allocation failed"
    end if

    ! 构造具有确定数值的 2 项 SOP 算符
    ! Term 1: c1 = 1.5, h1_1, h2_1, h3_1
    sop%coeffs(1) = 1.5_dp
    do j1 = 1, n1; do i1 = 1, n1
        sop%h1(1)%mat(i1, j1) = real(i1 + 2 * j1, dp) * 0.1_dp
    end do; end do
    do j2 = 1, n2; do i2 = 1, n2
        sop%h2(1)%mat(i2, j2) = real(i2 - j2, dp) * 0.2_dp
    end do; end do
    do j3 = 1, n3; do i3 = 1, n3
        sop%h3(1)%mat(i3, j3) = real(2 * i3 + j3, dp) * 0.15_dp
    end do; end do

    ! Term 2: c2 = -0.8
    sop%coeffs(2) = -0.8_dp
    do j1 = 1, n1; do i1 = 1, n1
        sop%h1(2)%mat(i1, j1) = real(i1 * j1, dp) * 0.05_dp
    end do; end do
    do j2 = 1, n2; do i2 = 1, n2
        sop%h2(2)%mat(i2, j2) = real(i2 + j2, dp) * 0.12_dp
    end do; end do
    do j3 = 1, n3; do i3 = 1, n3
        sop%h3(2)%mat(i3, j3) = real(i3 - 2 * j3, dp) * 0.08_dp
    end do; end do

    ! --------------------------------------------------------------------------
    ! Test 2: 矩阵自由算子作用与全 Kronecker 积稠密矩阵对比
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    ! 随机初始化测试向量
    do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
        in_vec(i1, i2, i3) = sin(real(i1 + 3 * i2 + 5 * i3, dp))
    end do; end do; end do

    ! 1. 矩阵自由 SOP 算子作用
    call sop_apply_3d(sop, in_vec, out_vec)

    ! 2. 构造全 Kronecker 积稠密参考矩阵
    h_kron = 0.0_dp
    do m = 1, sop%n_terms
        do j3 = 1, n3; do j2 = 1, n2; do j1 = 1, n1
            col = j1 + (j2 - 1) * n1 + (j3 - 1) * n1 * n2
            do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
                row = i1 + (i2 - 1) * n1 + (i3 - 1) * n1 * n2
                h_kron(row, col) = h_kron(row, col) + sop%coeffs(m) * &
                    sop%h1(m)%mat(i1, j1) * sop%h2(m)%mat(i2, j2) * sop%h3(m)%mat(i3, j3)
            end do; end do; end do
        end do; end do; end do
    end do

    dense_ref = matmul(h_kron, reshape(in_vec, (/ n_tot /)))
    max_diff = maxval(abs(reshape(out_vec, (/ n_tot /)) - dense_ref))
    rel_err = max_diff / max(maxval(abs(dense_ref)), 1.0e-30_dp)

    if (rel_err < 1.0e-12_dp) then
        print '(A, ES10.3)', " [PASS] Matrix-free SOP action matches dense Kronecker product: rel_err = ", rel_err
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] SOP action mismatch: rel_err = ", rel_err
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: 算子作用线性度检验 H(a*u + b*v) == a*Hu + b*Hv
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    do i3 = 1, n3; do i2 = 1, n2; do i1 = 1, n1
        in_vec2(i1, i2, i3) = cos(real(2 * i1 - i2 + 4 * i3, dp))
    end do; end do; end do

    call sop_apply_3d(sop, in_vec2, out_vec2)
    call sop_apply_3d(sop, 2.0_dp * in_vec - 1.5_dp * in_vec2, out_linear)

    max_diff = maxval(abs(out_linear - (2.0_dp * out_vec - 1.5_dp * out_vec2)))
    if (max_diff < 1.0e-12_dp) then
        print '(A, ES10.3)', " [PASS] Operator linearity strictly satisfied: max_diff = ", max_diff
        n_pass = n_pass + 1
    else
        print '(A, ES10.3)', " [FAIL] Operator linearity failed: max_diff = ", max_diff
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: 释放资源后状态检验
    ! --------------------------------------------------------------------------
    n_tests = n_tests + 1
    call destroy_sop_hamiltonian_3d(sop)
    if (sop%n_terms == 0 .and. .not. allocated(sop%coeffs)) then
        print '(A)', " [PASS] 3D SOP Hamiltonian successfully deallocated"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] 3D SOP deallocation failed"
    end if

    print '(A)', "=========================================================="
    print '(A, I2, A, I2, A)', "  SOP Hamiltonian Tests: ", n_pass, " / ", n_tests, " passed."
    print '(A)', "=========================================================="

    if (n_pass /= n_tests) error stop "Test failures detected in test_sop_hamiltonian."
end program test_sop_hamiltonian
