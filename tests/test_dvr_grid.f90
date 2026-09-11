!> \brief DVR 网格与 FGH 束缚态求解单元测试
program test_dvr_grid
    use general_module
    implicit none

    integer :: n_tests = 0
    integer :: n_passed = 0
    real(dp) :: tol = 1.0e-4_dp

    type(dvr_1d_t) :: dvr_sinc
    type(dvr_legendre_t) :: dvr_leg
    integer :: n_pts, stat, i, j
    real(dp) :: x_min, x_max, mass, omega, max_asym, sum_weights
    real(dp), allocatable :: v_pot(:), eig_vals(:), eig_vecs(:, :)

    print '(A)', "=================================================="
    print '(A)', "       GeneralModule Unit Tests: DVR Grid         "
    print '(A)', "=================================================="

    ! 1. Sinc-DVR 动能矩阵对称性校验
    x_min = -10.0_dp
    x_max = 10.0_dp
    n_pts = 64
    mass = 1.0_dp
    call dvr_sinc_init(x_min, x_max, n_pts, mass, dvr_sinc)

    max_asym = 0.0_dp
    do i = 1, n_pts
        do j = 1, n_pts
            max_asym = max(max_asym, abs(dvr_sinc%t_mat(i, j) - dvr_sinc%t_mat(j, i)))
        end do
    end do
    call assert_close("Sinc-DVR T_ij symmetry", 0.0_dp, max_asym, 1.0e-12_dp)

    ! 2. 谐振子 FGH 束缚态本征能级求解
    ! V(x) = 1/2 * m * omega^2 * x^2, m=1, omega=1 -> E_v = v + 0.5
    omega = 1.0_dp
    allocate(v_pot(n_pts))
    allocate(eig_vals(n_pts))
    allocate(eig_vecs(n_pts, n_pts))

    do i = 1, n_pts
        v_pot(i) = 0.5_dp * mass * (omega**2) * (dvr_sinc%x(i)**2)
    end do

    call fgh_solve_bound_states(dvr_sinc, v_pot, eig_vals, eig_vecs, stat)
    call assert_close("FGH status code", 0.0_dp, real(stat, dp), 1.0e-12_dp)

    ! 检验前 4 个本征态能级
    call assert_close("HO Ground State E0 = 0.5", 0.5_dp, eig_vals(1), tol)
    call assert_close("HO First Excited E1 = 1.5", 1.5_dp, eig_vals(2), tol)
    call assert_close("HO Second Excited E2 = 2.5", 2.5_dp, eig_vals(3), tol)
    call assert_close("HO Third Excited E3 = 3.5", 3.5_dp, eig_vals(4), tol)

    ! 3. Gauss-Legendre DVR 权重归一化: sum(w_i) = 2.0 (在 [-1, 1] 积分 1*dx = 2)
    call dvr_legendre_init(32, dvr_leg)
    sum_weights = sum(dvr_leg%weights)
    call assert_close("Legendre weights sum = 2.0", 2.0_dp, sum_weights, 1.0e-10_dp)

    ! 节点反对称性 x_i = -x_{N+1-i}
    max_asym = 0.0_dp
    do i = 1, 32
        max_asym = max(max_asym, abs(dvr_leg%x(i) + dvr_leg%x(33 - i)))
    end do
    call assert_close("Legendre nodes antisymmetry", 0.0_dp, max_asym, 1.0e-12_dp)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "DVR Grid Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some tests did not pass!"
        stop 1
    else
        print '(A)', "SUCCESS: All DVR grid tests passed."
    end if

contains

    subroutine assert_close(name, expected, actual, eps)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, eps
        real(dp) :: diff

        n_tests = n_tests + 1
        diff = abs(expected - actual)
        if (diff <= eps) then
            n_passed = n_passed + 1
            print '(A, A35, A)', " [PASS] ", name, ""
        else
            print '(A, A35, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_dvr_grid
