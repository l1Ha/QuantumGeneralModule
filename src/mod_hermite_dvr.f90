!> \brief Hermite DVR 模块 (谐振子基; 路线图 §12.1「更多正交多项式 DVR」项)
!> \details 物理学家 Hermite 多项式 H_n 与归一化 Hermite 函数
!>          psi_n(x) = (2^n n! sqrt(pi))^{-1/2} e^{-x^2/2} H_n(x)
!>          构成谐振子 H = p^2/(2m) + m omega^2 x^2 / 2 (hbar = 1) 的精确本征基。
!>          质量-频率标度坐标 q = sqrt(m omega) x 下基函数不含 m/omega,
!>          动能矩阵元由梯位算符严格给出 (无近似):
!>            <i|p_q^2|j> = (2j+1) d_{ij} - sqrt(j(j-1)) d_{i,j-2}
!>                          - sqrt((j+1)(j+2)) d_{i,j+2},
!>            T_ij = (omega/2) <i|p_q^2|j>  (p_phys = sqrt(m omega) p_q)。
!>          势能矩阵元经 Gauss-Hermite 求值 (Golub-Welsch 节点/权重):
!>            V_ij = sum_k w_k psit_i(q_k) V(q_k/sqrt(m omega)) psit_j(q_k),
!>          其中 psit_n(q) = psi_n(q) e^{+q^2/2} = H_n(q)/sqrt(2^n n! sqrt(pi))
!>          为去指数多项式部分 —— 权重 w_k 已含 e^{-q^2}, 若直接用含
!>          e^{-q^2/2} 的 psi 值则等价于以 e^{-2q^2} 为权, 全部矩阵元错误。
!>          该求值对 2N-1 次以下多项式势能精确 (V 以多项式乘入 e^{-q^2} 权)。
!> \author LiHao
module mod_hermite_dvr
    use mod_constants, only: dp, PI
    implicit none
    private

    public :: gauss_hermite_nodes_weights
    public :: hermite_function_values
    public :: hermite_ho_kinetic
    public :: hermite_dvr_hamiltonian

contains

    !> \brief 物理学家 Hermite 多项式 H_n(x) 与导数 H_n'(x) (0..n_max)
    pure subroutine hermite_poly_deriv(x, n_max, hn, dhn)
        real(dp), intent(in) :: x
        integer, intent(in) :: n_max
        real(dp), intent(out) :: hn(0:n_max), dhn(0:n_max)
        integer :: k
        hn(0) = 1.0_dp
        dhn(0) = 0.0_dp
        if (n_max >= 1) then
            hn(1) = 2.0_dp * x
            dhn(1) = 2.0_dp
        end if
        do k = 1, n_max - 1
            hn(k + 1) = 2.0_dp * x * hn(k) - 2.0_dp * real(k, dp) * hn(k - 1)
            dhn(k + 1) = 2.0_dp * real(k + 1, dp) * hn(k)
        end do
    end subroutine hermite_poly_deriv

    !> \brief Gauss-Hermite 节点与权重: ∫ e^{-x^2} f(x) dx ≈ Σ w_k f(x_k)
    !> \details **Golub–Welsch 算法** (稳健, 无求根): 单项化 Hermite 递推
    !>          p_{k+1} = x p_k - (k/2) p_{k-1} 给出 Jacobi 三对角矩阵
    !>          (对角 0, 次对角 sqrt(k/2)), 其本征值即节点, 本征矢量首分量
    !>          平方 × sqrt(pi) 即权重。早期 Newton 求根版本会跳根
    !>          (初值不当), 已由本实现替换。
    subroutine gauss_hermite_nodes_weights(n_ord, x_nodes, w_gauss)
        use mod_linear_algebra, only: diag_symmetric_matrix
        integer, intent(in) :: n_ord
        real(dp), intent(out) :: x_nodes(n_ord), w_gauss(n_ord)

        real(dp), allocatable :: jac(:, :), w_vec(:, :)
        real(dp) :: mu0
        integer :: k, st

        allocate (jac(n_ord, n_ord), w_vec(n_ord, n_ord))
        jac = 0.0_dp
        do k = 1, n_ord - 1
            jac(k, k + 1) = sqrt(0.5_dp * real(k, dp))
            jac(k + 1, k) = jac(k, k + 1)
        end do
        call diag_symmetric_matrix(n_ord, jac, x_nodes, w_vec, st)
        ! 权重: w_k = sqrt(pi) * v_{1k}^2  (Golub-Welsch, mu_0 = sqrt(pi))
        mu0 = sqrt(PI)
        do k = 1, n_ord
            w_gauss(k) = mu0 * w_vec(1, k)**2
        end do
        deallocate (jac, w_vec)
    end subroutine gauss_hermite_nodes_weights

    !> \brief 归一化 Hermite 函数 psi_n(x) = c_n e^{-x^2/2} H_n(x)
    pure subroutine hermite_function_values(x, n_max, psi)
        real(dp), intent(in) :: x
        integer, intent(in) :: n_max
        real(dp), intent(out) :: psi(0:n_max)
        real(dp) :: hn(0:n_max), dhn(0:n_max)
        integer :: k
        call hermite_poly_deriv(x, n_max, hn, dhn)
        ! 归一化因子在 log 域合成: ln norm = 0.5 (k ln2 + lgamma(k+1) + 0.5 ln pi)。
        ! 直接计算 2^k k! 会在 k >~ 165 溢出 (2^159 159! ~ 1e329),
        ! 合并进单个 exp 后连同 e^{-x^2/2} 一起求值, 同时少一次舍入。
        do k = 0, n_max
            psi(k) = hn(k) * exp(-0.5_dp * (real(k, dp) * log(2.0_dp) + log_gamma(real(k + 1, dp)) &
                                            + 0.5_dp * log(PI)) - 0.5_dp * x**2)
        end do
    end subroutine hermite_function_values

    !> \brief HO 动能矩阵 (梯位闭式): T_ij = (omega/2) <i|p_q^2|j>
    subroutine hermite_ho_kinetic(n_ord, omega, t_mat)
        integer, intent(in) :: n_ord
        real(dp), intent(in) :: omega
        real(dp), intent(out) :: t_mat(n_ord, n_ord)

        integer :: i, j
        real(dp) :: qj
        ! p_q = i(a^dag - a)/sqrt(2) => <j|p_q^2|j> = j + 1/2 (0-based j);
        !       <j-2|p_q^2|j> = -(1/2) sqrt(j(j-1)); <j+2|...> = -(1/2) sqrt(...)
        ! T_ij = (omega/2) <i|p_q^2|j>
        t_mat = 0.0_dp
        do j = 1, n_ord
            qj = real(j - 1, dp)
            do i = 1, n_ord
                if (i == j) t_mat(i, j) = 0.5_dp * omega * (qj + 0.5_dp)
                if (i == j - 2 .and. j >= 3) then
                    t_mat(i, j) = -0.25_dp * omega * sqrt(qj * (qj - 1.0_dp))
                    t_mat(j, i) = t_mat(i, j)
                end if
                if (i == j + 2) then
                    t_mat(i, j) = -0.25_dp * omega * sqrt((qj + 1.0_dp) * (qj + 2.0_dp))
                    t_mat(j, i) = t_mat(i, j)
                end if
            end do
        end do
    end subroutine hermite_ho_kinetic

    !> \brief Hermite DVR 哈密顿量组装 (势能由调用方在节点上提供)
    !> \param[in] n_ord 基维数
    !> \param[in] omega 频率 (标度坐标 q = sqrt(m omega) x; m 由调用方折算进 V)
    !> \param[in] v_at_nodes 势能在 q 空间节点上的值 (n_ord): 调用方先取节点
    !>            (gauss_hermite_nodes_weights) 再以 V(x = q/sqrt(m omega)) 求值
    !> \param[out] h_mat 哈密顿矩阵 (n_ord, n_ord)
    !> \param[out] x_nodes q 空间 Gauss-Hermite 节点
    subroutine hermite_dvr_hamiltonian(n_ord, omega, v_at_nodes, h_mat, x_nodes)
        integer, intent(in) :: n_ord
        real(dp), intent(in) :: omega
        real(dp), intent(in) :: v_at_nodes(n_ord)
        real(dp), intent(out) :: h_mat(n_ord, n_ord)
        real(dp), intent(out) :: x_nodes(n_ord)

        real(dp) :: w_gauss(n_ord), psi_k(0:n_ord - 1), t_mat(n_ord, n_ord)
        integer :: i, j, k

        call gauss_hermite_nodes_weights(n_ord, x_nodes, w_gauss)
        h_mat = 0.0_dp
        do k = 1, n_ord
            call hermite_function_values(x_nodes(k), n_ord - 1, psi_k)
            ! 关键: 权重 w_k 已含 e^{-q^2}; 函数值须去掉自身的 e^{-q^2/2}
            ! (乘 e^{+q^2/2}), 使被积函数恰为 权重 x 多项式 x V —— 否则
            ! 相当于用 e^{-2q^2} 权重求值, 全部矩阵元系统性错误。
            psi_k = psi_k * exp(0.5_dp * x_nodes(k)**2)
            do i = 1, n_ord
                do j = 1, n_ord
                    h_mat(i, j) = h_mat(i, j) + w_gauss(k) * psi_k(i - 1) &
                                  * v_at_nodes(k) * psi_k(j - 1)
                end do
            end do
        end do
        call hermite_ho_kinetic(n_ord, omega, t_mat)
        h_mat = h_mat + t_mat
    end subroutine hermite_dvr_hamiltonian

end module mod_hermite_dvr
