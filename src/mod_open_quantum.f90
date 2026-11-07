!> \brief 开放量子系统与 Lindblad 耗散主方程动力学模块
!> \details 提供密度矩阵表示、自发辐射弛豫 (T_1)、纯退相位 (T_2^*) 及热浴耦合 Jump 算符构建、
!>          Runge-Kutta 密度矩阵含时主方程推进器，以及量子纯度、相干度与冯·诺依曼熵诊断例程；
!>          Lindblad Liouvillian 谱方法 (多原子拓展第九步 / 开放系统谱方法):
!>          向量化超算符、稳态 (增广线性方程 + 闭式)、阻尼基谱分解与谱重构。
!> \author LiHao
module mod_open_quantum
    use mod_constants, only: dp, EYE, KB
    use mod_linear_algebra, only: diag_symmetric_matrix, inv_complex_matrix
    implicit none
    private

    public :: lindblad_dissipator
    public :: rk4_lindblad_step
    public :: calculate_quantum_purity
    public :: calculate_von_neumann_entropy
    public :: calculate_quantum_coherence
    public :: create_relaxation_jump_op
    public :: create_dephasing_jump_op
    public :: lindblad_supermatrix_2level
    public :: lindblad_steady_state_2level
    public :: lindblad_eigenmodes_2level
    public :: lindblad_spectral_reconstruct

contains

    !> \brief 构造单态自发跃迁弛豫算符 L = |i><j| (从激发态 j 衰变到基态 i)
    subroutine create_relaxation_jump_op(n, lower_idx, upper_idx, l_op)
        integer, intent(in) :: n, lower_idx, upper_idx
        complex(dp), intent(out) :: l_op(n, n)

        l_op = (0.0_dp, 0.0_dp)
        if (lower_idx >= 1 .and. lower_idx <= n .and. upper_idx >= 1 .and. upper_idx <= n) then
            l_op(lower_idx, upper_idx) = (1.0_dp, 0.0_dp)
        end if
    end subroutine create_relaxation_jump_op

    !> \brief 构造局域纯退相位算符 L = |k><k|
    subroutine create_dephasing_jump_op(n, state_idx, l_op)
        integer, intent(in) :: n, state_idx
        complex(dp), intent(out) :: l_op(n, n)

        l_op = (0.0_dp, 0.0_dp)
        if (state_idx >= 1 .and. state_idx <= n) then
            l_op(state_idx, state_idx) = (1.0_dp, 0.0_dp)
        end if
    end subroutine create_dephasing_jump_op

    !> \brief 计算单个 Lindblad 耗散超算符作用: D[L] rho = L * rho * L^dagger - 0.5 * (L^dagger * L * rho + rho * L^dagger * L)
    pure function lindblad_dissipator(l_op, rho) result(d_rho)
        complex(dp), intent(in) :: l_op(:, :), rho(:, :)
        complex(dp) :: d_rho(size(rho, 1), size(rho, 2))

        integer :: n
        complex(dp) :: l_dag(size(rho, 1), size(rho, 2))
        complex(dp) :: l_dag_l(size(rho, 1), size(rho, 2))
        complex(dp) :: term1(size(rho, 1), size(rho, 2))
        complex(dp) :: term2(size(rho, 1), size(rho, 2))
        complex(dp) :: term3(size(rho, 1), size(rho, 2))

        n = size(rho, 1)
        l_dag = conjg(transpose(l_op))
        l_dag_l = matmul(l_dag, l_op)

        term1 = matmul(l_op, matmul(rho, l_dag))
        term2 = matmul(l_dag_l, rho)
        term3 = matmul(rho, l_dag_l)

        d_rho = term1 - 0.5_dp * (term2 + term3)
    end function lindblad_dissipator

    !> \brief 4 阶 Runge-Kutta 推进含时 Lindblad 密度矩阵主方程:
    !>        d rho / dt = -i [H, rho] + sum_k gamma_k D[L_k] rho
    subroutine rk4_lindblad_step(rho, h_mat, l_ops, gamma_rates, dt)
        complex(dp), intent(inout) :: rho(:, :)
        complex(dp), intent(in) :: h_mat(:, :)
        complex(dp), intent(in) :: l_ops(:, :, :)  !< (n, n, n_channels)
        real(dp), intent(in) :: gamma_rates(:)     !< (n_channels)
        real(dp), intent(in) :: dt

        integer :: n, n_ops
        complex(dp), allocatable :: k1(:, :), k2(:, :), k3(:, :), k4(:, :), r_tmp(:, :)

        n = size(rho, 1)
        n_ops = size(gamma_rates)

        allocate(k1(n, n), k2(n, n), k3(n, n), k4(n, n), r_tmp(n, n))

        k1 = deriv_rho(rho, h_mat, l_ops, gamma_rates, n, n_ops)
        r_tmp = rho + 0.5_dp * dt * k1
        k2 = deriv_rho(r_tmp, h_mat, l_ops, gamma_rates, n, n_ops)
        r_tmp = rho + 0.5_dp * dt * k2
        k3 = deriv_rho(r_tmp, h_mat, l_ops, gamma_rates, n, n_ops)
        r_tmp = rho + dt * k3
        k4 = deriv_rho(r_tmp, h_mat, l_ops, gamma_rates, n, n_ops)

        rho = rho + (dt / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)

        ! 保持迹归一化 Tr(rho) = 1 与厄米对称性
        rho = 0.5_dp * (rho + conjg(transpose(rho)))
        call enforce_trace_conservation(rho)

        deallocate(k1, k2, k3, k4, r_tmp)
    end subroutine rk4_lindblad_step

    !> \brief 计算密度矩阵时间导数 d rho / dt
    pure function deriv_rho(r_in, h_mat, l_ops, gamma_rates, n, n_ops) result(drdt)
        integer, intent(in) :: n, n_ops
        complex(dp), intent(in) :: r_in(n, n), h_mat(n, n), l_ops(n, n, n_ops)
        real(dp), intent(in) :: gamma_rates(n_ops)
        complex(dp) :: drdt(n, n)

        integer :: k
        complex(dp) :: comm(n, n)

        ! -i [H, rho] = -i (H * rho - rho * H)
        comm = matmul(h_mat, r_in) - matmul(r_in, h_mat)
        drdt = -EYE * comm

        ! Lindblad 耗散求和
        do k = 1, n_ops
            if (gamma_rates(k) > 0.0_dp) then
                drdt = drdt + gamma_rates(k) * lindblad_dissipator(l_ops(:, :, k), r_in)
            end if
        end do
    end function deriv_rho

    !> \brief 强制保证 Tr(rho) = 1.0
    pure subroutine enforce_trace_conservation(rho)
        complex(dp), intent(inout) :: rho(:, :)
        integer :: n, i
        real(dp) :: tr

        n = size(rho, 1)
        tr = 0.0_dp
        do i = 1, n
            tr = tr + real(rho(i, i), dp)
        end do
        if (tr > 1.0e-12_dp) then
            rho = rho / tr
        end if
    end subroutine enforce_trace_conservation

    !> \brief 计算量子态纯度 Purity = Tr(rho^2) ∈ [1/n, 1]
    pure function calculate_quantum_purity(rho) result(purity)
        complex(dp), intent(in) :: rho(:, :)
        real(dp) :: purity

        complex(dp) :: rho2(size(rho, 1), size(rho, 2))
        integer :: i, n

        n = size(rho, 1)
        rho2 = matmul(rho, rho)
        purity = 0.0_dp
        do i = 1, n
            purity = purity + real(rho2(i, i), dp)
        end do
    end function calculate_quantum_purity

    !> \brief 计算 l1 范数量子相干度 C_l1 = sum_{i /= j} |rho_ij|
    pure function calculate_quantum_coherence(rho) result(coherence)
        complex(dp), intent(in) :: rho(:, :)
        real(dp) :: coherence

        integer :: i, j, n

        n = size(rho, 1)
        coherence = 0.0_dp
        do i = 1, n
            do j = 1, n
                if (i /= j) then
                    coherence = coherence + abs(rho(i, j))
                end if
            end do
        end do
    end function calculate_quantum_coherence

    !> \brief 计算冯·诺依曼熵 S = - Tr(rho * ln(rho)) = - sum_i p_i * ln(p_i)
    subroutine calculate_von_neumann_entropy(rho, entropy)
        complex(dp), intent(in) :: rho(:, :)
        real(dp), intent(out) :: entropy

        integer :: n, stat, i
        real(dp), allocatable :: r_real(:, :), eig_d(:), eig_z(:, :)
        real(dp) :: p_val

        n = size(rho, 1)
        entropy = 0.0_dp

        ! 提取实对称部分求谱
        allocate(r_real(n, n), eig_d(n), eig_z(n, n))
        do i = 1, n
            r_real(i, :) = real(rho(i, :), dp)
        end do

        call diag_symmetric_matrix(n, r_real, eig_d, eig_z, stat)

        do i = 1, n
            p_val = eig_d(i)
            if (p_val > 1.0e-12_dp) then
                entropy = entropy - p_val * log(p_val)
            end if
        end do

        deallocate(r_real, eig_d, eig_z)
    end subroutine calculate_von_neumann_entropy

    !> \brief 二能级热化 Lindblad 的向量化 Liouvillian 超矩阵 (4x4 复)
    !> \details 基矢排序 (ee, eg, ge, gg); H = delta |e><e|, 塌缩速率
    !>          gamma_down = gamma (n_th+1) (向下), gamma_up = gamma n_th (向上),
    !>          纯退相位 gamma_phi (sigma_z 约定, 相干衰减率 = Gamma/2 + 2 gamma_phi
    !>          其中 Gamma = gamma_down + gamma_up):
    !>            d rho_ee = -gamma_down rho_ee + gamma_up rho_gg
    !>            d rho_eg = -(Gamma/2 + 2 gamma_phi + i delta) rho_eg
    !>            d rho_ge = -(Gamma/2 + 2 gamma_phi - i delta) rho_ge
    !>            d rho_gg = +gamma_down rho_ee - gamma_up rho_gg
    !>          精确谱 {0, -Gamma, -(Gamma/2+2gamma_phi) -/+ i delta}。
    subroutine lindblad_supermatrix_2level(gamma_down, gamma_up, gamma_phi, &
                                           delta, l_super)
        real(dp), intent(in)  :: gamma_down, gamma_up, gamma_phi, delta
        complex(dp), intent(out) :: l_super(4, 4)

        real(dp) :: s_coh

        s_coh = 0.5_dp * (gamma_down + gamma_up) + 2.0_dp * gamma_phi
        l_super = (0.0_dp, 0.0_dp)
        l_super(1, 1) = -gamma_down
        l_super(1, 4) = gamma_up
        l_super(2, 2) = -s_coh - EYE * delta
        l_super(3, 3) = -s_coh + EYE * delta
        l_super(4, 1) = gamma_down
        l_super(4, 4) = -gamma_up
    end subroutine lindblad_supermatrix_2level

    !> \brief 二能级热化 Lindblad 稳态 (增广线性方程求解 + 闭式对照)
    !> \details 稳态满足 L rho_ss = 0, Tr = 1; 用增广实方程
    !>          [[L, 1], [1^T, 0]] [x; nu] = [0; 1] 经库内复矩阵求逆求解,
    !>          并与闭式 rho_ee = gamma_up/Gamma, rho_gg = gamma_down/Gamma
    !>          对照 (两者一致为验收)。
    subroutine lindblad_steady_state_2level(gamma_down, gamma_up, gamma_phi, &
                                            delta, rho_ss, stat)
        real(dp), intent(in)  :: gamma_down, gamma_up, gamma_phi, delta
        complex(dp), intent(out) :: rho_ss(2, 2)
        integer, intent(out)  :: stat

        complex(dp), allocatable :: l_super(:, :), aug(:, :), aug_inv(:, :)
        real(dp) :: gam_tot
        integer :: i, j, ierr

        allocate(l_super(4, 4), aug(5, 5), aug_inv(5, 5))
        call lindblad_supermatrix_2level(gamma_down, gamma_up, gamma_phi, &
                                         delta, l_super)
        ! 增广: [L, 1; 1^T, 0] (5x5 复), 右端 [0; 1]
        aug = (0.0_dp, 0.0_dp)
        aug(1:4, 1:4) = l_super
        do i = 1, 4
            aug(i, 5) = (1.0_dp, 0.0_dp)
            aug(5, i) = (1.0_dp, 0.0_dp)
        end do
        call inv_complex_matrix(5, aug, aug_inv, ierr)
        if (ierr /= 0) then
            stat = -1
            rho_ss = (0.0_dp, 0.0_dp)
            return
        end if
        ! x = aug_inv . [0,0,0,0; 1] = aug_inv(:, 5)
        rho_ss = (0.0_dp, 0.0_dp)
        rho_ss(1, 1) = aug_inv(1, 5)
        rho_ss(1, 2) = aug_inv(2, 5)
        rho_ss(2, 1) = aug_inv(3, 5)
        rho_ss(2, 2) = aug_inv(4, 5)
        ! 闭式对照
        gam_tot = gamma_down + gamma_up
        stat = 0
        if (abs(rho_ss(1, 1) - gamma_up / gam_tot) > 1.0e-10_dp .or. &
            abs(rho_ss(2, 2) - gamma_down / gam_tot) > 1.0e-10_dp) stat = -2
        deallocate(l_super, aug, aug_inv)
    end subroutine lindblad_steady_state_2level

    !> \brief 二能级 Lindblad 的阻尼基谱分解 (解析本征对 + 机器验证)
    !> \details 精确谱 lambda = {0, -Gamma, -(Gamma/2+2gamma_phi) -/+ i delta},
    !>          右本征密度矩阵 (阻尼基):
    !>            R_ss  = rho_ss                          (lambda = 0)
    !>            R_pop = |e><e| - |g><g|                  (lambda = -Gamma, 无迹)
    !>            R_eg  = |e><g|                           (lambda = -(s) - i delta)
    !>            R_ge  = |g><e|                           (lambda = -(s) + i delta)
    !>          验证: || L R_k - lambda_k R_k || < 1e-12 (机器精度)。
    !>          无一般复本征求解器依赖 —— 解析对由理论给出, 数值验证其谱关系。
    subroutine lindblad_eigenmodes_2level(gamma_down, gamma_up, gamma_phi, &
                                          delta, rho_ss, evals, modes, stat)
        real(dp), intent(in)  :: gamma_down, gamma_up, gamma_phi, delta
        complex(dp), intent(in)  :: rho_ss(2, 2)
        complex(dp), intent(out) :: evals(4)
        complex(dp), intent(out) :: modes(2, 2, 4)
        integer, intent(out)  :: stat

        complex(dp), allocatable :: l_super(:, :), prod(:, :)
        real(dp) :: s_coh
        integer :: k, i, j

        allocate(l_super(4, 4), prod(4, 1))
        call lindblad_supermatrix_2level(gamma_down, gamma_up, gamma_phi, &
                                         delta, l_super)
        s_coh = 0.5_dp * (gamma_down + gamma_up) + 2.0_dp * gamma_phi

        evals(1) = (0.0_dp, 0.0_dp)
        evals(2) = -(gamma_down + gamma_up)
        evals(3) = -s_coh - EYE * delta
        evals(4) = -s_coh + EYE * delta

        modes(:, :, 1) = rho_ss
        modes(:, :, 2) = (0.0_dp, 0.0_dp)
        modes(1, 1, 2) = (1.0_dp, 0.0_dp)
        modes(2, 2, 2) = -(1.0_dp, 0.0_dp)
        modes(:, :, 3) = (0.0_dp, 0.0_dp)
        modes(1, 2, 3) = (1.0_dp, 0.0_dp)
        modes(:, :, 4) = (0.0_dp, 0.0_dp)
        modes(2, 1, 4) = (1.0_dp, 0.0_dp)

        ! 机器验证: || L vec(R_k) - lambda_k vec(R_k) || = 0
        ! 向量化约定 (ee, eg, ge, gg) = [R11, R12, R21, R22] (行主序展开)
        stat = 0
        do k = 1, 4
            do i = 1, 2
                do j = 1, 2
                    prod((i - 1) * 2 + j, 1) = modes(i, j, k)
                end do
            end do
            prod = matmul(l_super, prod(1:4, 1:1)) &
                   - evals(k) * prod(1:4, 1:1)
            if (abs(sum(abs(prod))) > 1.0e-12_dp) stat = -k
        end do
        deallocate(l_super, prod)
    end subroutine lindblad_eigenmodes_2level

    !> \brief 阻尼基谱重构: rho(t) = rho_ss + sum_k c_k exp(lambda_k t) R_k
    !> \details c_eg = rho_eg(0), c_ge = rho_ge(0), c_pop = rho_ee(0) - rho_ee^ss
    !>          (R_ss 迹 1, R_pop 无迹, 对角分解系数即此)。谱重构与 RK4 推进
    !>          两条独立方法交叉验证 (测试断言)。
    subroutine lindblad_spectral_reconstruct(rho0, rho_ss, gamma_down, gamma_up, &
                                             gamma_phi, delta, t_now, rho_t)
        complex(dp), intent(in) :: rho0(2, 2), rho_ss(2, 2)
        real(dp), intent(in)  :: gamma_down, gamma_up, gamma_phi, delta, t_now
        complex(dp), intent(out) :: rho_t(2, 2)

        real(dp) :: s_coh, c_pop
        complex(dp) :: c_eg, c_ge, phase

        s_coh = 0.5_dp * (gamma_down + gamma_up) + 2.0_dp * gamma_phi
        c_pop = real(rho0(1, 1), dp) - real(rho_ss(1, 1), dp)
        c_eg = rho0(1, 2)
        c_ge = rho0(2, 1)

        phase = exp(-(s_coh + EYE * delta) * t_now)
        rho_t = rho_ss
        rho_t(1, 1) = rho_t(1, 1) + c_pop * exp(-(gamma_down + gamma_up) * t_now)
        rho_t(2, 2) = rho_t(2, 2) - c_pop * exp(-(gamma_down + gamma_up) * t_now)
        rho_t(1, 2) = rho_t(1, 2) + c_eg * phase
        rho_t(2, 1) = rho_t(2, 1) + c_ge * conjg(phase)
    end subroutine lindblad_spectral_reconstruct

end module mod_open_quantum
