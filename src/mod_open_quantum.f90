!> \brief 开放量子系统与 Lindblad 耗散主方程动力学模块
!> \details 提供密度矩阵表示、自发辐射弛豫 (T_1)、纯退相位 (T_2^*) 及热浴耦合 Jump 算符构建、
!>          Runge-Kutta 密度矩阵含时主方程推进器，以及量子纯度、相干度与冯·诺依曼熵诊断例程。
!> \author LiHao
module mod_open_quantum
    use mod_constants, only: dp, EYE, KB
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: lindblad_dissipator
    public :: rk4_lindblad_step
    public :: calculate_quantum_purity
    public :: calculate_von_neumann_entropy
    public :: calculate_quantum_coherence
    public :: create_relaxation_jump_op
    public :: create_dephasing_jump_op

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

end module mod_open_quantum
