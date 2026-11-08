!> \brief 开放量子系统与 Lindblad 耗散主方程动力学模块
!> \details 提供密度矩阵表示、自发辐射弛豫 (T_1)、纯退相位 (T_2^*) 及热浴耦合 Jump 算符构建、
!>          Runge-Kutta 密度矩阵含时主方程推进器，以及量子纯度、相干度与冯·诺依曼熵诊断例程；
!>          Lindblad Liouvillian 谱方法 (多原子拓展第九步 / 开放系统谱方法):
!>          向量化超算符、稳态 (增广线性方程 + 闭式)、阻尼基谱分解与谱重构;
!>          驱动耗散二能级: 光学 Bloch 稳态 (饱和闭式)、量子回归定理一阶
!>          关联函数与 Mollow 三重峰共振荧光谱 (多原子拓展第十步);
!>          通用 n 能级 Lindblad 稳态求解器 (多原子拓展第十一步);
!>          偏迹与纠缠熵 (量子信息工具箱, Araki-Lieb 不等式机器验证)。
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
    public :: optical_bloch_steady_state_driven
    public :: first_order_correlation_2level
    public :: mollow_spectrum_grid
    public :: lindblad_supermatrix_general
    public :: lindblad_steady_state_general
    public :: partial_trace_b
    public :: entanglement_entropy
    public :: lindblad_supermatrix_driven

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

    !> \brief 驱动耗散二能级的向量化 Liouvillian 超矩阵 (4x4 复, 行主序)
    !> \details 旋转框架 H = -delta sigma_z/2 + Omega sigma_x/2, 塌缩 Gamma
    !>          (sigma^-); 向量化约定 vec(i,j) = (i-1)*2 + j (行主序,
    !>          与第五步阻尼基一致)。稳态与 QRT 传播共用。
    subroutine lindblad_supermatrix_driven(omega, delta, gamma, l_super)
        real(dp), intent(in)  :: omega, delta, gamma
        complex(dp), intent(out) :: l_super(4, 4)

        complex(dp) :: h_rot(2, 2), l_m(2, 2), rho_p(2, 2)
        complex(dp) :: l_minus(2, 2)
        integer :: i, j, k, l

        h_rot = (0.0_dp, 0.0_dp)
        h_rot(1, 1) = -0.5_dp * delta
        h_rot(2, 2) = 0.5_dp * delta
        h_rot(1, 2) = 0.5_dp * omega
        h_rot(2, 1) = 0.5_dp * omega
        l_minus = (0.0_dp, 0.0_dp)
        l_minus(2, 1) = (1.0_dp, 0.0_dp)

        do j = 1, 2
            do i = 1, 2
                rho_p = (0.0_dp, 0.0_dp)
                rho_p(i, j) = (1.0_dp, 0.0_dp)
                l_m = -EYE * (matmul(h_rot, rho_p) - matmul(rho_p, h_rot)) &
                      + gamma * lindblad_dissipator(l_minus, rho_p)
                do l = 1, 2
                    do k = 1, 2
                        l_super((k - 1) * 2 + l, (i - 1) * 2 + j) = l_m(k, l)
                    end do
                end do
            end do
        end do
    end subroutine lindblad_supermatrix_driven

    !> \brief 驱动耗散二能级的稳态 (旋转框架, 增广求解 + 饱和闭式对照)
    !> \details 旋转框架 H = -delta sigma_z/2 + Omega sigma_x/2, 塌缩 Gamma
    !>          (sigma^-), 退相位 0 (Gamma 已含相干衰减 Gamma/2):
    !>          饱和参数 s = 2 Omega^2 / (Gamma^2 + 4 delta^2),
    !>          闭式 rho_ee^ss = s/(2(1+s))。
    !>          增广方程 [[L, 1], [1^T, 0]] 经库内复矩阵求逆求解并与闭式对照。
    subroutine optical_bloch_steady_state_driven(omega, delta, gamma, &
                                                 rho_ss, stat)
        real(dp), intent(in)  :: omega, delta, gamma
        complex(dp), intent(out) :: rho_ss(2, 2)
        integer, intent(out)  :: stat

        complex(dp), allocatable :: l_super(:, :), aug(:, :), aug_inv(:, :)
        real(dp) :: s_sat
        integer :: i, j, ierr

        allocate(l_super(4, 4), aug(5, 5), aug_inv(5, 5))
        call lindblad_supermatrix_driven(omega, delta, gamma, l_super)

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
            deallocate(l_super, aug, aug_inv)
            return
        end if
        rho_ss = (0.0_dp, 0.0_dp)
        rho_ss(1, 1) = aug_inv(1, 5)
        rho_ss(1, 2) = aug_inv(2, 5)
        rho_ss(2, 1) = aug_inv(3, 5)
        rho_ss(2, 2) = aug_inv(4, 5)
        s_sat = 2.0_dp * omega**2 / (gamma**2 + 4.0_dp * delta**2)
        stat = 0
        if (abs(real(rho_ss(1, 1), dp) - s_sat / (2.0_dp * (1.0_dp + s_sat))) > 1.0e-10_dp) &
            stat = -2
        deallocate(l_super, aug, aug_inv)
    end subroutine optical_bloch_steady_state_driven

    !> \brief 量子回归定理一阶关联函数 g1(tau) = <sigma^+(tau) sigma^-(0)>
    !> \details QRT: g1(tau) = tr[sigma^+ · e^{L tau}(sigma^- rho_ss)]; 传播用
    !>          库内 rk4_lindblad_step。弱驱动单 Lorentz 峰 (HWHM = Gamma/2),
    !>          强驱动 Mollow 三重峰 (边带 = Liouvillian 本征值虚部)。
    subroutine first_order_correlation_2level(omega, delta, gamma, rho_ss, &
                                              tau_max, n_tau, tau_grid, g1_grid, stat)
        real(dp), intent(in)  :: omega, delta, gamma, tau_max
        complex(dp), intent(in) :: rho_ss(2, 2)
        integer, intent(in)   :: n_tau
        real(dp), intent(out) :: tau_grid(n_tau), g1_grid(n_tau)
        integer, intent(out)  :: stat

        complex(dp), allocatable :: l_super(:, :), m_vec(:, :), k1(:, :), &
                                    k2(:, :), k3(:, :), k4(:, :), m_tmp(:, :)
        real(dp) :: dt
        complex(dp) :: m0(2, 2)
        integer :: i, istep, n_sub, r_i, c_i

        stat = 0
        if (n_tau < 2) then
            stat = -1
            return
        end if

        ! QRT: M(0) = sigma^- rho_ss; M_vec 传播 d vec(M)/dt = L vec(M)
        ! (超矩阵向量 RK4, 无厄米化/迹归一 —— M(tau) 不是密度矩阵!)
        allocate(l_super(4, 4), m_vec(4, 1), k1(4, 1), k2(4, 1), k3(4, 1), &
                 k4(4, 1), m_tmp(4, 1))
        call lindblad_supermatrix_driven(omega, delta, gamma, l_super)
        m0 = reshape((/ (0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), &
                       (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp) /), (/ 2, 2 /))
        m0 = matmul(m0, rho_ss)

        do i = 1, n_tau
            tau_grid(i) = tau_max * real(i - 1, dp) / real(n_tau - 1, dp)
        end do
        ! 初始向量 (行主序 vec(i,j) = (i-1)*2 + j)
        do r_i = 1, 2
            do c_i = 1, 2
                m_vec((r_i - 1) * 2 + c_i, 1) = m0(r_i, c_i)
            end do
        end do
        g1_grid(1) = m_vec(3, 1)              ! <sigma+ M> = M(2,1)

        n_sub = 40
        dt = tau_max / real(n_tau - 1, dp) / real(n_sub, dp)
        do i = 2, n_tau
            do istep = 1, n_sub
                k1 = matmul(l_super, m_vec)
                k2 = matmul(l_super, m_vec + 0.5_dp * dt * k1)
                k3 = matmul(l_super, m_vec + 0.5_dp * dt * k2)
                k4 = matmul(l_super, m_vec + dt * k3)
                m_vec = m_vec + (dt / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
            end do
            g1_grid(i) = m_vec(3, 1)          ! M(2,1): <sigma+ M>
        end do
        deallocate(l_super, m_vec, k1, k2, k3, k4, m_tmp)
    end subroutine first_order_correlation_2level

    !> \brief Mollow 共振荧光谱: Wiener-Khinchin 傅里叶变换 g1(tau)
    !> \details S(omega) ∝ Re int_0^inf g1(tau) exp(i omega tau) dtau
    !>          (g1 以 exp(-Gamma tau/2) 级衰减保证收敛)。
    !>          归一化: sum S domega = g1(0) = rho_ee^ss。
    subroutine mollow_spectrum_grid(omega, delta, gamma, rho_ss, tau_max, n_tau, &
                                    n_om, om_grid, s_grid, stat)
        real(dp), intent(in)  :: omega, delta, gamma, tau_max
        complex(dp), intent(in) :: rho_ss(2, 2)
        integer, intent(in)   :: n_tau, n_om
        real(dp), intent(out) :: om_grid(n_om), s_grid(n_om)
        integer, intent(out)  :: stat

        real(dp), allocatable :: tau_grid(:), g1_grid(:)
        real(dp) :: d_tau, om_lo, om_hi, d_om, re_part, dpi
        complex(dp) :: ph
        integer :: i, k

        stat = 0
        dpi = acos(-1.0_dp)
        if (n_tau < 2 .or. n_om < 2) then
            stat = -1
            return
        end if
        allocate(tau_grid(n_tau), g1_grid(n_tau))
        call first_order_correlation_2level(omega, delta, gamma, rho_ss, &
                                            tau_max, n_tau, tau_grid, g1_grid, stat)
        d_tau = tau_grid(2) - tau_grid(1)
        om_lo = -2.0_dp * (omega + 5.0_dp * gamma)
        om_hi = 2.0_dp * (omega + 5.0_dp * gamma)
        d_om = (om_hi - om_lo) / real(n_om - 1, dp)
        do k = 1, n_om
            om_grid(k) = om_lo + real(k - 1, dp) * d_om
            re_part = 0.0_dp
            do i = 1, n_tau
                ph = exp((0.0_dp, 1.0_dp) * om_grid(k) * tau_grid(i))
                re_part = re_part + real(g1_grid(i) * ph, dp) * &
                    merge(0.5_dp, 1.0_dp, i == 1 .or. i == n_tau)
            end do
            s_grid(k) = re_part * d_tau / dpi
        end do
        deallocate(tau_grid, g1_grid)
    end subroutine mollow_spectrum_grid

    !> \brief 通用 n 能级 Lindblad 向量化超矩阵 (n^2 x n^2 复)
    !> \details 从任意哈密顿量与任意塌缩通道列表构建完整 Liouvillian:
    !>          L(rho) = -i[H, rho] + sum_k gamma_k D[L_k] rho。
    !>          向量化约定: vec(rho)_{(i-1)*n+j} = rho(i,j) (行主序)。
    subroutine lindblad_supermatrix_general(h_mat, l_ops, gamma_rates, &
                                            n_dim, n_ops, l_super, stat)
        real(dp), intent(in)  :: h_mat(:, :)
        complex(dp), intent(in) :: l_ops(:, :, :)
        real(dp), intent(in)  :: gamma_rates(:)
        integer, intent(in)   :: n_dim, n_ops
        complex(dp), intent(out) :: l_super(n_dim * n_dim, n_dim * n_dim)
        integer, intent(out)  :: stat

        complex(dp) :: e_mat(n_dim, n_dim), l_m(n_dim, n_dim), rho_p(n_dim, n_dim)
        integer :: i, j, k, l, ch, vi, vj

        stat = 0
        if (n_dim < 2 .or. n_ops < 1) then
            stat = -1
            return
        end if
        l_super = (0.0_dp, 0.0_dp)

        do j = 1, n_dim
            do i = 1, n_dim
                e_mat = (0.0_dp, 0.0_dp)
                e_mat(i, j) = (1.0_dp, 0.0_dp)
                ! Hamiltonian 部分
                l_m = -EYE * (matmul(h_mat, e_mat) - matmul(e_mat, h_mat))
                ! Lindblad 耗散部分
                do ch = 1, n_ops
                    if (gamma_rates(ch) > 0.0_dp) then
                        l_m = l_m + gamma_rates(ch) * lindblad_dissipator( &
                            l_ops(:, :, ch), e_mat)
                    end if
                end do
                ! 散射到超矩阵 (行主序 vec(i,j) = (i-1)*n + j)
                do l = 1, n_dim
                    do k = 1, n_dim
                        vi = (k - 1) * n_dim + l
                        vj = (i - 1) * n_dim + j
                        l_super(vi, vj) = l_m(k, l)
                    end do
                end do
            end do
        end do
    end subroutine lindblad_supermatrix_general

    !> \brief 通用 n 能级 Lindblad 稳态 (增广线性方程 + 性质自检)
    !> \details 增广系统 [[L, I_vec], [I_vec^T, 0]] [x; nu] = [0; 1] 经库内
    !>          复矩阵求逆求解 (n^2 + 1 阶)。性质自检: Tr = 1、厄米性、
    !>          L rho_ss = 0 残差。
    !> \note   复杂度 O(n^6) (求逆) —— n <= 6 建议 (n^2+1 <= 37)。
    subroutine lindblad_steady_state_general(h_mat, l_ops, gamma_rates, &
                                             n_dim, n_ops, rho_ss, stat)
        real(dp), intent(in)  :: h_mat(:, :)
        complex(dp), intent(in) :: l_ops(:, :, :)
        real(dp), intent(in)  :: gamma_rates(:)
        integer, intent(in)   :: n_dim, n_ops
        complex(dp), intent(out), allocatable :: rho_ss(:, :)
        integer, intent(out)  :: stat

        complex(dp), allocatable :: l_super(:, :), aug(:, :), aug_inv(:, :)
        complex(dp), allocatable :: resid(:, :)
        integer :: n_aug, i, j, ierr
        real(dp) :: tr_val

        stat = 0
        if (n_dim < 2 .or. n_ops < 1) then
            stat = -1
            allocate(rho_ss(0, 0))
            return
        end if
        n_aug = n_dim * n_dim + 1
        allocate(l_super(n_dim * n_dim, n_dim * n_dim), aug(n_aug, n_aug), &
                 aug_inv(n_aug, n_aug), rho_ss(n_dim, n_dim), resid(n_dim, n_dim))

        call lindblad_supermatrix_general(h_mat, l_ops, gamma_rates, &
                                          n_dim, n_ops, l_super, stat)
        if (stat /= 0) return

        aug = (0.0_dp, 0.0_dp)
        aug(1:n_dim*n_dim, 1:n_dim*n_dim) = l_super
        do i = 1, n_dim
            aug((i - 1) * n_dim + i, n_aug) = (1.0_dp, 0.0_dp)
            aug(n_aug, (i - 1) * n_dim + i) = (1.0_dp, 0.0_dp)
        end do
        call inv_complex_matrix(n_aug, aug, aug_inv, ierr)
        if (ierr /= 0) then
            stat = -2
            deallocate(l_super, aug, aug_inv, rho_ss, resid)
            return
        end if
        rho_ss = (0.0_dp, 0.0_dp)
        do j = 1, n_dim
            do i = 1, n_dim
                rho_ss(i, j) = aug_inv((i - 1) * n_dim + j, n_aug)
            end do
        end do

        ! 性质自检: Tr = 1
        tr_val = 0.0_dp
        do i = 1, n_dim
            tr_val = tr_val + real(rho_ss(i, i), dp)
        end do
        if (abs(tr_val - 1.0_dp) > 1.0e-10_dp) stat = -3
        ! 厄米性
        do i = 1, n_dim
            do j = i + 1, n_dim
                if (abs(rho_ss(i, j) - conjg(rho_ss(j, i))) > 1.0e-10_dp) stat = -4
            end do
        end do
        ! L rho_ss = 0 残差
        do j = 1, n_dim
            do i = 1, n_dim
                resid(i, j) = sum(l_super((i-1)*n_dim+j, 1:n_dim*n_dim) &
                    * reshape(rho_ss, (/ n_dim * n_dim /)))
            end do
        end do
        if (abs(sum(abs(resid))) > 1.0e-10_dp) stat = -5

        deallocate(l_super, aug, aug_inv, resid)
    end subroutine lindblad_steady_state_general

    !> \brief 偏迹: rho_A = Tr_B(rho_AB) (dim_A × dim_B → dim_A × dim_A)
    !> \details 行主序约定: rho_AB((i_a-1)*dim_b+i_b, (j_a-1)*dim_b+j_b)
    !>          = <i_a i_b| rho_AB |j_a j_b>。
    subroutine partial_trace_b(rho_ab, dim_a, dim_b, rho_a)
        complex(dp), intent(in)  :: rho_ab(:, :)
        integer, intent(in)      :: dim_a, dim_b
        complex(dp), intent(out) :: rho_a(dim_a, dim_a)

        integer :: ia, ib, ja, jb

        rho_a = (0.0_dp, 0.0_dp)
        do ja = 1, dim_a
            do ia = 1, dim_a
                do ib = 1, dim_b
                    rho_a(ia, ja) = rho_a(ia, ja) + &
                        rho_ab((ia - 1) * dim_b + ib, (ja - 1) * dim_b + ib)
                end do
            end do
        end do
    end subroutine partial_trace_b

    !> \brief 纯态双粒子系统的纠缠熵 S = -Tr(rho_A ln rho_A)
    !> \details |psi> ∈ A ⊗ B 纯态 → S(rho_A) = S(rho_B) = 纠缠熵
    !>          (von Neumann 熵的约化密度矩阵); 通过 rho_A 的本征值计算。
    subroutine entanglement_entropy(psi, dim_a, dim_b, s_ent)
        complex(dp), intent(in)  :: psi(:)
        integer, intent(in)      :: dim_a, dim_b
        real(dp), intent(out)    :: s_ent

        complex(dp), allocatable :: rho_a(:, :)
        real(dp), allocatable :: ra(:, :), ed(:), ez(:, :)
        real(dp) :: p_val
        integer :: ia, ib, ja, stat

        allocate(rho_a(dim_a, dim_a), ra(dim_a, dim_a), ed(dim_a), ez(dim_a, dim_a))
        ! 从态矢量构建 rho_A = Tr_B(|psi><psi|)
        rho_a = (0.0_dp, 0.0_dp)
        do ja = 1, dim_a
            do ia = 1, dim_a
                do ib = 1, dim_b
                    rho_a(ia, ja) = rho_a(ia, ja) + &
                        conjg(psi((ia - 1) * dim_b + ib)) * psi((ja - 1) * dim_b + ib)
                end do
            end do
        end do

        ! von Neumann 熵
        ra = real(rho_a, dp)
        call diag_symmetric_matrix(dim_a, ra, ed, ez, stat)
        s_ent = 0.0_dp
        do ia = 1, dim_a
            p_val = ed(ia)
            if (p_val > 1.0e-14_dp) s_ent = s_ent - p_val * log(p_val)
        end do
        deallocate(rho_a, ra, ed, ez)
    end subroutine entanglement_entropy

end module mod_open_quantum
