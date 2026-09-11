!> \brief 切比雪夫多项式推进器与能谱滤波算子模块
!> \details 提供时间无关哈密顿量超大步长切比雪夫多项式精确演化（机器精度级范数守恒），
!>          以及 Schafer-Kulander 窗算子能谱投影滤波器（用于光电子动能谱 PES 精确提取）。
!> \author LiHao
module mod_chebyshev_propagator
    use mod_constants, only: dp, PI, TWOPI, EYE
    implicit none
    private

    public :: chebyshev_propagate_step
    public :: window_operator_pes

contains

    !> \brief 切比雪夫多项式展开单步推进演化: psi(t+dt) = exp(-i * H * dt) * psi(t)
    !> \param[inout] psi 复数波函数 (n_dim)
    !> \param[in] dt 时间步长 (a.u.)
    !> \param[in] e_min 哈密顿量谱范围估计下界
    !> \param[in] e_max 哈密顿量谱范围估计上界
    !> \param[in] h_mult 用户自定义哈密顿量矩阵-向量乘法子程序
    !> \param[in] order_cheb 切比雪夫展开截断阶数 (通常 15~40 即可达机器精度)
    subroutine chebyshev_propagate_step(psi, dt, e_min, e_max, h_mult, order_cheb)
        complex(dp), intent(inout) :: psi(:)
        real(dp), intent(in) :: dt, e_min, e_max
        integer, intent(in), optional :: order_cheb
        interface
            subroutine h_mult(v_in, v_out)
                use mod_constants, only: dp
                complex(dp), intent(in) :: v_in(:)
                complex(dp), intent(out) :: v_out(:)
            end subroutine h_mult
        end interface

        integer :: n, m, k
        real(dp) :: e_bar, delta_e, alpha
        real(dp), allocatable :: j_bessel(:)
        complex(dp), allocatable :: phi_0(:), phi_1(:), phi_2(:), psi_accum(:), h_phi(:)
        complex(dp) :: phase_factor

        n = size(psi)
        m = 25
        if (present(order_cheb)) m = max(5, order_cheb)

        ! 1. 谱范围平移归一化到 [-1, 1]: H_norm = (H - E_bar) / Delta_E
        e_bar = 0.5_dp * (e_max + e_min)
        delta_e = 0.5_dp * (e_max - e_min)
        alpha = delta_e * dt

        allocate(phi_0(n), phi_1(n), phi_2(n), psi_accum(n), h_phi(n), j_bessel(0:m))

        ! 计算第 1 类贝塞尔函数 J_k(alpha)
        call compute_bessel_j(alpha, m, j_bessel)

        ! 2. 初始切比雪夫多项式向量:
        ! T_0(H_norm) * psi = psi
        phi_0 = psi
        psi_accum = cmplx(j_bessel(0), 0.0_dp, kind=dp) * phi_0

        ! T_1(H_norm) * psi = -i * H_norm * psi
        call h_mult(phi_0, h_phi)
        phi_1 = -EYE * ((h_phi - e_bar * phi_0) / delta_e)
        psi_accum = psi_accum + 2.0_dp * cmplx(j_bessel(1), 0.0_dp, kind=dp) * phi_1

        ! 3. 递归递推 T_{k+1} = -2*i * H_norm * T_k + T_{k-1}
        do k = 2, m
            call h_mult(phi_1, h_phi)
            phi_2 = -2.0_dp * EYE * ((h_phi - e_bar * phi_1) / delta_e) + phi_0
            psi_accum = psi_accum + 2.0_dp * cmplx(j_bessel(k), 0.0_dp, kind=dp) * phi_2

            phi_0 = phi_1
            phi_1 = phi_2
        end do

        ! 4. 乘以谱平移总局域相位: exp(-i * E_bar * dt)
        phase_factor = exp(-EYE * e_bar * dt)
        psi = phase_factor * psi_accum

        deallocate(phi_0, phi_1, phi_2, psi_accum, h_phi, j_bessel)
    end subroutine chebyshev_propagate_step

    !> \brief 第一类贝塞尔函数 J_0 到 J_m 计算 (Miller 逆向递推算法)
    subroutine compute_bessel_j(x, m, j_arr)
        real(dp), intent(in) :: x
        integer, intent(in) :: m
        real(dp), intent(out) :: j_arr(0:m)

        integer :: l, l_start
        real(dp) :: j_prev, j_curr, j_next, sum_norm

        if (abs(x) < 1.0e-14_dp) then
            j_arr = 0.0_dp
            j_arr(0) = 1.0_dp
            return
        end if

        l_start = max(m + 15, int(abs(x) + 20.0_dp))
        j_next = 0.0_dp
        j_curr = 1.0e-30_dp
        sum_norm = 0.0_dp

        do l = l_start, 1, -1
            j_prev = (2.0_dp * real(l, dp) / x) * j_curr - j_next
            j_next = j_curr
            j_curr = j_prev

            if (l <= m) j_arr(l) = j_next
            if (mod(l, 2) == 0) then
                sum_norm = sum_norm + 2.0_dp * j_next
            end if
        end do

        j_arr(0) = j_curr
        sum_norm = sum_norm + j_curr

        ! 归一化: J_0 + 2*J_2 + 2*J_4 + ... = 1
        if (abs(sum_norm) > 0.0_dp) then
            j_arr = j_arr / sum_norm
        end if
    end subroutine compute_bessel_j

    !> \brief Schafer-Kulander 能量窗算子提取特定动能 E_k 的光电子概率
    !> \details W(E_k, gamma, n) = gamma^(2n) / [ (H - E_k)^(2n) + gamma^(2n) ]
    !> \param[in] psi 末态散射总波函数
    !> \param[in] dx 空间步长
    !> \param[in] e_k 目标电子能量 E_k
    !> \param[in] gamma 能量窗半高宽 (分辨率)
    !> \param[in] h_mat 对角化哈密顿本征能量 eig_vals (n) 与本征矢 eig_vecs (n, n)
    pure function window_operator_pes(psi, dx, e_k, gamma, eig_vals, eig_vecs) result(prob)
        complex(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: dx, e_k, gamma
        real(dp), intent(in) :: eig_vals(:)
        real(dp), intent(in) :: eig_vecs(:, :)
        real(dp) :: prob

        integer :: i, n, v
        complex(dp) :: c_proj
        real(dp) :: weight

        prob = 0.0_dp
        n = size(psi)

        do v = 1, n
            ! 投影系数 c_v = <phi_v | psi>
            c_proj = (0.0_dp, 0.0_dp)
            do i = 1, n
                c_proj = c_proj + eig_vecs(i, v) * psi(i) * dx
            end do

            ! 窗函数洛伦兹权重: gamma^2 / [ (E_v - E_k)^2 + gamma^2 ]
            weight = (gamma**2) / ((eig_vals(v) - e_k)**2 + gamma**2)
            prob = prob + (abs(c_proj)**2) * weight
        end do
    end function window_operator_pes

end module mod_chebyshev_propagator
