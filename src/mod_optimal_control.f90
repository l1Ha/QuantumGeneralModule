!> \brief 量子最优控制理论 (Quantum Optimal Control - Krotov 算法) 模块
!> \details 提供基于 Krotov 迭代算法的激光电场逆向设计求解器。
!>          根据指定目标态保真度与能量惩罚约束，自动反求实现任意量子态高精度跃迁的最优超快光场。
!> \author LiHao
module mod_optimal_control
    use mod_constants, only: dp, EYE, PI
    implicit none
    private

    public :: oct_fidelity
    public :: oct_shape_function
    public :: oct_krotov_step
    public :: oct_optimize_pulse

contains

    !> \brief 计算终态与目标态的跃迁保真度 Fidelity = |<phi_tgt | psi_final>|^2 ∈ [0, 1]
    pure function oct_fidelity(psi_final, phi_tgt) result(fid)
        complex(dp), intent(in) :: psi_final(:), phi_tgt(:)
        real(dp) :: fid
        complex(dp) :: overlap

        overlap = sum(conjg(phi_tgt(:)) * psi_final(:))
        fid = abs(overlap)**2
    end function oct_fidelity

    !> \brief 构造 Krotov 脉冲两端平滑开关形状约束函数 S(t) ∈ [0, 1]
    pure function oct_shape_function(t, t_total, t_ramp) result(s_val)
        real(dp), intent(in) :: t, t_total, t_ramp
        real(dp) :: s_val, t_r

        t_r = max(1.0e-12_dp, t_ramp)
        if (t <= 0.0_dp .or. t >= t_total) then
            s_val = 0.0_dp
        else if (t < t_r) then
            s_val = sin(0.5_dp * PI * t / t_r)**2
        else if (t > t_total - t_r) then
            s_val = sin(0.5_dp * PI * (t_total - t) / t_r)**2
        else
            s_val = 1.0_dp
        end if
    end function oct_shape_function

    !> \brief 执行一次完整的 Krotov 前向-后向迭代更新
    !> \param[in] n_states 系统量子态维度
    !> \param[in] nt 时间步数
    !> \param[in] dt 时间步长 (a.u.)
    !> \param[in] h_diag 无场哈密顿量本征对角元 (n_states)
    !> \param[in] dip_mat 跃迁偶极矩矩阵 (n_states, n_states)
    !> \param[in] psi_ini 初始态波函数 (n_states)
    !> \param[in] phi_tgt 目标态波函数 (n_states)
    !> \param[in] alpha_0 场强能量惩罚权重
    !> \param[in] t_ramp 脉冲两端开启/关闭上升沿 (a.u.)
    !> \param[inout] e_field 待优化的离散激光电场数组 (nt)
    !> \param[out] current_fid 本次迭代后的目标态保真度
    subroutine oct_krotov_step(n_states, nt, dt, h_diag, dip_mat, psi_ini, phi_tgt, &
                               alpha_0, t_ramp, e_field, current_fid)
        integer, intent(in) :: n_states, nt
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: h_diag(n_states)
        real(dp), intent(in) :: dip_mat(n_states, n_states)
        complex(dp), intent(in) :: psi_ini(n_states), phi_tgt(n_states)
        real(dp), intent(in) :: alpha_0, t_ramp
        real(dp), intent(inout) :: e_field(nt)
        real(dp), intent(out) :: current_fid

        complex(dp), allocatable :: chi_hist(:, :)
        complex(dp) :: d_mu, c_curr(n_states)
        real(dp) :: t, t_tot, s_t, delta_e, norm_sq
        integer :: it, k

        allocate(chi_hist(nt, n_states))
        t_tot = real(nt - 1, dp) * dt

        ! 1. 设置后向伴随状态终端条件: |chi(T)> = |phi_tgt>
        c_curr = phi_tgt
        chi_hist(nt, :) = c_curr

        ! 2. 后向反向推进 chi(t) (时间步为 -dt，辅以保范归一化防止累积发散)
        do it = nt, 2, -1
            call propagate_step(c_curr, h_diag, dip_mat, e_field(it), -dt)
            norm_sq = sum(abs(c_curr)**2)
            if (norm_sq > 0.0_dp) c_curr = c_curr / sqrt(norm_sq)
            chi_hist(it - 1, :) = c_curr
        end do

        ! 3. 前向并发推进新波包并原位更新激光电场:
        !    Delta E(t) = - (S(t) / alpha_0) * Im[ <chi(t) | mu | psi_new(t)> ]
        c_curr = psi_ini
        do it = 1, nt
            t = real(it - 1, dp) * dt
            s_t = oct_shape_function(t, t_tot, t_ramp)

            ! 计算偶极矩阵元期望 <chi(t) | mu | psi(t)>
            d_mu = (0.0_dp, 0.0_dp)
            do k = 1, n_states
                d_mu = d_mu + conjg(chi_hist(it, k)) * sum(dip_mat(k, :) * c_curr(:))
            end do

            delta_e = -(s_t / max(1.0e-3_dp, alpha_0)) * aimag(d_mu)
            e_field(it) = e_field(it) + delta_e

            if (it < nt) then
                call propagate_step(c_curr, h_diag, dip_mat, e_field(it), dt)
                norm_sq = sum(abs(c_curr)**2)
                if (norm_sq > 0.0_dp) c_curr = c_curr / sqrt(norm_sq)
            end if
        end do

        current_fid = oct_fidelity(c_curr, phi_tgt)
        deallocate(chi_hist)
    end subroutine oct_krotov_step

    !> \brief 高级最优控制主循环：迭代直至达到目标保真度或最大迭代次数
    subroutine oct_optimize_pulse(n_states, nt, dt, h_diag, dip_mat, psi_ini, phi_tgt, &
                                  alpha_0, t_ramp, max_iters, target_fid, e_field, &
                                  final_fid, actual_iters)
        integer, intent(in) :: n_states, nt
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: h_diag(n_states)
        real(dp), intent(in) :: dip_mat(n_states, n_states)
        complex(dp), intent(in) :: psi_ini(n_states), phi_tgt(n_states)
        real(dp), intent(in) :: alpha_0, t_ramp
        integer, intent(in) :: max_iters
        real(dp), intent(in) :: target_fid
        real(dp), intent(inout) :: e_field(nt)
        real(dp), intent(out) :: final_fid
        integer, intent(out) :: actual_iters

        integer :: iter
        real(dp) :: fid

        final_fid = 0.0_dp
        actual_iters = 0

        do iter = 1, max_iters
            call oct_krotov_step(n_states, nt, dt, h_diag, dip_mat, psi_ini, phi_tgt, &
                                 alpha_0, t_ramp, e_field, fid)
            actual_iters = iter
            final_fid = fid
            if (fid >= target_fid) exit
        end do
    end subroutine oct_optimize_pulse

    !> \brief 4 阶 Runge-Kutta 单步推进子程序 (支持正向 dt > 0 与反向 dt < 0)
    pure subroutine propagate_step(c_vec, diag_e, dip_m, e_t, dt_step)
        complex(dp), intent(inout) :: c_vec(:)
        real(dp), intent(in) :: diag_e(:), dip_m(:, :), e_t, dt_step

        complex(dp) :: k1(size(c_vec)), k2(size(c_vec)), k3(size(c_vec)), k4(size(c_vec))
        complex(dp) :: c_tmp(size(c_vec))

        k1 = deriv_c(c_vec, diag_e, dip_m, e_t)
        c_tmp = c_vec + 0.5_dp * dt_step * k1
        k2 = deriv_c(c_tmp, diag_e, dip_m, e_t)
        c_tmp = c_vec + 0.5_dp * dt_step * k2
        k3 = deriv_c(c_tmp, diag_e, dip_m, e_t)
        c_tmp = c_vec + dt_step * k3
        k4 = deriv_c(c_tmp, diag_e, dip_m, e_t)

        c_vec = c_vec + (dt_step / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end subroutine propagate_step

    pure function deriv_c(c_in, diag_e, dip_m, e_t) result(dcdt)
        complex(dp), intent(in) :: c_in(:)
        real(dp), intent(in) :: diag_e(:), dip_m(:, :), e_t
        complex(dp) :: dcdt(size(c_in))
        integer :: i_idx, j_idx, dim_n
        complex(dp) :: h_c

        dim_n = size(c_in)
        do i_idx = 1, dim_n
            h_c = diag_e(i_idx) * c_in(i_idx)
            do j_idx = 1, dim_n
                if (dip_m(i_idx, j_idx) /= 0.0_dp) then
                    h_c = h_c - e_t * dip_m(i_idx, j_idx) * c_in(j_idx)
                end if
            end do
            dcdt(i_idx) = -EYE * h_c
        end do
    end function deriv_c

end module mod_optimal_control
