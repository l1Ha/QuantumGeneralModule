!> \brief 热统计力学与玻尔兹曼系综平均算法模块
!> \details 计算有限温度 T 下刚体转动与振动能级的玻尔兹曼热初态布居权重，
!>          实现宏观可观测期望值（定向度、排列度、电离产额）的系综加权统计平均，
!>          并提供开放系统玻色-爱因斯坦（Bose-Einstein）热浴关联函数。
!> \author LiHao
module mod_thermal_ensemble
    use mod_constants, only: dp, KB, K2AU
    implicit none
    private

    public :: boltzmann_rotational_weights
    public :: boltzmann_vibrational_weights
    public :: thermal_average_1d
    public :: thermal_average_2d
    public :: bose_einstein_factor

contains

    !> \brief 计算转动能级 J 的玻尔兹曼统计权重: w_J = (2J+1)*exp(-B*J*(J+1)/k_B*T) / Z_rot
    !> \param[in] temp_k 气体温度 (K)
    !> \param[in] b_rot 分子转动常数 B_e (a.u.)
    !> \param[in] j_max 最大转动量子数截断
    !> \param[out] weights 归一化统计权重向量 (0:j_max)
    !> \param[out] z_rot 配分函数
    subroutine boltzmann_rotational_weights(temp_k, b_rot, j_max, weights, z_rot)
        real(dp), intent(in) :: temp_k
        real(dp), intent(in) :: b_rot
        integer, intent(in) :: j_max
        real(dp), intent(out) :: weights(0:j_max)
        real(dp), intent(out), optional :: z_rot

        integer :: j
        real(dp) :: beta, energy_j, z_sum

        if (temp_k <= 1.0e-5_dp) then
            ! T = 0 K 纯态基态
            weights = 0.0_dp
            weights(0) = 1.0_dp
            if (present(z_rot)) z_rot = 1.0_dp
            return
        end if

        ! k_B * T 转换为原子单位
        beta = 1.0_dp / (temp_k * K2AU)
        z_sum = 0.0_dp

        do j = 0, j_max
            energy_j = b_rot * real(j * (j + 1), dp)
            weights(j) = real(2 * j + 1, dp) * exp(-beta * energy_j)
            z_sum = z_sum + weights(j)
        end do

        ! 归一化: sum w_J = 1
        if (z_sum > 0.0_dp) then
            weights = weights / z_sum
        end if

        if (present(z_rot)) z_rot = z_sum
    end subroutine boltzmann_rotational_weights

    !> \brief 计算简谐振动能级 v 的玻尔兹曼统计权重
    subroutine boltzmann_vibrational_weights(temp_k, omega_e, v_max, weights)
        real(dp), intent(in) :: temp_k
        real(dp), intent(in) :: omega_e
        integer, intent(in) :: v_max
        real(dp), intent(out) :: weights(0:v_max)

        integer :: v
        real(dp) :: beta, z_sum

        if (temp_k <= 1.0e-5_dp) then
            weights = 0.0_dp
            weights(0) = 1.0_dp
            return
        end if

        beta = 1.0_dp / (temp_k * K2AU)
        z_sum = 0.0_dp

        do v = 0, v_max
            weights(v) = exp(-beta * omega_e * real(v, dp))
            z_sum = z_sum + weights(v)
        end do

        if (z_sum > 0.0_dp) then
            weights = weights / z_sum
        end if
    end subroutine boltzmann_vibrational_weights

    !> \brief 对一维各态可观测量进行加权热平均: <A>_T = sum_J w_J * A_J
    pure function thermal_average_1d(observables, weights) result(avg_val)
        real(dp), intent(in) :: observables(0:)
        real(dp), intent(in) :: weights(0:)
        real(dp) :: avg_val
        integer :: j, n

        avg_val = 0.0_dp
        n = min(ubound(observables, 1), ubound(weights, 1))
        do j = 0, n
            avg_val = avg_val + weights(j) * observables(j)
        end do
    end function thermal_average_1d

    !> \brief 对含时二维矩阵演化曲线进行加权热平均: <A>(t)_T = sum_J w_J * A_J(t)
    subroutine thermal_average_2d(obs_time_series, weights, avg_series)
        real(dp), intent(in) :: obs_time_series(:, 0:)  !< 矩阵 (n_steps, 0:j_max)
        real(dp), intent(in) :: weights(0:)             !< 权重 (0:j_max)
        real(dp), intent(out) :: avg_series(:)          !< 热平均时序 (n_steps)
        integer :: i, j, n_steps, n_states

        n_steps = size(obs_time_series, 1)
        n_states = min(ubound(obs_time_series, 2), ubound(weights, 1))

        avg_series = 0.0_dp
        do j = 0, n_states
            do i = 1, n_steps
                avg_series(i) = avg_series(i) + weights(j) * obs_time_series(i, j)
            end do
        end do
    end subroutine thermal_average_2d

    !> \brief 玻色-爱因斯坦热库分布因子: n_BE(omega, T) = 1 / (exp(omega / k_B*T) - 1)
    pure function bose_einstein_factor(omega_au, temp_k) result(n_be)
        real(dp), intent(in) :: omega_au
        real(dp), intent(in) :: temp_k
        real(dp) :: n_be
        real(dp) :: x

        if (temp_k <= 1.0e-5_dp .or. omega_au <= 1.0e-12_dp) then
            n_be = 0.0_dp
            return
        end if

        x = omega_au / (temp_k * K2AU)
        if (x > 80.0_dp) then
            n_be = 0.0_dp
        else if (x < 1.0e-5_dp) then
            n_be = 1.0_dp / x  ! Rayleigh-Jeans 极限
        else
            n_be = 1.0_dp / (exp(x) - 1.0_dp)
        end if
    end function bose_einstein_factor

end module mod_thermal_ensemble
