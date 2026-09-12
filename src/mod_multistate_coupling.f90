!> \brief 多势能面非绝热耦合核动力学模块
!> \details 提供双通道非绝热耦合薛定谔方程（2-Channel Split-Operator）精确演化、
!>          Landau-Zener 经典跃迁几率计算以及各电子态布居数与分支比监测。
!> \author LiHao
module mod_multistate_coupling
    use mod_constants, only: dp, PI, TWOPI, EYE
    use mod_linear_algebra, only: fft_1d
    implicit none
    private

    public :: propagate_split_operator_2channel
    public :: landau_zener_probability
    public :: calculate_channel_populations

contains

    !> \brief 双态非绝热耦合波包二阶辛对称分裂算符步进
    !> \details 演化方程: i * d/dt [psi1, psi2]^T = [ T_kin * I + V_diabatic(R) ] * [psi1, psi2]^T
    !>          其中势能算子采用 2x2 解析酉矩阵指数精确计算，确保波包绝对幺正守恒。
    !> \param[inout] psi1 第 1 态波函数向量 (n_pts)
    !> \param[inout] psi2 第 2 态波函数向量 (n_pts)
    !> \param[in] v11 第 1 态绝热势 V_11(R) (n_pts)
    !> \param[in] v22 第 2 态绝热势 V_22(R) (n_pts)
    !> \param[in] v12 非绝热耦合矩阵元 V_12(R) (n_pts)
    !> \param[in] dx 坐标空间网格步长 (a.u.)
    !> \param[in] mass 核约化质量 (a.u.)
    !> \param[in] dt 时间步长 (a.u.)
    subroutine propagate_split_operator_2channel(psi1, psi2, v11, v22, v12, dx, mass, dt)
        complex(dp), intent(inout) :: psi1(:), psi2(:)
        real(dp), intent(in) :: v11(:), v22(:), v12(:)
        real(dp), intent(in) :: dx, mass, dt

        integer :: n, i
        real(dp) :: dp_k, p_val, half_dt
        real(dp), allocatable :: p_kin(:)
        complex(dp), allocatable :: exp_t(:)

        n = size(psi1)
        half_dt = 0.5_dp * dt
        allocate(p_kin(n), exp_t(n))

        ! 动量空间相位算子: exp(-i * p^2/(2*m) * dt)
        dp_k = TWOPI / (real(n, dp) * dx)
        do i = 1, n
            if (i <= n / 2) then
                p_val = real(i - 1, dp) * dp_k
            else
                p_val = real(i - 1 - n, dp) * dp_k
            end if
            p_kin(i) = (p_val**2) / (2.0_dp * mass)
            exp_t(i) = exp(-EYE * p_kin(i) * dt)
        end do

        ! 1. 势能前半步: exp(-i * V * dt / 2) 作用于两通道
        call apply_2x2_potential_exponential(psi1, psi2, v11, v22, v12, half_dt)

        ! 2. 动量空间演化: 变换到动量空间 -> 相位乘积 -> 逆变换回坐标空间
        call fft_1d(psi1, -1)
        call fft_1d(psi2, -1)

        psi1 = psi1 * exp_t
        psi2 = psi2 * exp_t

        call fft_1d(psi1, 1)
        call fft_1d(psi2, 1)

        ! 3. 势能后半步: exp(-i * V * dt / 2)
        call apply_2x2_potential_exponential(psi1, psi2, v11, v22, v12, half_dt)

        deallocate(p_kin, exp_t)
    end subroutine propagate_split_operator_2channel

    !> \brief 解析计算 2x2 对称势能矩阵指数作用:
    !>        exp(-i*V*dt) = exp(-i*v_bar*dt) * [ cos(Omega*dt)*I -
    !>                       i*sin(Omega*dt)/Omega * (Delta*sigma_z + v12*sigma_x) ]
    subroutine apply_2x2_potential_exponential(psi1, psi2, v11, v22, v12, dt_step)
        complex(dp), intent(inout) :: psi1(:), psi2(:)
        real(dp), intent(in) :: v11(:), v22(:), v12(:)
        real(dp), intent(in) :: dt_step

        integer :: i, n
        real(dp) :: v_bar, delta_v, omega, c_val, s_val
        complex(dp) :: phase_common, u11, u22, u12, p1_old, p2_old

        n = size(psi1)
        do i = 1, n
            v_bar = 0.5_dp * (v11(i) + v22(i))
            delta_v = 0.5_dp * (v11(i) - v22(i))
            omega = sqrt(delta_v**2 + v12(i)**2)
            phase_common = exp(-EYE * v_bar * dt_step)

            if (omega < 1.0e-14_dp) then
                u11 = phase_common
                u22 = phase_common
                u12 = (0.0_dp, 0.0_dp)
            else
                c_val = cos(omega * dt_step)
                s_val = sin(omega * dt_step) / omega
                u11 = phase_common * cmplx(c_val, -delta_v * s_val, kind=dp)
                u22 = phase_common * cmplx(c_val, delta_v * s_val, kind=dp)
                u12 = phase_common * cmplx(0.0_dp, -v12(i) * s_val, kind=dp)
            end if

            p1_old = psi1(i)
            p2_old = psi2(i)

            psi1(i) = u11 * p1_old + u12 * p2_old
            psi2(i) = u12 * p1_old + u22 * p2_old
        end do
    end subroutine apply_2x2_potential_exponential

    !> \brief 计算 Landau-Zener 经典非绝热跃迁概率:
    !>        P_LZ = exp( - 2 * pi * V_12^2 / (v * |d(V1 - V2)/dR|) )
    pure function landau_zener_probability(v12_crossing, velocity, delta_slope) result(p_lz)
        real(dp), intent(in) :: v12_crossing, velocity, delta_slope
        real(dp) :: p_lz
        real(dp) :: denom

        denom = abs(velocity * delta_slope)
        if (denom <= 1.0e-15_dp) then
            p_lz = 0.0_dp
        else
            p_lz = exp(-TWOPI * (v12_crossing**2) / denom)
        end if
    end function landau_zener_probability

    !> \brief 计算通道 1 与通道 2 上的波包总几率与分支比
    pure subroutine calculate_channel_populations(psi1, psi2, dx, pop1, pop2, ratio)
        complex(dp), intent(in) :: psi1(:), psi2(:)
        real(dp), intent(in) :: dx
        real(dp), intent(out) :: pop1, pop2, ratio
        real(dp) :: sum1, sum2
        integer :: i, n

        sum1 = 0.0_dp
        sum2 = 0.0_dp
        n = size(psi1)

        do i = 1, n
            sum1 = sum1 + (abs(psi1(i))**2) * dx
            sum2 = sum2 + (abs(psi2(i))**2) * dx
        end do

        pop1 = sum1
        pop2 = sum2
        if (sum1 + sum2 > 1.0e-15_dp) then
            ratio = sum2 / (sum1 + sum2)
        else
            ratio = 0.0_dp
        end if
    end subroutine calculate_channel_populations

end module mod_multistate_coupling
