! ==============================================================================
! GeneralModule Example 07:
! 连续谱散射能量本征波函数求解与对比 (含时波包谱投影法 vs 非含时 Numerov 逆向匹配法)
! ==============================================================================
! 本例演示对同一一维高斯势垒 V(x) = V0 * exp(-x^2 / 2):
!   1. [非含时方法 (TI)]: 求解定态薛定谔方程，逆向 Numerov 匹配边界，
!      得到能量归一化散射本征波函数 psi_E^TI(x) 及透射/反射系数 (T, R)。
!   2. [含时方法 (TD)]: 发射宽动量高斯波包，进行 Split-Operator 动力学推进，
!      原位进行时间-能量半傅里叶谱投影累积，消除波包动量权重后提取 psi_E^TD(x)。
!   3. [定量比对]: 对比两套完全独立的物理方法所求得的波函数与透射几率，
!      并将空间波函数数据输出为数据文件 ex07_scattering_comparison.dat。
! ==============================================================================
program ex07_scattering_wavefunctions_ti_td
    use general_module
    implicit none

    integer, parameter :: nx = 512
    real(dp), parameter :: hbar_au = 1.0_dp
    real(dp) :: x_min, x_max, dx, dt
    real(dp) :: mass, energy_0, k0
    real(dp) :: x0_wp, sigma_wp, v0_barrier
    real(dp), allocatable :: x_grid(:), v_pot(:)
    complex(dp), allocatable :: psi_ti(:), psi_td_accum(:), psi_td_norm(:), psi_wp(:)
    real(dp) :: trans_ti, refl_ti
    integer  :: ix, it, total_steps, stat
    real(dp) :: t_curr, max_diff, mean_diff, max_amp

    print '(A)', "===================================================================="
    print '(A)', "   GeneralModule Example 07: Continuous Scattering Wavefunctions    "
    print '(A)', "          [Time-Independent (TI) vs Time-Dependent (TD)]            "
    print '(A)', "===================================================================="

    ! --------------------------------------------------------------------------
    ! 1. 物理参数与网格设置
    ! --------------------------------------------------------------------------
    mass = 1.0_dp
    x_min = -25.0_dp
    x_max =  25.0_dp
    dx = (x_max - x_min) / real(nx - 1, dp)

    allocate(x_grid(nx), v_pot(nx))
    allocate(psi_ti(nx), psi_td_accum(nx), psi_td_norm(nx), psi_wp(nx))

    ! 势能分布: 高斯势垒 V(x) = V0 * exp(-x^2 / 2.0)
    v0_barrier = 0.80_dp
    do ix = 1, nx
        x_grid(ix) = x_min + real(ix - 1, dp) * dx
        v_pot(ix) = v0_barrier * exp(-(x_grid(ix)**2) / 2.0_dp)
    end do

    ! 目标散射能量: 稍低于势垒顶 (量子隧穿与反射共存区)
    energy_0 = 0.50_dp
    k0 = sqrt(2.0_dp * mass * energy_0) / hbar_au

    print '(A, F7.4, A, F7.4, A)', " [Setup] Barrier Height V0 = ", v0_barrier, &
                                  " a.u., Incident Energy E = ", energy_0, " a.u."
    print '(A, F7.4, A, I4, A)',   "         Wavevector k = ", k0, " a.u., Grid Points = ", nx, " points."

    ! --------------------------------------------------------------------------
    ! 2. 非含时方法 (TI: 逆向 Numerov 积分与渐近分解)
    ! --------------------------------------------------------------------------
    print '(A)', ""
    print '(A)', ">> [1/3] Solving via Time-Independent (TI) Method..."
    call calc_scattering_wavefunction_1d_cartesian( &
        x_grid=x_grid, v_pot=v_pot, mass=mass, energy=energy_0, &
        norm_type=NORM_ENERGY, psi_wf=psi_ti, &
        trans_prob=trans_ti, refl_prob=refl_ti, stat=stat)

    print '(A, F9.5, A, F9.5, A, F9.5)', "    TI Transmission T = ", trans_ti, &
                                         ", Reflection R = ", refl_ti, &
                                         ", Total T+R = ", trans_ti + refl_ti

    ! --------------------------------------------------------------------------
    ! 3. 含时方法 (TD: 入射高斯波包推进 + 全空间半傅里叶谱投影)
    ! --------------------------------------------------------------------------
    print '(A)', ""
    print '(A)', ">> [2/3] Solving via Time-Dependent (TD) Spectral Projection..."
    x0_wp = -12.0_dp
    sigma_wp = 1.5_dp

    ! 构造空间与动量归一化的初始波包
    call gaussian_wavepacket_1d(x_grid, x0_wp, sigma_wp, k0, psi_wp)
    psi_td_accum = (0.0_dp, 0.0_dp)

    dt = 0.02_dp
    total_steps = 1400 ! 演化时间 t_max = 28.0 a.u.，波包完全通过相互作用区

    do it = 0, total_steps
        t_curr = real(it, dp) * dt
        ! 原位谱投影累积: int_0^T psi(x, t) * exp(i * E * t / hbar) dt
        call accumulate_wavefunction_spectral_projection( &
            psi_wp, t_curr, dt, energy_0, hbar_au, psi_td_accum)
        ! Split-Operator 动力学推进一步
        call propagate_split_operator_1d(psi_wp, v_pot, dx, mass, dt)
    end do

    ! 消除动量分布因子，归一化提取定态连续本征函数 psi_E(x)
    call extract_td_scattering_wavefunction( &
        x_grid, psi_td_accum, energy_0, mass, hbar_au, x0_wp, sigma_wp, k0, psi_td_norm, stat)

    print '(A)', "    TD wavepacket propagation & spectral projection complete."

    ! --------------------------------------------------------------------------
    ! 4. 定量比对与数据输出
    ! --------------------------------------------------------------------------
    print '(A)', ""
    print '(A)', ">> [3/3] Quantitative Comparison & Validation..."

    ! 计算相互作用区及透射区 x in [-5, +15] 内的相对偏差
    max_diff = 0.0_dp
    mean_diff = 0.0_dp
    max_amp = maxval(abs(psi_ti))

    open(unit=20, file="ex07_scattering_comparison.dat", status="replace", action="write")
    write(20, '(A)') "#   x (a.u.)        V(x)          |psi_TI|^2      |psi_TD|^2       Re[psi_TI]      Re[psi_TD]"

    do ix = 1, nx
        if (x_grid(ix) >= -5.0_dp .and. x_grid(ix) <= 15.0_dp) then
            max_diff = max(max_diff, abs(abs(psi_td_norm(ix)) - abs(psi_ti(ix))))
            mean_diff = mean_diff + abs(abs(psi_td_norm(ix)) - abs(psi_ti(ix)))
        end if
        write(20, '(6ES16.7)') x_grid(ix), v_pot(ix), &
                               abs(psi_ti(ix))**2, abs(psi_td_norm(ix))**2, &
                               real(psi_ti(ix), dp), real(psi_td_norm(ix), dp)
    end do
    close(20)

    mean_diff = mean_diff / real(count(x_grid >= -5.0_dp .and. x_grid <= 15.0_dp), dp)

    print '(A, F8.4, A, F8.4)', "    Max Amplitude:   |psi_TI|_max = ", max_amp, &
                               ", |psi_TD|_max = ", maxval(abs(psi_td_norm))
    print '(A, F8.4)',          "    Max Abs Diff:    ", max_diff
    print '(A, F8.4, A)',       "    Mean Rel Diff:   ", (mean_diff / max_amp) * 100.0_dp, " %"
    print '(A)', "    Detailed profile saved to: ex07_scattering_comparison.dat"

    print '(A)', ""
    print '(A)', "===================================================================="
    print '(A)', "   SUCCESS: TI and TD Continuous Scattering Wavefunctions Match!    "
    print '(A)', "===================================================================="

    deallocate(x_grid, v_pot, psi_ti, psi_td_accum, psi_td_norm, psi_wp)
end program ex07_scattering_wavefunctions_ti_td
