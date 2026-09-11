!> \brief 示例 3: 一维量子波包 Split-Operator 动力学演化与复吸收边界 (CAP)
!> \details 模拟高斯波包穿过势垒发生透射与反射，并在边界被复吸收势平滑吸收，演示概率通量与范数守恒。
program ex03_split_operator_1d
    use general_module
    implicit none

    integer, parameter :: nx = 512
    integer, parameter :: n_time_steps = 1000
    real(dp) :: x_min, x_max, dx, mass, dt, t
    real(dp) :: x_0, p_0, sigma_x, v_barrier, barrier_width
    real(dp) :: x_grid(nx), v_pot(nx)
    complex(dp) :: psi(nx)
    type(absorbing_boundary_t) :: cap
    real(dp) :: norm_tot, norm_inside, flux
    integer :: i, step, file_unit, snap_unit

    print '(A)', "=========================================================="
    print '(A)', "  Example 03: 1D Split-Operator Wavepacket with CAP       "
    print '(A)', "=========================================================="

    ! 1. 空间网格设置 [-40, 60] a.u.
    x_min = -40.0_dp
    x_max = 60.0_dp
    dx = (x_max - x_min) / real(nx, dp)
    mass = 1.0_dp   ! 电子质量 1 a.u.

    do i = 1, nx
        x_grid(i) = x_min + real(i - 1, dp) * dx
    end do

    ! 2. 势垒设置: 高斯势垒 V(x) = V_0 * exp(-(x / w)^2)
    v_barrier = 0.1_dp      ! 势垒高度 0.1 Hartree (~2.7 eV)
    barrier_width = 2.0_dp
    do i = 1, nx
        v_pot(i) = v_barrier * exp(-(x_grid(i) / barrier_width)**2)
    end do

    ! 3. 初始化复吸收边界 (在 x > 40 a.u. 设置吸收)
    call cap_init(r_start=40.0_dp, r_end=58.0_dp, strength=0.2_dp, &
                  cap_type=CAP_SIN2, cap_obj=cap)

    ! 4. 准备初始高斯波包 (从 x_0 = -20 向正方向发射，初始动量 p_0)
    x_0 = -20.0_dp
    p_0 = 0.5_dp        ! 动能 E_k = p_0^2 / 2m = 0.125 Hartree (略高于势垒)
    sigma_x = 2.5_dp

    do i = 1, nx
        psi(i) = (1.0_dp / (PI * sigma_x**2)**0.25_dp) * &
                 exp(-0.5_dp * ((x_grid(i) - x_0) / sigma_x)**2) * &
                 exp(EYE * p_0 * x_grid(i))
    end do

    ! 初始归一化校验
    norm_tot = sum(abs(psi)**2) * dx
    print '(A, F10.6)', "Initial wavepacket norm: ", norm_tot

    ! 5. 打开结果输出文件
    open(newunit=file_unit, file="wavepacket_dynamics.dat", status="replace", action="write")
    write(file_unit, '(A)') "# Time(au)  TotalNorm  NormInside(x<40)  FluxAtDetection"

    open(newunit=snap_unit, file="wavepacket_snapshots.dat", status="replace", action="write")
    write(snap_unit, '(A)') "# x(au)  |Psi_t0|^2  |Psi_tMid|^2  |Psi_tFinal|^2"

    ! 保存 t = 0 快照
    do i = 1, nx
        write(snap_unit, '(ES16.8)', advance='no') abs(psi(i))**2
    end do
    write(snap_unit, *)

    ! 6. 时间推进循环 (dt = 0.1 a.u.)
    dt = 0.1_dp
    do step = 1, n_time_steps
        t = real(step, dp) * dt

        ! A. 分裂算符单步推进
        call propagate_split_operator_1d(psi, v_pot, dx, mass, dt)

        ! B. 复吸收边界衰减
        call cap_apply_mask(x_grid, dt, cap, psi)

        ! C. 动力学可观测量监测 (每 10 步记录一次)
        if (mod(step, 10) == 0) then
            norm_tot = sum(abs(psi)**2) * dx
            norm_inside = calculate_norm_inside(x_grid, psi, cap%r_start)
            flux = calculate_probability_flux(x_grid, mass, psi, nx / 2)
            write(file_unit, '(4ES16.8)') t, norm_tot, norm_inside, flux
        end if
    end do

    close(file_unit)
    close(snap_unit)
    print '(A, F10.6)', "Final wavepacket norm inside [x<40]: ", norm_inside
    print '(A)', "Dynamics output saved to: wavepacket_dynamics.dat"

end program ex03_split_operator_1d
