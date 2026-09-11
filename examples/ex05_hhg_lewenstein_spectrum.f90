!> \brief 示例 5: 强场原子高次谐波 (HHG) 发射谱与半经典截止能模拟
!> \details 模拟氩原子在 800 nm 强激光场中的阿秒高次谐波发射，计算 Ehrenfest 偶极加速度并提取 HHG 平台区与截止能。
program ex05_hhg_lewenstein_spectrum
    use general_module
    implicit none

    type(atom_config_t) :: atom
    type(pulse_config_t) :: laser
    integer, parameter :: nt = 4096
    real(dp) :: t_min, t_max, dt, t
    real(dp) :: t_arr(nt), a_field(nt), e_field(nt), dip_acc(nt)
    real(dp) :: omega_arr(nt/2), hhg_spec(nt/2)
    real(dp) :: ip_au, up_val, e_cutoff, order_cut
    integer :: i, file_unit

    print '(A)', "=========================================================="
    print '(A)', "  Example 05: High Harmonic Generation (HHG) Simulation   "
    print '(A)', "=========================================================="

    ! 1. 选取目标靶原子与强场参数 (Argon 原子, 800 nm 激光场)
    call get_atom_config("Ar", atom)
    ip_au = atom%ip_au

    laser%shape_type = PULSE_SIN2
    laser%field_peak = 0.053_dp           ! 峰值场强 ~1.0e14 W/cm^2
    laser%freq_central = 0.057_dp         ! 800 nm 基频
    laser%duration = 30.0_dp * FS2AU      ! 30 fs 脉冲宽度
    laser%t_center = 0.0_dp
    laser%cep_phase = 0.0_dp

    ! 2. 强场理论参数分析
    up_val = ponderomotive_energy(laser%field_peak, laser%freq_central)
    e_cutoff = hhg_cutoff_energy(ip_au, laser%field_peak, laser%freq_central)
    order_cut = e_cutoff / laser%freq_central

    print '(A, A)',          " Target Atom:            ", trim(atom%name)
    print '(A, F10.4, A)',   " Ionization Potential:   ", ip_au * AU2EV, " eV"
    print '(A, F10.4, A)',   " Laser Frequency (800nm):", laser%freq_central * AU2EV, " eV"
    print '(A, F10.4, A)',   " Ponderomotive Energy Up:", up_val * AU2EV, " eV"
    print '(A, F10.4, A)',   " Keldysh Parameter:      ", keldysh_parameter(laser%freq_central, laser%field_peak, ip_au), ""
    print '(A, F10.2, A, F6.1)', " Theoretical HHG Cutoff: ", e_cutoff * AU2EV, " eV (Harmonic Order ~", order_cut, ")"
    print '(A)', "----------------------------------------------------------"

    ! 3. 时间网格与外场序列采样
    t_min = -laser%duration
    t_max = laser%duration
    dt = (t_max - t_min) / real(nt, dp)

    do i = 1, nt
        t = t_min + real(i - 1, dp) * dt
        t_arr(i) = t
        e_field(i) = pulse_electric_field(t, laser)
        a_field(i) = pulse_vector_potential(t, laser)
    end do

    ! 4. 计算强场 SFA 鞍点偶极发射 (Lewenstein 模型)
    print '(A)', "Calculating Lewenstein SFA dipole moment time-series..."
    do i = 1, nt
        t = t_arr(i)
        ! 使用过去窗口点计算瞬时偶极
        dip_acc(i) = lewenstein_sfa_dipole(t, t_arr(1:i), a_field(1:i), e_field(1:i), ip_au, dt)
    end do

    ! 5. 傅里叶变换提取 HHG 谐波功率发射谱
    call hhg_power_spectrum(dip_acc, dt, omega_arr, hhg_spec)

    ! 6. 保存数据至数据文件
    open(newunit=file_unit, file="hhg_spectrum.dat", status="replace", action="write")
    write(file_unit, '(A)') "# Harmonic_Order  Energy(eV)  Power_Spectrum"
    do i = 1, min(nt / 2, 120)
        write(file_unit, '(3ES16.8)') omega_arr(i) / laser%freq_central, &
                                      omega_arr(i) * AU2EV, &
                                      hhg_spec(i)
    end do
    close(file_unit)

    print '(A)', "HHG Spectrum saved to: hhg_spectrum.dat"
    print '(A)', "=========================================================="

end program ex05_hhg_lewenstein_spectrum
