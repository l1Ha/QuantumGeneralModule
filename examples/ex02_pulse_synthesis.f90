!> \brief 示例 2: 超快强激光脉冲时域合成与频域谱分析
!> \details 演示线性啁啾脉冲、双色不对称场（omega + 2*omega）及太赫兹脉冲串的生成与 FFT 频谱转换。
program ex02_pulse_synthesis
    use general_module
    implicit none

    integer, parameter :: nt = 2048
    real(dp) :: t_min, t_max, dt, t_val
    real(dp) :: t_arr(nt)
    complex(dp) :: e_chirp(nt), e_twocolor(nt), e_thz(nt)
    type(pulse_config_t) :: cfg_chirp, cfg_twocolor, cfg_thz
    integer :: i, file_unit

    print '(A)', "=========================================================="
    print '(A)', "  Example 02: Ultrafast Laser Pulse Synthesis & FFT      "
    print '(A)', "=========================================================="

    ! 1. 时间网格设置 (-200 fs 到 +200 fs)
    t_min = -200.0_dp * FS2AU
    t_max = 200.0_dp * FS2AU
    dt = (t_max - t_min) / real(nt, dp)

    do i = 1, nt
        t_arr(i) = t_min + real(i - 1, dp) * dt
    end do

    ! 2. 配置线性啁啾激光脉冲 (Chirped Pulse)
    cfg_chirp%shape_type = PULSE_CHIRP
    cfg_chirp%field_peak = 0.03_dp           ! 峰值场强 ~3.1e13 W/cm^2
    cfg_chirp%freq_central = 0.057_dp        ! 800 nm 基频
    cfg_chirp%duration = 40.0_dp * FS2AU     ! 40 fs FWHM
    cfg_chirp%t_center = 0.0_dp
    cfg_chirp%chirp_rate = 1.0e-5_dp         ! 正啁啾系数
    cfg_chirp%cep_phase = 0.0_dp

    ! 3. 配置双色不对称合成场 (Two-color omega + 2omega)
    cfg_twocolor%shape_type = PULSE_TWOCOLOR
    cfg_twocolor%field_peak = 0.03_dp
    cfg_twocolor%freq_central = 0.057_dp
    cfg_twocolor%duration = 50.0_dp * FS2AU
    cfg_twocolor%t_center = 0.0_dp
    cfg_twocolor%two_color_ratio = 0.3_dp    ! 2*omega 强度为 30%
    cfg_twocolor%two_color_phase = HALFPI    ! 相对相位 pi/2 产生最大空间不对称度

    ! 4. 配置太赫兹单周期脉冲串 (THz Train)
    cfg_thz%shape_type = PULSE_THZ_TRAIN
    cfg_thz%field_peak = 0.005_dp
    cfg_thz%freq_central = 0.002_dp          ! 低频 THz (~13 THz)
    cfg_thz%duration = 80.0_dp * FS2AU
    cfg_thz%t_center = -100.0_dp * FS2AU
    cfg_thz%train_count = 3                  ! 3 个脉冲序列
    cfg_thz%train_delay = 100.0_dp * FS2AU   ! 间隔 100 fs

    ! 5. 采样电场时域信号
    do i = 1, nt
        t_val = t_arr(i)
        e_chirp(i) = cmplx(pulse_electric_field(t_val, cfg_chirp), 0.0_dp, kind=dp)
        e_twocolor(i) = cmplx(pulse_electric_field(t_val, cfg_twocolor), 0.0_dp, kind=dp)
        e_thz(i) = cmplx(pulse_electric_field(t_val, cfg_thz), 0.0_dp, kind=dp)
    end do

    ! 6. 导出时域波形数据
    open(newunit=file_unit, file="synthesized_pulses.dat", status="replace", action="write")
    write(file_unit, '(A)') "# Time(fs)  E_Chirp(au)  E_TwoColor(au)  E_THz(au)"
    do i = 1, nt
        write(file_unit, '(4ES16.8)') t_arr(i) * AU2FS, real(e_chirp(i), dp), &
            real(e_twocolor(i), dp), real(e_thz(i), dp)
    end do
    close(file_unit)
    print '(A)', "Time-domain pulses saved to: synthesized_pulses.dat"

    ! 7. 频谱变换展示 (调用内置 1D FFT)
    call fft_1d(e_chirp, -1)
    open(newunit=file_unit, file="spectrum_chirp.dat", status="replace", action="write")
    write(file_unit, '(A)') "# Frequency(au)  Power_Spectrum(|E|^2)"
    do i = 1, nt / 2
        write(file_unit, '(2ES16.8)') real(i - 1, dp) * TWOPI / (real(nt, dp) * dt), &
            abs(e_chirp(i))**2
    end do
    close(file_unit)
    print '(A)', "FFT spectrum saved to: spectrum_chirp.dat"

end program ex02_pulse_synthesis
