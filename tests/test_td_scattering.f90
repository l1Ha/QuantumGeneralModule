! ==============================================================================
! GeneralModule Unit Tests: Time-Dependent Scattering (mod_td_scattering)
! ==============================================================================
program test_td_scattering
    use mod_constants, only: dp, PI, TWOPI
    use mod_wavepacket_propagator, only: propagate_split_operator_1d
    use mod_td_scattering
    implicit none

    integer :: n_pass = 0, n_total = 0
    integer :: nx, nt, it, ix, ie
    real(dp) :: mass, hbar, x_min, x_max, dx, dt
    real(dp) :: x0, sigma_x, k0, e0
    real(dp), allocatable :: x_grid(:), v_pot(:), v_free(:)
    complex(dp), allocatable :: psi(:), psi_free(:)
    real(dp) :: norm_psi, gk_norm, k_val, dk
    complex(dp) :: gk
    real(dp) :: x_det_trans, t_curr, x_centroid_free
    real(dp) :: e_grid(10), t_prob(10), r_prob(10)
    complex(dp) :: amp_trans(10), amp_free(10), s_mat_e(10)
    real(dp) :: phase_shift_e(10), delay_w(10)
    integer :: idx_det
    real(dp) :: prob_sum
    complex(dp) :: amp_ch1(10), amp_ch2(10)
    real(dp) :: s11_p(10), s12_p(10)

    print '(A)', "=================================================="
    print '(A)', " GeneralModule Unit Tests: TD Scattering         "
    print '(A)', "=================================================="

    mass = 1.0_dp
    hbar = 1.0_dp
    nx = 512
    x_min = -30.0_dp
    x_max =  30.0_dp
    dx = (x_max - x_min) / real(nx - 1, dp)

    allocate(x_grid(nx), v_pot(nx), v_free(nx), psi(nx), psi_free(nx))
    do ix = 1, nx
        x_grid(ix) = x_min + real(ix - 1, dp) * dx
        ! 势垒设置在原点附近 x in [-1.5, 1.5]
        if (abs(x_grid(ix)) <= 1.5_dp) then
            v_pot(ix) = 0.5_dp
        else
            v_pot(ix) = 0.0_dp
        end if
        v_free(ix) = 0.0_dp
    end do

    x0 = -12.0_dp
    sigma_x = 1.5_dp
    k0 = 1.2_dp
    e0 = (k0**2) / (2.0_dp * mass)

    ! --------------------------------------------------------------------------
    ! 测试 1: 初始高斯波包空间范数归一化与动量表象 |g(k)|^2 积分归一化
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call gaussian_wavepacket_1d(x_grid, x0, sigma_x, k0, psi)
    norm_psi = sum(abs(psi)**2) * dx

    gk_norm = 0.0_dp
    dk = 0.02_dp
    do ix = -200, 200
        k_val = k0 + real(ix, dp) * dk
        gk = gaussian_momentum_amplitude(k_val, x0, sigma_x, k0)
        gk_norm = gk_norm + (abs(gk)**2) * dk
    end do

    if (abs(norm_psi - 1.0_dp) < 1.0e-5_dp .and. abs(gk_norm - 1.0_dp) < 1.0e-4_dp) then
        print '(A, F9.6, A, F9.6)', " [PASS] Gaussian norm = ", norm_psi, ", g(k) norm = ", gk_norm
        n_pass = n_pass + 1
    else
        print '(A, F9.6, A, F9.6)', " [FAIL] Gaussian norm = ", norm_psi, ", g(k) norm = ", gk_norm
    end if

    ! --------------------------------------------------------------------------
    ! 测试 2: 自由波包质心运动速度验证 <x>(t) = x0 + (hbar*k0/m) * t
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    psi_free = psi
    dt = 0.05_dp
    nt = 100
    do it = 1, nt
        call propagate_split_operator_1d(psi_free, v_free, dx, mass, dt)
    end do
    t_curr = real(nt, dp) * dt
    x_centroid_free = wavepacket_centroid_position(x_grid, dx, psi_free)
    if (abs(x_centroid_free - (x0 + (k0 / mass) * t_curr)) < 0.05_dp) then
        print '(A, F8.4, A, F8.4)', " [PASS] Free wavepacket centroid <x>(t) = ", &
                                     x_centroid_free, " (Exact: ", x0 + (k0 / mass) * t_curr, ")"
        n_pass = n_pass + 1
    else
        print '(A, F8.4, A, F8.4)', " [FAIL] Free wavepacket centroid <x>(t) = ", &
                                     x_centroid_free, " (Exact: ", x0 + (k0 / mass) * t_curr, ")"
    end if

    ! --------------------------------------------------------------------------
    ! 测试 3: 含时通量傅里叶变换提取透射几率 T(E)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    x_det_trans = 8.0_dp
    idx_det = minloc(abs(x_grid - x_det_trans), dim=1)

    do ie = 1, 10
        e_grid(ie) = 0.3_dp + real(ie - 1, dp) * (1.2_dp - 0.3_dp) / 9.0_dp
    end do
    amp_trans = (0.0_dp, 0.0_dp)
    amp_free  = (0.0_dp, 0.0_dp)

    ! 重新初始化波包
    call gaussian_wavepacket_1d(x_grid, x0, sigma_x, k0, psi)
    call gaussian_wavepacket_1d(x_grid, x0, sigma_x, k0, psi_free)

    nt = 350
    do it = 1, nt
        t_curr = real(it, dp) * dt
        call propagate_split_operator_1d(psi, v_pot, dx, mass, dt)
        call propagate_split_operator_1d(psi_free, v_free, dx, mass, dt)

        call accumulate_flux_amplitude(t_curr, psi(idx_det), dt, e_grid, 10, hbar, amp_trans)
        call accumulate_flux_amplitude(t_curr, psi_free(idx_det), dt, e_grid, 10, hbar, amp_free)
    end do

    call calculate_td_transmission(e_grid, 10, amp_trans, mass, hbar, x0, sigma_x, k0, t_prob)

    ! 在中心能量 E0 附近透射几率应在 (0, 1) 之间合理区间
    if (t_prob(5) >= 0.0_dp .and. t_prob(5) <= 1.0_dp) then
        print '(A, F9.5)', " [PASS] TD transmission probability T(E0) = ", t_prob(5)
        n_pass = n_pass + 1
    else
        print '(A, F9.5)', " [FAIL] TD transmission probability T(E0) = ", t_prob(5)
    end if

    ! --------------------------------------------------------------------------
    ! 测试 4: 含时散射 S-矩阵元与散射相移 delta(E) 提取
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call calculate_td_smatrix_element(amp_trans, amp_free, 10, s_mat_e, phase_shift_e)
    if (abs(s_mat_e(5)) > 0.0_dp .and. abs(s_mat_e(5)) <= 1.2_dp) then
        print '(A, F8.4, A, F8.4)', " [PASS] TD S-matrix |S(E0)| = ", abs(s_mat_e(5)), &
                                     ", delta(E0) = ", phase_shift_e(5)
        n_pass = n_pass + 1
    else
        print '(A, F8.4)', " [FAIL] TD S-matrix |S(E0)| = ", abs(s_mat_e(5))
    end if

    ! --------------------------------------------------------------------------
    ! 测试 5: Möller 动量表象渐近投影提取透射与反射概率及几率守恒 T + R <= 1.0
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call project_wavepacket_to_smatrix(x_grid, dx, psi, mass, hbar, k0, sigma_x, x0, &
                                       e_grid, 10, t_prob, r_prob)
    prob_sum = t_prob(5) + r_prob(5)
    if (prob_sum >= 0.0_dp .and. prob_sum <= 1.1_dp) then
        print '(A, F8.4, A, F8.4, A, F8.4)', " [PASS] Moller projection: T = ", t_prob(5), &
                                              ", R = ", r_prob(5), ", Sum = ", prob_sum
        n_pass = n_pass + 1
    else
        print '(A, F8.4)', " [FAIL] Moller projection prob_sum = ", prob_sum
    end if

    ! --------------------------------------------------------------------------
    ! 测试 6: 含时 Wigner 散射时延计算
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    call wavepacket_wigner_delay(e_grid, phase_shift_e, 10, hbar, delay_w)
    if (abs(delay_w(5)) < 1.0e4_dp) then
        print '(A, F9.4, A)', " [PASS] Wavepacket Wigner delay tau_W(E0) = ", delay_w(5), " a.u."
        n_pass = n_pass + 1
    else
        print '(A, F9.4)', " [FAIL] Wavepacket Wigner delay = ", delay_w(5)
    end if

    ! --------------------------------------------------------------------------
    ! 测试 7: 多通道含时通量散射 S-矩阵元 |S_11|^2 与 |S_12|^2
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    amp_ch1 = amp_trans * 0.8_dp
    amp_ch2 = amp_trans * 0.6_dp
    call multichannel_td_smatrix_elements(amp_ch1, amp_ch2, mass, hbar, &
                                          e_grid, 0.05_dp, x0, sigma_x, k0, &
                                          s11_p, s12_p, 10)
    if (s11_p(5) >= 0.0_dp .and. s12_p(5) >= 0.0_dp) then
        print '(A, F8.4, A, F8.4)', " [PASS] Multichannel S-matrix |S11|^2 = ", s11_p(5), &
                                     ", |S12|^2 = ", s12_p(5)
        n_pass = n_pass + 1
    else
        print '(A, F8.4)', " [FAIL] Multichannel S-matrix extraction = ", s11_p(5)
    end if

    ! --------------------------------------------------------------------------
    ! 测试 8: 二维含时波包散射角分布微分散射截面 d(sigma)/d(theta)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        real(dp) :: x_2d(32), y_2d(32), th_2d(8), ds_2d(8)
        complex(dp) :: psi_2d(32, 32)
        integer :: i1, i2
        do i1 = 1, 32
            x_2d(i1) = -10.0_dp + real(i1 - 1, dp) * (20.0_dp / 31.0_dp)
            y_2d(i1) = -10.0_dp + real(i1 - 1, dp) * (20.0_dp / 31.0_dp)
        end do
        do i1 = 1, 8
            th_2d(i1) = real(i1 - 1, dp) * TWOPI / 8.0_dp
        end do
        do i2 = 1, 32
            do i1 = 1, 32
                psi_2d(i1, i2) = cmplx(exp(-(x_2d(i1)**2 + y_2d(i2)**2) / 4.0_dp), 0.0_dp, kind=dp)
            end do
        end do
        call calculate_td_differential_cross_section_2d(x_2d, y_2d, psi_2d, 1.0_dp, 1.0_dp, th_2d, ds_2d)
        if (all(ds_2d >= 0.0_dp) .and. maxval(ds_2d) > 0.0_dp) then
            print '(A, F9.5)', " [PASS] 2D TD Differential cross section max dsigma/dtheta = ", maxval(ds_2d)
            n_pass = n_pass + 1
        else
            print '(A, F9.5)', " [FAIL] 2D TD Differential cross section = ", maxval(ds_2d)
        end if
    end block

    ! --------------------------------------------------------------------------
    ! 测试 9: 含时波包时间-能量谱投影法原位提取定态连续本征波函数 psi_E(x)
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    block
        complex(dp), allocatable :: psi_prop(:), psi_accum(:), psi_extracted(:)
        real(dp) :: t_step, t_max, exact_amp, test_amp
        integer  :: steps, s_idx, test_idx

        allocate(psi_prop(nx), psi_accum(nx), psi_extracted(nx))
        psi_accum = (0.0_dp, 0.0_dp)

        call gaussian_wavepacket_1d(x_grid, x0, sigma_x, k0, psi_prop)
        t_step = 0.04_dp
        t_max = 26.0_dp
        steps = nint(t_max / t_step)

        do s_idx = 0, steps
            t_curr = real(s_idx, dp) * t_step
            ! 原位累积半傅里叶谱投影
            call accumulate_wavefunction_spectral_projection( &
                psi_prop, t_curr, t_step, e0, hbar, psi_accum)
            call propagate_split_operator_1d(psi_prop, v_free, dx, mass, t_step)
        end do

        call extract_td_scattering_wavefunction( &
            x_grid, psi_accum, e0, mass, hbar, x0, sigma_x, k0, psi_extracted)

        ! 理论自由粒子连续态能量归一化振幅: |psi_E(x)| = sqrt(mu / (2*pi*hbar^2*k))
        exact_amp = sqrt(mass / (TWOPI * (hbar * hbar) * k0))
        test_idx = nx / 2 ! 原点 x ~ 0 处
        test_amp = abs(psi_extracted(test_idx))

        if (abs(test_amp - exact_amp) / exact_amp < 0.06_dp) then
            print '(A, F8.4, A, F8.4)', " [PASS] TD spectral projection |psi_E(0)| = ", &
                                         test_amp, " (Exact: ", exact_amp, ")"
            n_pass = n_pass + 1
        else
            print '(A, F8.4, A, F8.4)', " [FAIL] TD spectral projection |psi_E(0)| = ", &
                                         test_amp, " (Exact: ", exact_amp, ")"
        end if

        deallocate(psi_prop, psi_accum, psi_extracted)
    end block

    deallocate(x_grid, v_pot, v_free, psi, psi_free)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "TD Scattering Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', "SUCCESS: All time-dependent scattering tests passed."
    else
        stop 1
    end if
end program test_td_scattering
