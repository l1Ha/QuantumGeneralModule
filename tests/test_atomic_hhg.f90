!> \brief 强场原子物理、高次谐波与多态动力学扩展模块单元测试
program test_atomic_hhg
    use general_module
    implicit none

    integer :: n_tests = 0
    integer :: n_passed = 0
    real(dp) :: tol = 1.0e-5_dp

    ! 强场原子参数测试
    type(atom_config_t) :: atom_h, atom_ar
    real(dp) :: gamma_k, up_val, e_cut, v_pot, dv_pot, w_adk, w_adk_strong

    ! 激光场与矢量势测试
    type(pulse_config_t) :: p_cfg
    real(dp) :: ex, ey, a_vec

    ! 二维波包演化测试
    integer, parameter :: nx2 = 32, ny2 = 32
    complex(dp) :: psi_2d(nx2, ny2)
    real(dp) :: v_pot_2d(nx2, ny2)
    real(dp) :: dx2, dy2, norm_2d_init, norm_2d_final
    integer :: ix, iy, step

    ! 多态耦合与 Landau-Zener 测试
    real(dp) :: p_lz, pop1, pop2, ratio
    complex(dp) :: psi1(64), psi2(64)
    real(dp) :: v11(64), v22(64), v12(64)

    print '(A)', "=================================================="
    print '(A)', "  GeneralModule Extended Tests: Atomic & Dynamics "
    print '(A)', "=================================================="

    ! 1. 强场原子参数与截断定律校验
    call get_atom_config("H", atom_h)
    call assert_close("Hydrogen Ip = 0.5 a.u.", 0.5_dp, atom_h%ip_au, 1.0e-8_dp)

    call get_atom_config("Ar", atom_ar)
    call assert_close("Argon Ip = 0.579 a.u.", 0.57915_dp, atom_ar%ip_au, 1.0e-4_dp)

    ! Keldysh 参数: omega=0.057 (~800nm), E_0=0.05 -> gamma = 0.057 * 1.0 / 0.05 = 1.14
    gamma_k = keldysh_parameter(0.057_dp, 0.05_dp, 0.5_dp)
    call assert_close("Keldysh gamma = 1.14", 1.14_dp, gamma_k, 1.0e-4_dp)

    ! Ponderomotive energy: Up = E0^2 / (4*w^2) = 0.0025 / (4 * 0.057^2) = 0.19236 Hartree
    up_val = ponderomotive_energy(0.05_dp, 0.057_dp)
    call assert_close("Ponderomotive Up", 0.19236688_dp, up_val, 1.0e-4_dp)

    ! HHG 截断法则: Ip + 3.17 * Up = 0.5 + 3.1725955 * 0.19236688 = 1.1103 Hartree
    e_cut = hhg_cutoff_energy(0.5_dp, 0.05_dp, 0.057_dp)
    call assert_close("HHG Cutoff Law Ip + 3.17Up", 1.110306_dp, e_cut, 1.0e-4_dp)

    ! 2. 软核库仑势及其一阶导数
    v_pot = soft_core_coulomb_potential(0.0_dp, 1.0_dp, 1.0_dp)
    call assert_close("Soft-core V(0) = -1.0", -1.0_dp, v_pot, 1.0e-8_dp)

    dv_pot = soft_core_coulomb_derivative(0.0_dp, 1.0_dp, 1.0_dp)
    call assert_close("Soft-core dV/dx(0) = 0.0", 0.0_dp, dv_pot, 1.0e-8_dp)

    ! 3. ADK 隧穿电离率单调性
    w_adk = adk_ionization_rate(0.03_dp, 0.5_dp)
    w_adk_strong = adk_ionization_rate(0.06_dp, 0.5_dp)
    if (w_adk > 0.0_dp .and. w_adk_strong > w_adk) then
        call assert_close("ADK rate positive & monotonic", 1.0_dp, 1.0_dp, tol)
    else
        call assert_close("ADK rate positive & monotonic", 1.0_dp, 0.0_dp, tol)
    end if

    ! 4. 激光椭圆偏振与矢量势
    p_cfg%shape_type = PULSE_GAUSSIAN
    p_cfg%field_peak = 0.05_dp
    p_cfg%freq_central = 0.057_dp
    p_cfg%duration = 50.0_dp * FS2AU
    p_cfg%ellipticity = 0.5_dp
    p_cfg%t_center = 0.0_dp

    call pulse_electric_field_2d(0.0_dp, p_cfg, ex, ey)
    ! t=0 时 cos(0)=1, sin(0)=0 -> Ex > 0, Ey = 0
    call assert_close("Elliptical field Ey at t=0", 0.0_dp, ey, 1.0e-8_dp)

    ! 矢量势 A(t) 在 t=0 (CEP=0) 时 sin(0)=0 -> A(0) = 0
    a_vec = pulse_vector_potential(0.0_dp, p_cfg)
    call assert_close("Vector potential A(0)=0", 0.0_dp, a_vec, 1.0e-8_dp)

    ! 5. Landau-Zener 跃迁几率
    ! P_LZ = exp(-2*pi*V12^2 / (v * |dF|))
    p_lz = landau_zener_probability(0.1_dp, 1.0_dp, 1.0_dp)
    call assert_close("Landau-Zener probability", exp(-TWOPI * 0.01_dp), p_lz, 1.0e-6_dp)

    ! 6. 双态耦合布居数守恒与演化
    do ix = 1, 64
        psi1(ix) = cmplx(exp(-0.5_dp * (real(ix - 32, dp) / 3.0_dp)**2), 0.0_dp, kind=dp)
        psi2(ix) = (0.0_dp, 0.0_dp)
        v11(ix) = 0.1_dp
        v22(ix) = 0.2_dp
        v12(ix) = 0.05_dp
    end do
    call calculate_channel_populations(psi1, psi2, 0.2_dp, pop1, pop2, ratio)
    call assert_close("Initial State 2 population = 0", 0.0_dp, pop2, 1.0e-6_dp)

    ! 演化 10 步验证总几率守恒
    do step = 1, 10
        call propagate_split_operator_2channel(psi1, psi2, v11, v22, v12, 0.2_dp, 1.0_dp, 0.05_dp)
    end do
    call assert_close("2-Channel total population conservation", 1.0_dp, &
                      (pop1 + pop2) / (sum(abs(psi1)**2 + abs(psi2)**2) * 0.2_dp), 1.0e-6_dp)

    ! 7. 二维 Split-Operator 波包范数守恒
    dx2 = 0.5_dp
    dy2 = 0.5_dp
    do iy = 1, ny2
        do ix = 1, nx2
            v_pot_2d(ix, iy) = 0.02_dp * (real(ix - nx2/2, dp)**2 + real(iy - ny2/2, dp)**2)
            psi_2d(ix, iy) = cmplx(exp(-0.2_dp * (real(ix - nx2/2, dp)**2 + real(iy - ny2/2, dp)**2)), 0.0_dp, kind=dp)
        end do
    end do
    norm_2d_init = sum(abs(psi_2d)**2) * dx2 * dy2

    do step = 1, 10
        call propagate_split_operator_2d(psi_2d, v_pot_2d, dx2, dy2, 1.0_dp, 1.0_dp, 0.05_dp)
    end do
    norm_2d_final = sum(abs(psi_2d)**2) * dx2 * dy2
    call assert_close("2D Split-Operator norm conservation", norm_2d_init, norm_2d_final, 1.0e-5_dp)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Extended Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some extended tests did not pass!"
        stop 1
    else
        print '(A)', "SUCCESS: All extended tests passed."
    end if

contains

    subroutine assert_close(name, expected, actual, eps)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, eps
        real(dp) :: diff

        n_tests = n_tests + 1
        diff = abs(expected - actual)
        if (diff <= eps) then
            n_passed = n_passed + 1
            print '(A, A35, A)', " [PASS] ", name, ""
        else
            print '(A, A35, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_atomic_hhg
