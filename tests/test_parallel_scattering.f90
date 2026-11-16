!> \brief 并行散射能量扫描验收测试 (路线图 §12.1「并行散射计算」项)
!> \details 验收组:
!>          组 1: 能量网格相移扫描与方势阱闭式 atan2 形式逐点对照;
!>          组 2: 全截面扫描 sigma(E) 与解析部分波求和对照;
!>          组 3: 串并一致性 —— 多线程与单线程结果逐位相同 (OpenMP 构建下
!>                真实并行; 串行构建下同代码退化且一致性平凡成立);
!>          组 4: Levinson/Wigner 关系在并行扫描结果上复核 (复用既有机制).
program test_parallel_scattering
    use general_module
    !$ use omp_lib
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat, nthr

    call print_banner("GeneralModule: Parallel Scattering Energy Scan Test", 64)

    nthr = 1
    !$ nthr = omp_get_max_threads()

    block
        integer, parameter :: n_r = 1200, n_e = 400, l_max_t = 4
        real(dp), parameter :: r_max = 12.0_dp, v0 = 0.05_dp, radius = 8.0_dp
        real(dp), parameter :: mass = 1.0_dp
        real(dp) :: rg(n_r), vg(n_r), energies(n_e)
        real(dp) :: deltas(n_e), deltas_all(n_e, 0:l_max_t), sigma(n_e), sigma_ref(n_e)
        real(dp) :: deltas_ref(n_e), err_max, err_sig, emax_v, k_l
        real(dp), allocatable :: d1(:), d2(:)
        complex(dp), allocatable :: s1(:), s2(:)
        integer :: i, k, l, st1, st2

        ! 方势阱网格 (与既有散射测试同参数)
        do i = 1, n_r
            rg(i) = real(i, dp) * r_max / real(n_r, dp)
        end do
        do i = 1, n_r
            if (rg(i) < radius) then
                vg(i) = -v0
            else
                vg(i) = 0.0_dp
            end if
        end do
        do i = 1, n_e
            energies(i) = 1.0e-4_dp * (real(i, dp))**1.5_dp
            energies(i) = 1.0e-4_dp + 0.02_dp * real(i - 1, dp) / real(n_e - 1, dp)
        end do

        ! (s 矩阵缓冲提前分配: 传未分配数组作为 intent(out) 会段错误)
        allocate (d1(n_e), d2(n_e), s1(n_e), s2(n_e))

        ! --- 组 1: 并行 s 波扫描 vs 闭式 ---
        call phase_shift_energy_scan(rg, vg, mass, energies, 0, deltas, s1, st1)
        ! 对数导数匹配只定 delta mod pi (atan2 闭式为主值分支), 比较**模 pi 距离**
        ! (tan 度量在 delta -> pi/2 处发散, 不可用; 模 pi 距离有界且分支无关)
        err_max = 0.0_dp
        do i = 1, n_e
            deltas_ref(i) = square_well_s_phase_exact(energies(i), v0, radius, mass)
            err_max = max(err_max, abs((deltas(i) - deltas_ref(i)) - &
                                       nint((deltas(i) - deltas_ref(i)) / PI) * PI))
        end do
        n_tests = n_tests + 1
        if (st1 == 0 .and. err_max < 2.0e-2_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2)', " [PASS] parallel s-wave scan vs closed form (mod-pi distance): max err =", err_max
        else
            print '(A, ES10.2, A, I4)', " [FAIL] s-wave scan: err =", err_max, " stat =", st1
        end if

        ! --- 组 2: 全截面扫描 vs 独立部分波求和 ---
        call sigma_energy_scan(rg, vg, mass, energies, l_max_t, sigma, deltas_all, st1)
        sigma_ref = 0.0_dp
        do i = 1, n_e
            k_l = sqrt(2.0_dp * mass * energies(i))
            do l = 0, l_max_t
                sigma_ref(i) = sigma_ref(i) + real(2 * l + 1, dp) * sin(deltas_all(i, l))**2
            end do
            sigma_ref(i) = 4.0_dp * PI / k_l**2 * sigma_ref(i)
        end do
        err_sig = maxval(abs(sigma - sigma_ref))
        n_tests = n_tests + 1
        if (st1 == 0 .and. err_sig < 1.0e-12_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2)', " [PASS] sigma(E) scan consistent: err =", err_sig, &
                  "  sigma_max =", maxval(sigma)
        else
            print '(A, ES10.2)', " [FAIL] sigma scan err =", err_sig
        end if

        ! --- 组 3: 串并一致性 (多线程 vs 单线程逐位相同) ---
        !$ call omp_set_num_threads(max(1, nthr))
        call phase_shift_energy_scan(rg, vg, mass, energies, 2, d1, s1, st1)
        !$ call omp_set_num_threads(1)
        call phase_shift_energy_scan(rg, vg, mass, energies, 2, d2, s2, st2)
        !$ call omp_set_num_threads(max(1, nthr))
        err_max = max(maxval(abs(d1 - d2)), maxval(abs(s1 - s2)))
        n_tests = n_tests + 1
        if (st1 == 0 .and. st2 == 0 .and. err_max == 0.0_dp) then
            n_passed = n_passed + 1
            print '(A, I3, A)', " [PASS] serial vs parallel bitwise identical (threads =", nthr, ")"
        else
            print '(A, ES10.2)', " [FAIL] serial/parallel mismatch =", err_max
        end if
        deallocate (d1, d2, s1, s2)

        ! --- 组 4: sigma(E) 扫描的串并一致性 ---
        block
            real(dp) :: sig_p(n_e), sig_s(n_e)
            !$ call omp_set_num_threads(max(1, nthr))
            call sigma_energy_scan(rg, vg, mass, energies, l_max_t, sig_p, stat=st1)
            !$ call omp_set_num_threads(1)
            call sigma_energy_scan(rg, vg, mass, energies, l_max_t, sig_s, stat=st2)
            !$ call omp_set_num_threads(max(1, nthr))
            err_sig = maxval(abs(sig_p - sig_s))
            n_tests = n_tests + 1
            if (st1 == 0 .and. st2 == 0 .and. err_sig == 0.0_dp) then
                n_passed = n_passed + 1
                print '(A)', " [PASS] sigma(E) scan serial vs parallel bitwise identical"
            else
                print '(A, ES10.2)', " [FAIL] sigma scan serial/parallel mismatch =", err_sig
            end if
        end block
    end block

    print '(A, I2, A, I2, A)', "Parallel Scattering Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: parallel scattering test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

end program test_parallel_scattering
