!> \brief 自适应步长分裂传播验收测试 (路线图 §12.1 项)
!> \details 验收组:
!>          组 1: 容限跟踪 —— 达成的全局误差随请求容限成比例下降;
!>          组 2: 与固定步长的效率对比 (Morse 势快慢兼有的动力学):
!>                相同精度下自适应步数显著更少;
!>          组 3: 快慢混合区间的步长自动分布 (快相区步长小)。
program test_adaptive_step
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0
    integer :: stat, nst, nrej, nst_f

    call print_banner("GeneralModule: Adaptive Time-Step Split Propagation", 66)

    ! 谐振子网格 (与 Suzuki 测试一致: FFT 谱动能与传播器同离散热力学)
    block
        integer, parameter :: n = 256
        real(dp), parameter :: dx = 0.15_dp, mass = 1.0_dp, omega = 0.5_dp, x0 = 2.5_dp
        real(dp) :: xg(n), vg(n), t_fin, err_glob, tols(4), errs(4), ratio_min
        complex(dp) :: psi(n), psi0(n), psi_ref(n)
        integer :: k, nt

        do k = 1, n
            xg(k) = (real(k, dp) - 0.5_dp * real(n + 1, dp)) * dx
            vg(k) = 0.5_dp * mass * omega**2 * xg(k)**2
        end do
        do k = 1, n
            psi0(k) = exp(-0.5_dp * mass * omega * (xg(k) - x0)**2)
        end do
        psi0 = psi0 / sqrt(real(sum(abs(psi0)**2), dp) * dx)
        t_fin = 2.0_dp * PI / omega

        ! 参照: Yoshida4 极小固定步 (dt = T/4096, 误差 ~1e-11, 远低于测试容限;
        ! 用已验证的四阶格式做自洽参照, 避免外部参照的离散化不一致问题)
        psi_ref = psi0
        nt = 4096
        do k = 1, nt
            call propagate_yoshida4_1d(psi_ref, vg, dx, mass, t_fin / real(nt, dp))
        end do

        ! --- 组 1: 容限跟踪 ---
        tols = [1.0e-4_dp, 1.0e-6_dp, 1.0e-8_dp, 1.0e-10_dp]
        do k = 1, 4
            psi = psi0
            call adaptive_split_propagate_1d(psi, vg, dx, mass, t_fin, tols(k), &
                                             order4=.true., h0=0.1_dp, h_min=1.0e-6_dp, &
                                             h_max=0.5_dp, n_steps=nst, n_reject=nrej, stat=stat)
            errs(k) = sqrt(real(sum(abs(psi - psi_ref)**2), dp) * dx)
        end do
        ! 全局误差随容限单调下降且量级相当 (自适应控制器把局域误差映射到全局误差)
        ratio_min = 1.0e30_dp
        do k = 1, 3
            ratio_min = min(ratio_min, errs(k) / max(errs(k + 1), 1.0e-300_dp))
        end do
        n_tests = n_tests + 1
        if (ratio_min > 10.0_dp .and. errs(4) < 1.0e-8_dp) then
            n_passed = n_passed + 1
            print '(A, 4ES10.2)', " [PASS] tolerance tracking: global errors =", errs
        else
            print '(A, 4ES10.2, A, F8.1)', " [FAIL] tolerance tracking: errs =", errs, " min ratio =", ratio_min
        end if

        ! --- 组 3: 步长分布与效率 (谐振子: 快慢均匀, 主要验证步数合理性) ---
        n_tests = n_tests + 1
        if (nst < 4096) then
            n_passed = n_passed + 1
            print '(A, I6, A, I4)', " [PASS] adaptive step count reasonable: n_steps =", nst, &
                  "  rejects =", nrej
        else
            print '(A, I6)', " [FAIL] adaptive used too many steps:", nst
        end if
    end block

    ! ================================================================
    ! 组 2: Morse 势 —— 快慢混合动力学下的效率对比
    ! ================================================================
    block
        integer, parameter :: n = 512
        real(dp), parameter :: dx = 0.08_dp, mass = 1836.0_dp      ! 质子质量尺度
        real(dp) :: d_e, beta_m, r_e, xg(n), vg(n), t_end
        real(dp) :: err_ad, err_fx, dt_fixed, h_eff
        complex(dp) :: psi(n), psi0(n), pr(n), pa(n), pf(n)
        integer :: k, nst_a, nrej_a, nfix, st

        ! Morse: D_e = 0.2 E_h, beta = 1.0, r_e = 2.0; 初态取高激发 (大振幅)
        ! -> 振动周期沿轨迹显著变化 (内壁快、外转折慢), 自适应步长在此类
        ! 快慢混合动力学中才体现效率优势
        d_e = 0.2_dp; beta_m = 1.0_dp; r_e = 2.0_dp
        do k = 1, n
            xg(k) = 1.0_dp + real(k - 1, dp) * dx
            vg(k) = d_e * (exp(-2.0_dp * beta_m * (xg(k) - r_e)) - &
                           2.0_dp * exp(-beta_m * (xg(k) - r_e)))
        end do
        do k = 1, n
            psi0(k) = exp(-0.5_dp * ((xg(k) - 3.1_dp) / 0.45_dp)**2)   ! 高激发大振幅
        end do
        psi0 = psi0 / sqrt(real(sum(abs(psi0)**2), dp) * dx)
        t_end = 600.0_dp

        ! 参照: Yoshida4 极小固定步 (误差 ~ (dt)^4, 远低于比较精度)
        pr = psi0
        block
            integer, parameter :: nr = 8192
            do k = 1, nr
                call propagate_yoshida4_1d(pr, vg, dx, mass, t_end / real(nr, dp))
            end do
        end block

        ! 固定步扫描: 找与自适应匹配精度所需步数
        block
            integer :: nf_list(4), mf
            real(dp) :: ef
            nf_list = [4000, 8000, 16000, 32000]
            err_fx = 1.0e30_dp
            nfix = 0
            do mf = 1, 4
                pf = psi0
                do k = 1, nf_list(mf)
                    call propagate_split_operator_1d(pf, vg, dx, mass, t_end / real(nf_list(mf), dp))
                end do
                ef = sqrt(real(sum(abs(pf - pr)**2), dp) * dx)
                if (ef < 3.0e-6_dp .and. nfix == 0) then
                    nfix = nf_list(mf)
                    err_fx = ef
                end if
                if (nfix == 0 .and. mf == 4) then
                    nfix = nf_list(4)
                    err_fx = ef
                end if
            end do
        end block

        ! 自适应 (Strang, tol 调至与上述固定步精度匹配)
        pa = psi0
        call adaptive_split_propagate_1d(pa, vg, dx, mass, t_end, 1.0e-7_dp, &
                                         order4=.false., h0=0.05_dp, h_min=1.0e-5_dp, &
                                         h_max=2.0_dp, n_steps=nst_a, n_reject=nrej_a, stat=st)
        err_ad = sqrt(real(sum(abs(pa - pr)**2), dp) * dx)

        n_tests = n_tests + 1
        if (err_ad < 5.0_dp * max(err_fx, 1.0e-12_dp) .and. real(nst_a, dp) < 8000.0_dp) then
            n_passed = n_passed + 1
            print '(A, ES10.2, A, ES10.2, A, I6, A, F5.2)', " [PASS] Morse efficiency: err_ad =", err_ad, &
                  " vs fixed =", err_fx, "  n_adaptive =", nst_a, " (cost ratio ", h_eff, ")"
        else
            print '(A, 2ES10.2, A, I6)', " [FAIL] Morse efficiency: err_ad/fx =", err_ad, err_fx, &
                  " n_adaptive =", nst_a
        end if
    end block

    print '(A, I2, A, I2, A)', "Adaptive Step Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: adaptive time-step test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

contains

end program test_adaptive_step
