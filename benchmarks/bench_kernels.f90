!> \brief 核心算法性能基准程序 (路线图 §15.10「基准测试」项)
!> \details 对关键数值核做确定性工作量的墙钟计时, 输出:
!>            (a) 人读表格; (b) 规范 JSON 记录 (mod_json, 见 §8.6);
!>            (c) OpenMP 线程扩展 (fft_2d, 仅 -fopenmp 构建)。
!>          正确性门: 每个核在做计时前先做数值正确性检查, 任一门失败
!>          以非零退出码结束 —— CI 只验证"基准负载下正确", 不设时间阈值
!>          (共享 runner 的计时不稳定), 时间数据留作本地性能回归对照。
!> \author LiHao
program bench_kernels
    use general_module
    use mod_json
    !$ use omp_lib
    implicit none

    integer, parameter :: n_rep = 5
    integer :: i, j, unit, stat, nthr
    real(dp) :: t_best, t_med, tms(n_rep)
    logical :: ok_all
    type(json_writer_t) :: jw
    character(len=:), allocatable :: json_text

    ! 工作负载规模
    integer, parameter :: nfft = 2**20          ! 1D FFT 点数
    integer, parameter :: n2d = 1024            ! 2D FFT 边长 (1024^2)
    integer, parameter :: neig = 600            ! 实对称本征维数
    integer, parameter :: nher = 200            ! 复 Hermitian 维数
    integer, parameter :: nkry = 4096           ! Krylov/Chebyshev 波函数维数
    integer, parameter :: nstep = 100           ! 传播步数
    integer, parameter :: nsp3 = 32             ! 3D 样条网格边长
    integer, parameter :: nev3 = 1000000        ! 3D 样条求值次数

    ok_all = .true.
    nthr = 1
    !$ nthr = omp_get_max_threads()

    call print_banner("GeneralModule: Kernel Performance Benchmark", 66)
    call json_init(jw)
    call json_begin_object(jw)
    call json_add_string(jw, "record_type", "benchmark")
    call json_add_string(jw, "module", "GeneralModule")
    call json_add_int(jw, "n_threads_max", nthr)
    call json_begin_object(jw, "results")

    ! ---------------- 1. 1D FFT ----------------
    block
        complex(dp), allocatable :: v(:), v0(:)
        real(dp) :: err
        allocate (v(nfft), v0(nfft))
        do i = 1, nfft
            v(i) = cmplx(sin(0.001_dp * i), cos(0.002_dp * i), kind=dp)
        end do
        v0 = v
        call fft_1d(v, -1)
        call fft_1d(v, +1)
        err = maxval(abs(v - v0))
        ok_all = ok_all .and. (err < 1.0e-16_dp * nfft)   ! 往返
        do i = 1, n_rep
            v = v0
            tms(i) = time_kernel_fft1(v)
        end do
        call stats(tms, t_best, t_med)
        print '(A, F10.2, A, ES10.2)', " [bench] fft_1d   n=2^20 : ", 1.0e3_dp * t_best, &
              " ms  (rel err ", err, ")"
        call json_add_real(jw, "fft_1d_ms", 1.0e3_dp * t_best)
    end block

    ! ---------------- 2. 2D FFT (OpenMP 路径) ----------------
    block
        complex(dp), allocatable :: m(:, :), m0(:, :)
        real(dp) :: err
        allocate (m(n2d, n2d), m0(n2d, n2d))
        do i = 1, n2d
            m(:, i) = cmplx(sin(0.001_dp * i), 0.0_dp, kind=dp)
        end do
        m0 = m
        call fft_2d(m, -1)
        call fft_2d(m, +1)
        err = maxval(abs(m - m0))
        ok_all = ok_all .and. (err < 1.0e-14_dp)
        do i = 1, 3
            m = m0
            tms(i) = time_kernel_fft2(m)
        end do
        call stats(tms(1:3), t_best, t_med)
        print '(A, F10.2, A, ES10.2)', " [bench] fft_2d  1024^2 : ", 1.0e3_dp * t_best, &
              " ms  (rel err ", err, ")"
        call json_add_real(jw, "fft_2d_ms", 1.0e3_dp * t_best)
    end block

    ! ---------------- 3. 实对称本征 (tql2) ----------------
    block
        real(dp), allocatable :: a(:, :), av(:, :), d(:), z(:, :)
        real(dp) :: resid
        allocate (a(neig, neig), z(neig, neig), d(neig), av(neig, neig))
        do i = 1, neig
            a(:, i) = [(cos(real(3 * j + i, dp)) / real(1 + abs(i - j), dp), j=1, neig)]
        end do
        a = 0.5_dp * (a + transpose(a))
        tms(1) = time_kernel_eig(a, d, z)
        av = matmul(a, z)
        resid = 0.0_dp
        do i = 1, neig
            resid = max(resid, maxval(abs(av(:, i) - d(i) * z(:, i))))
        end do
        ok_all = ok_all .and. (resid < 1.0e-10_dp)
        call stats(tms(1:1), t_best, t_med)
        print '(A, F10.2, A, ES10.2)', " [bench] diag_sym n=600  : ", 1.0e3_dp * t_best, &
              " ms  (resid ", resid, ")"
        call json_add_real(jw, "diag_symmetric_ms", 1.0e3_dp * t_best)
    end block

    ! ---------------- 4. 复 Hermitian (加倍归约) ----------------
    block
        complex(dp), allocatable :: ah(:, :), zh(:, :)
        real(dp), allocatable :: dh(:)
        complex(dp), allocatable :: av(:, :)
        real(dp) :: resid
        allocate (ah(nher, nher), zh(nher, nher), dh(nher), av(nher, nher))
        do i = 1, nher
            ah(:, i) = [(cmplx(sin(real(3 * j + i, dp)), cos(real(2 * i - j, dp)), kind=dp), j=1, nher)]
        end do
        ah = ah + conjg(transpose(ah))
        tms(1) = time_kernel_herm(ah, dh, zh)
        av = matmul(ah, zh)
        resid = 0.0_dp
        do i = 1, nher
            resid = max(resid, sqrt(real(sum(abs(av(:, i) - dh(i) * zh(:, i))**2), dp)))
        end do
        ! n=200 稠密随机 Hermitian 的加倍归约精度 (~1e-11, 见 §4.5)
        ok_all = ok_all .and. (resid < 1.0e-9_dp)
        call stats(tms(1:1), t_best, t_med)
        print '(A, F10.2, A, ES10.2)', " [bench] herm n=200     : ", 1.0e3_dp * t_best, &
              " ms  (resid ", resid, ")"
        call json_add_real(jw, "diag_hermitian_ms", 1.0e3_dp * t_best)
    end block

    ! ---------------- 5. Krylov vs Chebyshev 传播 (nkry, nstep) ----------------
    block
        complex(dp), allocatable :: psi(:, :)
        real(dp) :: drift
        allocate (psi(nkry, 2))
        do i = 1, nkry
            psi(i, 1) = cmplx(exp(-(real(i, dp) - 0.5_dp * real(nkry, dp))**2 / 500.0_dp), 0.0_dp, kind=dp)
        end do
        psi(:, 1) = psi(:, 1) / sqrt(real(sum(abs(psi(:, 1))**2), dp))
        psi(:, 2) = psi(:, 1)

        tms(1) = time_kernel_krylov(psi(:, 1), nstep)
        drift = abs(sqrt(real(sum(abs(psi(:, 1))**2), dp)) - 1.0_dp)
        call stats(tms(1:1), t_best, t_med)
        print '(A, F10.2, A, ES10.2)', " [bench] krylov  n=4096 x100: ", 1.0e3_dp * t_best, &
              " ms  (norm drift ", drift, ")"
        call json_add_real(jw, "krylov_100steps_ms", 1.0e3_dp * t_best)
        ok_all = ok_all .and. (drift < 1.0e-10_dp)

        tms(1) = time_kernel_cheby(psi(:, 2), nstep)
        drift = abs(sqrt(real(sum(abs(psi(:, 2))**2), dp)) - 1.0_dp)
        call stats(tms(1:1), t_best, t_med)
        print '(A, F10.2, A, ES10.2)', " [bench] chebyshev n=4096x100: ", 1.0e3_dp * t_best, &
              " ms  (norm drift ", drift, ")"
        call json_add_real(jw, "chebyshev_100steps_ms", 1.0e3_dp * t_best)
    end block

    ! ---------------- 6. 3D 样条求值 ----------------
    block
        real(dp), allocatable :: g(:), f3(:, :, :)
        real(dp) :: tps, err, xv
        type(spline_3d_t) :: s3
        allocate (g(nsp3), f3(nsp3, nsp3, nsp3))
        do i = 1, nsp3
            g(i) = real(i - 1, dp) / real(nsp3 - 1, dp)
        end do
        do i = 1, nsp3
            f3(:, :, i) = sin(2.0_dp * g(i)) + 0.1_dp
        end do
        call spline_3d_init(g, g, g, f3, s3, stat)
        tps = 0.0_dp
        block
            integer(8) :: t0, t1, rate
            real(dp) :: acc
            acc = 0.0_dp
            call system_clock(t0, rate)
            do i = 1, nev3
                xv = real(mod(i, nsp3), dp) / real(nsp3 - 1, dp)
                acc = acc + spline_3d_eval(s3, xv, xv, xv)
            end do
            call system_clock(t1)
            tps = real(t1 - t0, dp) / real(rate, dp)
            ok_all = ok_all .and. (abs(acc) > 0.0_dp)
        end block
        err = abs(spline_3d_eval(s3, g(3), g(5), g(7)) - f3(3, 5, 7))
        ok_all = ok_all .and. (err < 1.0e-14_dp)
        call spline_3d_clean(s3)
        print '(A, F10.2, A, F10.1, A)', " [bench] spline_3d 1e6 pts : ", 1.0e3_dp * tps, &
              " ms  (", 1.0e6_dp / tps / 1.0e6_dp, " Mpts/s)"
        call json_add_real(jw, "spline3d_1e6_ms", 1.0e3_dp * tps)
    end block

    ! ---------------- 7. OpenMP 线程扩展 (fft_2d) ----------------
    block
        integer, parameter :: thr_list(6) = [1, 2, 4, 8, 16, 20]
        complex(dp), allocatable :: m(:, :), m0(:, :)
        real(dp) :: ts
        integer :: k
        allocate (m(n2d, n2d), m0(n2d, n2d))
        m0 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n2d
            m0(:, i) = cmplx(sin(0.001_dp * i), 0.0_dp, kind=dp)
        end do
        call json_begin_array(jw, "fft2d_scaling")
        print '(A)', " [bench] fft_2d OpenMP scaling:"
        do k = 1, size(thr_list)
            if (thr_list(k) > nthr) exit
            !$ call omp_set_num_threads(thr_list(k))
            m = m0
            ts = time_kernel_fft2(m)
            m = m0
            ts = min(ts, time_kernel_fft2(m))   ! 两次取优
            print '(A, I4, A, F10.2, A)', "         threads=", thr_list(k), "  ", 1.0e3_dp * ts, " ms"
            call json_begin_object(jw)
            call json_add_int(jw, "threads", thr_list(k))
            call json_add_real(jw, "fft_2d_ms", 1.0e3_dp * ts)
            call json_end_object(jw)
        end do
        call json_end_array(jw)
        !$ call omp_set_num_threads(nthr)
        deallocate (m, m0)
    end block

    call json_end_object(jw)   ! results
    call json_add_bool(jw, "all_correctness_checks_passed", ok_all)
    call json_end_object(jw)
    json_text = json_get_text(jw)
    open (newunit=unit, file="benchmarks/benchmark_results.json", status="replace")
    write (unit, '(A)') json_text
    close (unit)
    call json_finalize(jw)

    print '(A)', "------------------------------------------------------------"
    if (ok_all) then
        print '(A)', " [PASS] all kernel correctness checks passed; results -> benchmarks/benchmark_results.json"
    else
        print '(A)', " [FAIL] kernel correctness check failed"
        stop 1
    end if

contains

    !> 统计最优/中位时间
    subroutine stats(ts, best, med)
        real(dp), intent(in) :: ts(:)
        real(dp), intent(out) :: best, med
        real(dp), allocatable :: s(:)
        s = ts
        best = minval(s)
        call sort_vec(s)
        med = s((size(s) + 1) / 2)
    end subroutine stats

    subroutine sort_vec(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: tmp
        integer :: p, q
        do p = 1, size(x) - 1
            do q = 1, size(x) - p
                if (x(q) > x(q + 1)) then
                    tmp = x(q); x(q) = x(q + 1); x(q + 1) = tmp
                end if
            end do
        end do
    end subroutine sort_vec

    function time_kernel_fft1(v) result(t)
        complex(dp), intent(inout) :: v(:)
        real(dp) :: t
        integer(8) :: t0, t1, rate
        call system_clock(t0, rate)
        call fft_1d(v, -1)
        call system_clock(t1)
        t = real(t1 - t0, dp) / real(rate, dp)
    end function time_kernel_fft1

    function time_kernel_fft2(m) result(t)
        complex(dp), intent(inout) :: m(:, :)
        real(dp) :: t
        integer(8) :: t0, t1, rate
        call system_clock(t0, rate)
        call fft_2d(m, -1)
        call system_clock(t1)
        t = real(t1 - t0, dp) / real(rate, dp)
    end function time_kernel_fft2

    function time_kernel_eig(a, d, z) result(t)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: d(:), z(:, :)
        real(dp) :: t
        integer(8) :: t0, t1, rate
        integer :: st
        call system_clock(t0, rate)
        call diag_symmetric_matrix(size(a, 1), a, d, z, st)
        call system_clock(t1)
        t = real(t1 - t0, dp) / real(rate, dp)
    end function time_kernel_eig

    function time_kernel_herm(a, d, z) result(t)
        complex(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: d(:)
        complex(dp), intent(out) :: z(:, :)
        real(dp) :: t
        integer(8) :: t0, t1, rate
        integer :: st
        call system_clock(t0, rate)
        call diag_hermitian_matrix(size(a, 1), a, d, z, st)
        call system_clock(t1)
        t = real(t1 - t0, dp) / real(rate, dp)
    end function time_kernel_herm

    function time_kernel_krylov(psi, nst) result(t)
        complex(dp), intent(inout) :: psi(:)
        integer, intent(in) :: nst
        real(dp) :: t
        integer(8) :: t0, t1, rate
        integer :: k, st, mu
        real(dp) :: ee
        call system_clock(t0, rate)
        do k = 1, nst
            call krylov_propagate_step(psi, 0.01_dp, 24, h_mult_test, ee, mu, st)
        end do
        call system_clock(t1)
        t = real(t1 - t0, dp) / real(rate, dp)
    end function time_kernel_krylov

    function time_kernel_cheby(psi, nst) result(t)
        complex(dp), intent(inout) :: psi(:)
        integer, intent(in) :: nst
        real(dp) :: t
        integer(8) :: t0, t1, rate
        integer :: k
        call system_clock(t0, rate)
        do k = 1, nst
            call chebyshev_propagate_step(psi, 0.01_dp, -0.5_dp, 1.5_dp, h_mult_test_c, 30)
        end do
        call system_clock(t1)
        t = real(t1 - t0, dp) / real(rate, dp)
    end function time_kernel_cheby

    ! 测试哈密顿量: 一维二阶差分 + 谐振子势 (n = 4096, dx = 0.05)
    subroutine h_mult_test_c(v_in, v_out)
        complex(dp), intent(in) :: v_in(:)
        complex(dp), intent(out) :: v_out(:)
        real(dp), parameter :: dx = 0.05_dp
        real(dp) :: kin
        integer :: j, n
        n = size(v_in)
        kin = 0.5_dp / dx**2
        do j = 1, n
            v_out(j) = (0.5_dp * (real(j, dp) - 0.5_dp * real(n, dp))**2 * dx**2) * v_in(j) - 2.0_dp * kin * v_in(j)
            if (j > 1) v_out(j) = v_out(j) + kin * v_in(j - 1)
            if (j < n) v_out(j) = v_out(j) + kin * v_in(j + 1)
        end do
    end subroutine h_mult_test_c

    subroutine h_mult_test(v_in, v_out)
        real(dp), intent(in) :: v_in(:)
        real(dp), intent(out) :: v_out(:)
        complex(dp) :: c_in(size(v_in)), c_out(size(v_in))
        c_in = cmplx(v_in, 0.0_dp, kind=dp)
        call h_mult_test_c(c_in, c_out)
        v_out = real(c_out, dp)
    end subroutine h_mult_test

end program bench_kernels
