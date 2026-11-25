!> \brief GPU/OpenACC 热核基准: sop_apply_3d (3D 和积作用) 与 fft_1d
!> \details 计时对比 CPU (gfortran [+OpenMP]) 与 GPU (nvfortran -acc=gpu)。
!>          用法: gfortran/ nvfortran 编译后运行; 输出各核平均耗时。
program bench_gpu
    use mod_constants, only: dp, PI
    use mod_linear_algebra, only: fft_1d
    use mod_sop_hamiltonian
    implicit none

    integer, parameter :: n3d = 48, n_terms = 8, n_fft = 2**22, reps = 5
    type(sop_hamiltonian_3d_t) :: sop
    real(dp), allocatable :: vin(:, :, :), vout(:, :, :)
    complex(dp), allocatable :: cv(:)
    real(dp) :: t0, t1, t_sop, t_fft
    integer :: m, i1, i2, k, r

    allocate (vin(n3d, n3d, n3d), vout(n3d, n3d, n3d))
    call random_number(vin)
    call init_sop_hamiltonian_3d(sop, n3d, n3d, n3d, n_terms)
    do m = 1, n_terms
        sop%coeffs(m) = 1.0_dp / real(m, dp)
        do k = 1, n3d
            sop%h1(m)%mat(k, k) = 2.0_dp + real(m, dp) * 0.1_dp
            sop%h2(m)%mat(k, k) = 2.0_dp + real(m, dp) * 0.1_dp
            sop%h3(m)%mat(k, k) = 2.0_dp + real(m, dp) * 0.1_dp
            if (k < n3d) then
                sop%h1(m)%mat(k, k + 1) = -0.5_dp
                sop%h1(m)%mat(k + 1, k) = -0.5_dp
                sop%h2(m)%mat(k, k + 1) = -0.5_dp
                sop%h2(m)%mat(k + 1, k) = -0.5_dp
                sop%h3(m)%mat(k, k + 1) = -0.5_dp
                sop%h3(m)%mat(k + 1, k) = -0.5_dp
            end if
        end do
    end do
    ! 预热 + 计时
    call sop_apply_3d(sop, vin, vout)
    t1 = walltime()
    do r = 1, reps
        call sop_apply_3d(sop, vin, vout)
    end do
    t_sop = (walltime() - t1) / real(reps, dp)

    allocate (cv(n_fft))
    cv = cmplx(0.1_dp, 0.05_dp, kind=dp)
    call fft_1d(cv, -1)
    t1 = walltime()
    do r = 1, reps
        call fft_1d(cv, -1)
    end do
    t_fft = (walltime() - t1) / real(reps, dp)

    write (*, '(A, I4, A, I2, A, F10.4, A)') "SOP 3D n=", n3d, " M=", n_terms, &
          "  t=", t_sop * 1e3_dp, " ms"
    write (*, '(A, I9, A, F10.4, A)') "FFT n=", n_fft, "  t=", t_fft * 1e3_dp, " ms"
    deallocate (vin, vout, cv)

contains
    function walltime() result(t)
        real(dp) :: t
        integer :: cnt, rate
        call system_clock(cnt, rate)
        t = real(cnt, dp) / real(rate, dp)
    end function walltime
end program bench_gpu
