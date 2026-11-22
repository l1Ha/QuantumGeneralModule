!> \brief 量子计算机上的潘宁电离动力学 (参考 Zhuang et al., JPCL 2026, 17, 9513)
!> \details 将 Zikun Zhuang, Chengdong Yang, Yuchen Wang, Dong H. Zhang,
!>          Bin Zhao 的光解离量子算法四步工作流移植到潘宁电离:
!>          1) 二进制网格态初始化 (N=2^nq 点, 坐标寄存器 nq qubits);
!>          2) 二阶分裂算符传播 (DVR 势相位 + QFT 动能相位);
!>          3) 非幺正复光学势 e^{-Gamma(R) dt/2} 经**膨胀(dilation)** ancilla:
!>             U_dil |0>_a |psi> = |0>_a (M|psi>) + |1>_a (sqrt(I-M^2)|psi>),
!>             测 |0>_a 投影到 M|psi>, 概率 p0 = ||M psi||^2
!>             (论文 eq. (6)-(8); 本库的 Gamma(R) 就是物理吸收,
!>              与光解离把 CAP 当人为边界不同 —— 后选择失败率即电离率);
!>          4) 读出: 存活概率 p_surv = prod p0, P_ion = 1 - p_surv;
!>             可选 Hadamard test 自相关 C(t) = <psi0|psi(t)> (论文 eq. 10).
!>          态矢模拟器精确复现该电路 (无采样噪声), 用于可行性验收:
!>          与库内经典复势传播/轨迹电离概率对表 (`mod_penning_...`)。
!>          资源: nq + 1 (膨胀) (+1 若取 C(t)) qubits, 深度 ~ 2 * n_steps * (QFT + 相位 + 旋转).
!> \author LiHao
module mod_quantum_hardware
    use mod_constants, only: dp, PI, AU2EV
    use mod_penning_associative_ionization, only: penning_system_t, &
        eval_penning_potentials
    implicit none
    private

    public :: qhw_penning_result_t
    public :: qhw_penning_on_circuit

    !> \brief 量子电路模拟结果 (潘宁电离)
    type, public :: qhw_penning_result_t
        real(dp) :: p_ion          !< 电离概率 = 1 - 存活概率
        real(dp) :: p_surv         !< 存活概率 (全步后选择连乘)
        real(dp) :: p_step_min     !< 单步最小后选择概率 (诊断)
        integer  :: n_grid         !< 网格点数 N = 2^nq
        integer  :: n_qubits       !< 坐标寄存器 qubit 数 nq
        integer  :: n_ancilla      !< ancilla 数 (膨胀 1; 含 Hadamard 测试 2)
        integer  :: n_steps        !< 传播步数
        real(dp) :: dt_au          !< 时间步 (a.u.)
        real(dp) :: e_coll_ev      !< 碰撞能 (eV)
    end type qhw_penning_result_t

contains

    !> \brief 在态矢模拟器上执行论文四步工作流 (潘宁电离)
    !> \param[in] sys 潘宁体系 (势能与宽度参数)
    !> \param[in] e_coll_ev 碰撞能 (eV)
    !> \param[in] nq 坐标寄存器 qubit 数 (N = 2^nq 点; 建议 5-8)
    !> \param[in] n_steps 总时间步数 (dt = t_total / n_steps)
    !> \param[in] t_total_eff 有效传播时长 (a.u.; 应覆盖穿越吸收区)
    !> \param[in] r_min, r_max 网格区间 (Angstrom)
    !> \param[in] r0_ang, sig_ang 入射波包中心与宽度 (Angstrom)
    !> \param[in] use_hadamard 是否同时模拟 Hadamard test 自相关 (资源 +1 ancilla)
    subroutine qhw_penning_on_circuit(sys, e_coll_ev, nq, n_steps, t_total_eff, &
                                      r_min, r_max, r0_ang, sig_ang, use_hadamard, &
                                      res, stat)
        type(penning_system_t), intent(in) :: sys
        real(dp), intent(in) :: e_coll_ev, t_total_eff, r_min, r_max, r0_ang, sig_ang
        integer, intent(in) :: nq, n_steps
        logical, intent(in), optional :: use_hadamard
        type(qhw_penning_result_t), intent(out) :: res
        integer, intent(out) :: stat

        complex(dp), allocatable :: psi(:), psi0(:), spec(:)
        real(dp), allocatable :: r_g(:), v_star_ev(:), v_plus_ev(:), gam_ev(:)
        real(dp), allocatable :: v_au(:), gam_au(:), kin(:), amp(:), abs_g(:)
        real(dp) :: dr, k0, mu_au, v_rel, p_step, nrm2, dt
        integer :: n, i, sj, s, nq_eff
        logical :: had

        stat = 0
        had = .false.
        if (present(use_hadamard)) had = use_hadamard
        nq_eff = nq
        n = 2**nq_eff
        res%n_grid = n
        res%n_qubits = nq_eff
        res%n_ancilla = 1 + merge(1, 0, had)
        res%n_steps = n_steps
        res%e_coll_ev = e_coll_ev
        dt = t_total_eff / real(max(n_steps, 1), dp)
        res%dt_au = dt

        allocate (psi(n), psi0(n), spec(n))
        allocate (r_g(n), v_star_ev(n), v_plus_ev(n), gam_ev(n))
        allocate (v_au(n), gam_au(n), kin(n), amp(n), abs_g(n))

        dr = (r_max - r_min) / real(n - 1, dp)
        do i = 1, n
            r_g(i) = r_min + real(i - 1, dp) * dr
            call eval_penning_potentials(sys, r_g(i), v_star_ev(i), &
                                         v_plus_ev(i), gam_ev(i))
            v_au(i) = v_star_ev(i) / AU2EV
            gam_au(i) = gam_ev(i) / AU2EV
        end do
        mu_au = sys%mu_au
        k0 = sqrt(2.0_dp * mu_au * e_coll_ev / AU2EV)   ! 1/bohr, 向 R 减小方向

        ! ---- 1) 态初始化: 二进制网格编码的高斯波包 ----
        psi0 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n
            psi0(i) = exp(-((r_g(i) - r0_ang) / sig_ang)**2) &
                      * cmplx(cos(k0 * r_g(i) / 0.5291772109_dp), &
                             -sin(k0 * r_g(i) / 0.5291772109_dp), kind=dp)
        end do
        psi0 = psi0 / sqrt(real(sum(abs(psi0)**2), dp))
        psi = psi0

        ! ---- 动能谱 (FBR): k_j = 2 pi (j-1-N/2) / (N dR) ----
        do i = 1, n
            if (i <= n / 2) then
                kin(i) = 2.0_dp * PI * real(i - 1, dp) / (real(n, dp) * dr)
            else
                kin(i) = 2.0_dp * PI * real(i - 1 - n, dp) / (real(n, dp) * dr)
            end if
            kin(i) = kin(i) / 0.5291772109_dp               ! 1/bohr
            kin(i) = kin(i)**2 / (2.0_dp * mu_au)          ! a.u.
        end do
        ! 膨胀矩阵对角: M = exp(-Gamma dt/2) -> 振幅 sqrt(exp(-Gamma dt))
        amp = sqrt(exp(-gam_au * dt))
        abs_g = gam_au

        ! ---- 2) 分裂算符传播 + 3) 膨胀后选择 ----
        res%p_surv = 1.0_dp
        res%p_step_min = 1.0_dp
        do s = 1, max(n_steps, 1)
            ! 二阶分裂: e^{-iV dt/2} e^{-iT dt} e^{-iV dt/2}
            do i = 1, n
                psi(i) = psi(i) * cmplx(cos(0.5_dp * v_au(i) * dt), &
                                       -sin(0.5_dp * v_au(i) * dt), kind=dp)
            end do
            spec = psi
            call fft_complex(spec)
            do i = 1, n
                spec(i) = spec(i) * cmplx(cos(kin(i) * dt), -sin(kin(i) * dt), kind=dp)
            end do
            call fft_complex_inv(spec)
            psi = spec
            do i = 1, n
                psi(i) = psi(i) * cmplx(cos(0.5_dp * v_au(i) * dt), &
                                       -sin(0.5_dp * v_au(i) * dt), kind=dp)
            end do
            ! 膨胀 / 后选择
            do i = 1, n
                psi(i) = psi(i) * amp(i)
            end do
            nrm2 = real(sum(abs(psi)**2), dp)
            if (nrm2 > 1.0e-300_dp) then
                p_step = nrm2
            else
                p_step = 1.0e-300_dp
                psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
            end if
            res%p_surv = res%p_surv * p_step
            res%p_step_min = min(res%p_step_min, p_step)
            if (nrm2 > 1.0e-300_dp) psi = psi / sqrt(nrm2)
        end do
        res%p_ion = 1.0_dp - res%p_surv
        ! Hadamard test 端点 (可选): C(t_end) = <psi0|psi> (态矢模拟给出期望值)
        if (had) then
            block
                complex(dp) :: cend
                cend = sum(conjg(psi0) * psi)
                write (*, '(A, 2F10.6)') "  [qhw] C(Te) = ", real(cend, dp), aimag(cend)
            end block
        end if

        deallocate (psi, psi0, spec, r_g, v_star_ev, v_plus_ev, gam_ev)
        deallocate (v_au, gam_au, kin, amp, abs_g)
        stat = 0
    end subroutine qhw_penning_on_circuit

    !> \brief 复数基 2 FFT (迭代 Cooley-Tukey: 比特反序 + 蝶形, 长度 2^n)
    subroutine fft_complex(v)
        complex(dp), intent(inout) :: v(:)
        complex(dp) :: t, w
        real(dp) :: ang
        integer :: n, i, j, len, half, k, m
        n = size(v)
        if (n <= 1) return
        ! 比特反序
        j = 0
        do i = 2, n
            m = n / 2
            do while (j >= m .and. m >= 1)
                j = j - m
                m = m / 2
            end do
            j = j + m
            if (i < j + 1) then
                t = v(i)
                v(i) = v(j + 1)
                v(j + 1) = t
            end if
        end do
        ! 蝶形: len = 2, 4, ..., n
        len = 2
        do while (len <= n)
            half = len / 2
            ang = -2.0_dp * PI / real(len, dp)
            do k = 0, half - 1
                w = cmplx(cos(ang * real(k, dp)), sin(ang * real(k, dp)), kind=dp)
                do i = k + 1, n, len
                    t = w * v(i + half)
                    v(i + half) = v(i) - t
                    v(i) = v(i) + t
                end do
            end do
            len = len * 2
        end do
    end subroutine fft_complex

    !> \brief 逆 FFT (共轭法 + 1/N)
    subroutine fft_complex_inv(v)
        complex(dp), intent(inout) :: v(:)
        integer :: n
        n = size(v)
        v = conjg(v)
        call fft_complex(v)
        v = conjg(v) / real(n, dp)
    end subroutine fft_complex_inv

end module mod_quantum_hardware
