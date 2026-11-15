!> \brief Lanczos 短迭代 Krylov 传播器模块 (SIL, §14.1 路线图项实现)
!> \details 对实对称、时间无关哈密顿量 H, 以 m 步 Lanczos 全重正交三对角化
!>          构造 Krylov 子空间 K_m(H, |psi>) = span{|psi>, H|psi>, ..., H^{m-1}|psi>},
!>          在子空间内对三对角矩阵 T_m 精确求矩阵指数, 实现
!>            psi(t + dt) = exp(-i H dt) psi(t) ≈ V_m exp(-i T_m dt) e_1 ||psi||.
!>          关键性质:
!>            (1) 范数守恒到机器精度 (V_m 正交 + T_m 指数幺正);
!>            (2) 误差随 m 超指数衰减;
!>            (3) 不变子空间精确终止 —— 若初态张成 d 维不变子空间,
!>                Lanczos 在第 d 步 beta_{d+1} = 0, 传播变为严格精确
!>                (二能级 d = 2, 任意本征态 d = 1);
!>            (4) 后验误差估计 err ≈ beta_{m+1} |[exp(-i T_m dt)]_{m,1}|。
!>          复波函数按实部/虚部分别传播 (H 实对称 ⇒ 两分量独立演化)。
!>          边界: 本模块要求 H 实对称 (无磁场、无复吸收势); 复势
!>          (CAP) 情形请用切比雪夫/劈裂算符路径 (§6.7)。
!> \author LiHao
module mod_krylov_propagator
    use mod_constants, only: dp
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: lanczos_tridiag
    public :: krylov_propagate_step

contains

    !> \brief m 步 Lanczos 三对角化 (两遍全重正交, 数值稳定)
    !> \param[in]  n_dim    向量维度
    !> \param[in]  m_max    最大 Krylov 步数
    !> \param[in]  v_start  起始实向量 (非零)
    !> \param[in]  r_mult   实对称哈密顿量矩阵-向量乘法子程序
    !> \param[out] alpha    三对角对角元 alpha(1:m_eff)
    !> \param[out] beta     次对角元 beta(2:m_eff), beta(1) = 0;
    !>                      beta(m_eff+1) 为终止时剩余范数 (不变子空间时为 0)
    !> \param[out] v_basis  正交归一闪 Krylov 基 (n_dim, m_max), 前 m_eff 列有效
    !> \param[out] m_eff    实际有效步数 (含不变子空间精确终止)
    !> \param[out] stat     状态码 (0 正常; 1 起始向量为零)
    subroutine lanczos_tridiag(n_dim, m_max, v_start, r_mult, alpha, beta, &
                               v_basis, m_eff, stat)
        integer, intent(in) :: n_dim, m_max
        real(dp), intent(in) :: v_start(n_dim)
        interface
            subroutine r_mult(v_in, v_out)
                use mod_constants, only: dp
                real(dp), intent(in) :: v_in(:)
                real(dp), intent(out) :: v_out(:)
            end subroutine r_mult
        end interface
        real(dp), intent(out) :: alpha(m_max), beta(m_max + 1)
        real(dp), intent(out) :: v_basis(n_dim, m_max)
        integer, intent(out) :: m_eff, stat

        real(dp), allocatable :: w(:), v_prev(:)
        real(dp) :: nrm, proj, tol_break
        integer :: i, j

        stat = 0
        alpha = 0.0_dp
        beta = 0.0_dp
        v_basis = 0.0_dp

        allocate(w(n_dim), v_prev(n_dim))
        v_prev = 0.0_dp

        nrm = sqrt(dot_product(v_start, v_start))
        if (nrm <= tiny(1.0_dp)) then
            stat = 1
            m_eff = 0
            deallocate(w, v_prev)
            return
        end if
        v_basis(:, 1) = v_start / nrm
        tol_break = 1.0e-14_dp * nrm

        m_eff = m_max
        do j = 1, m_max
            call r_mult(v_basis(:, j), w)
            alpha(j) = dot_product(v_basis(:, j), w)
            w = w - alpha(j) * v_basis(:, j)
            if (j > 1) w = w - beta(j) * v_prev
            ! 两遍经典 Gram-Schmidt 全重正交 (V_m^T V_m = I 到机器精度)
            do i = 1, j
                proj = dot_product(v_basis(:, i), w)
                w = w - proj * v_basis(:, i)
            end do
            do i = 1, j
                proj = dot_product(v_basis(:, i), w)
                w = w - proj * v_basis(:, i)
            end do
            nrm = sqrt(dot_product(w, w))
            beta(j + 1) = nrm
            if (nrm <= tol_break) then
                ! 不变子空间: Krylov 空间维数耗尽, 该子空间内传播严格精确
                m_eff = j
                exit
            end if
            if (j < m_max) then
                v_basis(:, j + 1) = w / nrm
                v_prev = v_basis(:, j)
            end if
        end do

        deallocate(w, v_prev)
    end subroutine lanczos_tridiag

    !> \brief Lanczos 短迭代传播单步: psi(t+dt) = exp(-i H dt) psi(t)
    !> \param[inout] psi     复数波函数 (n_dim), 原位更新
    !> \param[in]    dt      时间步长 (a.u.)
    !> \param[in]    m_krylov Krylov 子空间最大维数 (通常 10~30; 不变子空间自动提前终止)
    !> \param[in]    r_mult  实对称哈密顿量矩阵-向量乘法子程序
    !> \param[out]   err_est 后验误差估计 beta_{m+1}|[exp(-iT dt)]_{m,1}| (可选)
    !> \param[out]   m_used  实际使用的 Krylov 维数 (可选)
    !> \param[out]   stat    状态码 (0 正常; 2 Lanczos/对角化异常)
    subroutine krylov_propagate_step(psi, dt, m_krylov, r_mult, err_est, m_used, stat)
        complex(dp), intent(inout) :: psi(:)
        real(dp), intent(in) :: dt
        integer, intent(in) :: m_krylov
        interface
            subroutine r_mult(v_in, v_out)
                use mod_constants, only: dp
                real(dp), intent(in) :: v_in(:)
                real(dp), intent(out) :: v_out(:)
            end subroutine r_mult
        end interface
        real(dp), intent(out), optional :: err_est
        integer, intent(out), optional :: m_used
        integer, intent(out) :: stat

        real(dp), allocatable :: v_re(:), v_im(:)
        complex(dp), allocatable :: part_c(:)
        real(dp) :: nrm_re, nrm_im, err_re, err_im
        integer :: n, mm, m_eff_re, m_eff_im, st_lan, st_diag

        stat = 0
        n = size(psi)
        mm = max(2, m_krylov)
        allocate(v_re(n), v_im(n), part_c(n))
        v_re = real(psi, dp)
        v_im = aimag(psi)
        nrm_re = sqrt(dot_product(v_re, v_re))
        nrm_im = sqrt(dot_product(v_im, v_im))
        err_re = 0.0_dp
        err_im = 0.0_dp
        m_eff_re = 0
        m_eff_im = 0

        ! 实部
        if (nrm_re > tiny(1.0_dp)) then
            call lanczos_apply_expm(n, mm, v_re, nrm_re, dt, r_mult, part_c, &
                                    err_re, m_eff_re, st_lan, st_diag)
            if (st_lan /= 0 .or. st_diag /= 0) stat = 2
            psi = part_c
        else
            psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end if

        ! 虚部 (独立演化, 乘 i 后叠加)
        if (nrm_im > tiny(1.0_dp)) then
            call lanczos_apply_expm(n, mm, v_im, nrm_im, dt, r_mult, part_c, &
                                    err_im, m_eff_im, st_lan, st_diag)
            if (st_lan /= 0 .or. st_diag /= 0) stat = 2
            psi = psi + cmplx(0.0_dp, 1.0_dp, kind=dp) * part_c
        end if

        if (present(err_est)) err_est = max(err_re, err_im)
        if (present(m_used)) m_used = max(m_eff_re, m_eff_im)
        deallocate(v_re, v_im, part_c)
    end subroutine krylov_propagate_step

    !> \brief 对实向量 v 施加 exp(-i H dt): out = nrm * V_m exp(-i T_m dt) e_1
    subroutine lanczos_apply_expm(n, mm, v, nrm, dt, r_mult, out_c, err_est, m_eff, &
                                  st_lan, st_diag)
        integer, intent(in) :: n, mm
        real(dp), intent(in) :: v(n), nrm, dt
        interface
            subroutine r_mult(v_in, v_out)
                use mod_constants, only: dp
                real(dp), intent(in) :: v_in(:)
                real(dp), intent(out) :: v_out(:)
            end subroutine r_mult
        end interface
        complex(dp), intent(out) :: out_c(n)
        real(dp), intent(out) :: err_est
        integer, intent(out) :: m_eff, st_lan, st_diag

        real(dp), allocatable :: alpha(:), beta(:), v_basis(:, :)
        real(dp), allocatable :: tri(:, :), lam(:), umat(:, :)
        complex(dp), allocatable :: coeff(:)
        complex(dp) :: acc
        integer :: i, j

        allocate(alpha(mm), beta(mm + 1), v_basis(n, mm))
        call lanczos_tridiag(n, mm, v, r_mult, alpha, beta, v_basis, m_eff, st_lan)
        if (st_lan /= 0) then
            out_c = cmplx(v, 0.0_dp, kind=dp)
            err_est = 0.0_dp
            st_diag = 0
            deallocate(alpha, beta, v_basis)
            return
        end if

        ! 谱分解 T_m = U diag(lam) U^T (tql2):
        ! exp(-i T dt) e_1 = U [e^{-i lam_j dt} U(1, j)]_j
        allocate(tri(m_eff, m_eff), lam(m_eff), umat(m_eff, m_eff), coeff(m_eff))
        tri = 0.0_dp
        do i = 1, m_eff
            tri(i, i) = alpha(i)
            if (i < m_eff) then
                tri(i, i + 1) = beta(i + 1)
                tri(i + 1, i) = beta(i + 1)
            end if
        end do
        call diag_symmetric_matrix(m_eff, tri, lam, umat, st_diag)
        if (st_diag /= 0) then
            out_c = cmplx(v, 0.0_dp, kind=dp)
            err_est = 0.0_dp
            deallocate(alpha, beta, v_basis, tri, lam, umat, coeff)
            return
        end if

        do i = 1, m_eff
            coeff(i) = exp(cmplx(0.0_dp, -lam(i) * dt, kind=dp)) * umat(1, i)
        end do

        ! 后验误差估计: beta_{m+1} * |[exp(-iT dt)]_{m,1}|
        acc = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, m_eff
            acc = acc + umat(m_eff, i) * coeff(i)
        end do
        err_est = beta(m_eff + 1) * abs(acc)
        if (beta(m_eff + 1) <= 0.0_dp) err_est = 0.0_dp

        ! Krylov 系数 c = U * coeff, 则 exp(-iT dt) e_1 = c;
        ! out = nrm * V_m * c = nrm * sum_i c_i * v_basis(:, i)
        out_c = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, m_eff
            acc = cmplx(0.0_dp, 0.0_dp, kind=dp)
            do j = 1, m_eff
                acc = acc + umat(i, j) * coeff(j)
            end do
            out_c = out_c + acc * cmplx(v_basis(:, i), 0.0_dp, kind=dp)
        end do
        out_c = nrm * out_c

        deallocate(alpha, beta, v_basis, tri, lam, umat, coeff)
    end subroutine lanczos_apply_expm

end module mod_krylov_propagator
