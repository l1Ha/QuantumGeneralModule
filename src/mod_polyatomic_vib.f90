!> \brief 多原子分子 N 模振动: 乘积 DVR 变分, 谐振/Morse 非谐与双线性耦合
!>         (多原子拓展第三、四步)
!> \details 质量权重简正坐标下的多原子振动哈密顿量
!>          H = sum_k [ T_k (+) I ] + V(Q),
!>          V(Q) = sum_k (1/2) w_k^2 Q_k^2            (poly_vib_product_levels, 谐振)
!>               | sum_k D_k (1-exp(-a_k Q_k))^2      (poly_vib_morse_levels,  Morse 非谐)
!>               + sum_{k<l} lam_kl Q_k Q_l,
!>          T_k 为各模的 Sinc-DVR 动能 (复用 mod_dvr_grid), 势能与双线性耦合
!>          在乘积格点上是对角的, 全哈密顿量经 Kronecker 和组装后用库内 QL
!>          对称求解器对角化。解析验收基准: 谐振极限 sum_k w_k(v_k+1/2) 与
!>          双模耦合简正频率 coupled2_normal_freqs; Morse 单模精确谱
!>          G(v) = w0(v+1/2) - w0^2 (v+1/2)^2/(4D) (morse_bound_energy);
!>          Coriolis 简并模对精确谱 E(v,l) = w(v+1) - g*l
!>          (poly_vib_coriolis_levels, 多原子拓展第五步);
!>          一般力场 (三阶+四阶力常数张量, 展开式系数)
!>          poly_vib_force_field_levels (多原子拓展第六步)。
!> \note   多原子拓展步序: (3) 谐振 + 双线性耦合; (4) Morse 型非谐力场;
!>          (5) Watson Coriolis 振动角动量 (简并模对 z 轴耦合可解核心);
!>          (6) 一般力场 (三阶/四阶力常数张量)。
!>          完整振转 Coriolis (J 耦合) 属后续增量。
!> \author LiHao
module mod_polyatomic_vib
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_dvr_grid, only: dvr_1d_t, dvr_sinc_init
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: coupled2_normal_freqs
    public :: poly_vib_product_levels
    public :: poly_vib_morse_levels
    public :: morse_bound_energy
    public :: morse_bound_count
    public :: poly_vib_coriolis_levels
    public :: poly_vib_force_field_levels
    integer, parameter, public :: POLY_MAX_MODES = 3

contains

    !> \brief 基矢线性索引: idx(i1, i2, i3) = ((i3-1) n + (i2-1)) n + i1 (1-based)
    pure integer function basis_idx(i1, i2, i3, n_pts) result(idx)
        integer, intent(in) :: i1, i2, i3, n_pts
        idx = ((i3 - 1) * n_pts + (i2 - 1)) * n_pts + i1
    end function basis_idx

    !> \brief 双模双线性耦合谐振子的解析简正频率
    !> \details 动力学矩阵 K = [[w1^2, lam], [lam, w2^2]] 的本征值开方:
    !>          Om^2 = ((w1^2 + w2^2) -/+ sqrt((w1^2 - w2^2)^2 + 4 lam^2)) / 2。
    !>          lam = 0 时退化为 Om = w1, w2。
    pure subroutine coupled2_normal_freqs(w1, w2, lam, om_lo, om_hi)
        real(dp), intent(in)  :: w1, w2, lam
        real(dp), intent(out) :: om_lo, om_hi
        real(dp) :: s, disc

        s = (w1 * w1 + w2 * w2) * 0.5_dp
        disc = sqrt(max((w1 * w1 - w2 * w2)**2 + 4.0_dp * lam * lam, 0.0_dp)) * 0.5_dp
        om_lo = sqrt(max(s - disc, 1.0e-12_dp))
        om_hi = sqrt(s + disc)
    end subroutine coupled2_normal_freqs

    !> \brief Morse 振子 (质量权重坐标) 的精确束缚能级 (从势阱底部计)
    !> \details V(Q) = D_e (1 - exp(-a Q))^2 (阱底为 0, 解离极限 D_e) 的精确谱
    !>          (因式分解/Poeschl-Teller):
    !>          G(v) = w0 (v+1/2) - w0^2 (v+1/2)^2 / (4 D_e),  w0 = a sqrt(2 D_e)。
    !>          二次项系数 w0^2/(4 D_e) 即 Dunham 展开的 w_e x_e (书中 7.10 节的
    !>          谐振微扰结果在这里严格化)。束缚能级数见 morse_bound_count。
    pure function morse_bound_energy(v, d_e, a_morse) result(e_v)
        real(dp), intent(in) :: v, d_e, a_morse
        real(dp) :: e_v, w0, vh

        w0 = a_morse * sqrt(2.0_dp * d_e)
        vh = v + 0.5_dp
        e_v = w0 * vh - w0 * w0 * vh * vh / (4.0_dp * d_e)
    end function morse_bound_energy

    !> \brief Morse 振子束缚能级总数: floor(sqrt(2 D_e)/a - 1/2) + 1
    pure function morse_bound_count(d_e, a_morse) result(n_bound)
        real(dp), intent(in) :: d_e, a_morse
        integer :: n_bound

        n_bound = int(floor(sqrt(2.0_dp * d_e) / a_morse - 0.5_dp)) + 1
    end function morse_bound_count

    !> \brief 向乘积哈密顿量累加第 k 模的 Kronecker 和动能 I⊗..⊗T_k⊗..⊗I
    !> \details 模 0 为最快指标 (basis_idx 约定); 低维时高阶指标循环自动收缩
    !>          为单次迭代 (basis_idx 把越界指标折叠到 1)。
    subroutine add_mode_kinetic(n_modes, k, n_pts, t_mat, h_mat)
        integer, intent(in) :: n_modes, k, n_pts
        real(dp), intent(in) :: t_mat(n_pts, n_pts)
        real(dp), intent(inout) :: h_mat(:, :)
        integer :: i, ip, j2, j3

        select case (k)
        case (1)
            do j3 = 1, merge(n_pts, 1, n_modes >= 3)
                do j2 = 1, merge(n_pts, 1, n_modes >= 2)
                    do i = 1, n_pts
                        do ip = 1, n_pts
                            h_mat(basis_idx(ip, j2, j3, n_pts), basis_idx(i, j2, j3, n_pts)) = &
                                h_mat(basis_idx(ip, j2, j3, n_pts), basis_idx(i, j2, j3, n_pts)) &
                                + t_mat(ip, i)
                        end do
                    end do
                end do
            end do
        case (2)
            do j3 = 1, merge(n_pts, 1, n_modes >= 3)
                do i = 1, n_pts
                    do j2 = 1, n_pts
                        do ip = 1, n_pts
                            h_mat(basis_idx(i, ip, j3, n_pts), basis_idx(i, j2, j3, n_pts)) = &
                                h_mat(basis_idx(i, ip, j3, n_pts), basis_idx(i, j2, j3, n_pts)) &
                                + t_mat(ip, j2)
                        end do
                    end do
                end do
            end do
        case (3)
            do i = 1, n_pts
                do j2 = 1, n_pts
                    do j3 = 1, n_pts
                        do ip = 1, n_pts
                            h_mat(basis_idx(i, j2, ip, n_pts), basis_idx(i, j2, j3, n_pts)) = &
                                h_mat(basis_idx(i, j2, ip, n_pts), basis_idx(i, j2, j3, n_pts)) &
                                + t_mat(ip, j3)
                        end do
                    end do
                end do
            end do
        end select
    end subroutine add_mode_kinetic

    !> \brief N 模 (N = 1..3) 谐振 + 双线性耦合的乘积 DVR 变分能级
    !> \details H = sum_k kron(T_k, I) + 对角势 1/2 sum_k w_k^2 x_k^2
    !>          + sum_{k<l} lam_kl x_k x_l (乘积格点上对角)。
    !>          x_k 在 [-x_half(k), +x_half(k)] 上取 n_pts 个 Sinc-DVR 点,
    !>          总基 = n_pts^N, 直接稠密对角化 (限 n_pts^N <= ~3400)。
    !>          验收: lambda = 0 时能级为 sum_k w_k (v_k + 1/2) (解析);
    !>          双模耦合时 ZPE = (Om_lo + Om_hi)/2 (解析, 见 coupled2_normal_freqs)。
    subroutine poly_vib_product_levels(n_modes, omegas, lam, x_half, n_pts, &
                                       e_levels, stat)
        integer, intent(in)  :: n_modes                !< 模数 1..3
        real(dp), intent(in) :: omegas(n_modes)        !< 各模谐振频率 (a.u.)
        real(dp), intent(in) :: lam(n_modes, n_modes)  !< 双线性耦合 (上三角有效)
        real(dp), intent(in) :: x_half(n_modes)        !< 各模半盒宽 (a.u.)
        integer, intent(in)  :: n_pts                  !< 每模格点数
        real(dp), intent(out), allocatable :: e_levels(:)
        integer, intent(out) :: stat

        type(dvr_1d_t) :: d1, d2, d3
        real(dp), allocatable :: h_mat(:, :), vecs(:, :)
        real(dp) :: pot
        integer  :: n_tot, i1, i1p, j2, j2p, j3, j3p
        integer  :: ierr

        stat = -1
        if (n_modes < 1 .or. n_modes > 3) return
        if (n_pts < 5 .or. n_pts**n_modes > 3400) return
        if (any(omegas(1:n_modes) <= 0.0_dp)) return
        if (any(x_half(1:n_modes) <= 0.0_dp)) return

        call dvr_sinc_init(-x_half(1), x_half(1), n_pts, 1.0_dp, d1)
        if (n_modes >= 2) call dvr_sinc_init(-x_half(2), x_half(2), n_pts, 1.0_dp, d2)
        if (n_modes >= 3) call dvr_sinc_init(-x_half(3), x_half(3), n_pts, 1.0_dp, d3)

        n_tot = n_pts**n_modes
        allocate(h_mat(n_tot, n_tot), vecs(n_tot, n_tot))
        h_mat = 0.0_dp

        ! Kronecker 和动能: 各模 T_k (+) I (公共助手, Morse 与谐振共用)
        call add_mode_kinetic(n_modes, 1, n_pts, d1%t_mat, h_mat)
        if (n_modes >= 2) call add_mode_kinetic(n_modes, 2, n_pts, d2%t_mat, h_mat)
        if (n_modes >= 3) call add_mode_kinetic(n_modes, 3, n_pts, d3%t_mat, h_mat)

        ! 对角势: 谐振 + 双线性耦合 (乘积格点上对角)
        do j3 = 1, merge(n_pts, 1, n_modes >= 3)
            do j2 = 1, merge(n_pts, 1, n_modes >= 2)
                do i1 = 1, n_pts
                    pot = 0.5_dp * omegas(1)**2 * d1%x(i1)**2
                    if (n_modes >= 2) then
                        pot = pot + 0.5_dp * omegas(2)**2 * d2%x(j2)**2 &
                            + lam(1, 2) * d1%x(i1) * d2%x(j2)
                    end if
                    if (n_modes >= 3) then
                        pot = pot + 0.5_dp * omegas(3)**2 * d3%x(j3)**2 &
                            + lam(1, 3) * d1%x(i1) * d3%x(j3) + lam(2, 3) * d2%x(j2) * d3%x(j3)
                    end if
                    h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                        h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) + pot
                end do
            end do
        end do

        allocate(e_levels(n_tot))
        call diag_symmetric_matrix(n_tot, h_mat, e_levels, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine poly_vib_product_levels

    !> \brief N 模 (N = 1..3) Morse 非谐 + 双线性耦合的乘积 DVR 变分能级
    !>         (多原子拓展第四步)
    !> \details H = sum_k kron(T_k, I) + 对角势 sum_k D_k (1 - exp(-a_k Q_k))^2
    !>          + sum_{k<l} lam_kl Q_k Q_l (乘积格点上对角)。
    !>          动能组装与谐振版本共用 add_mode_kinetic; Morse 非谐把能级间隔
    !>          从常数 w0 压缩为线性递减 (E_{v+1} - E_v = w0 - w0^2 (v+1)/2D_e)。
    !>          验收: 单模 lambda = 0 时 E_v 与精确谱 morse_bound_energy 对照,
    !>          束缚能级数 = morse_bound_count, 且 E_1 - E_0 < w0 (非谐红移);
    !>          多模 lambda = 0 时谱可分离, E(v1..vN) = sum_k G_k(v_k) (解析)。
    !> \warning 双线性耦合项与饱和 Morse 墙不相容: V 在 Q_k = -Q_l 方向
    !>          (Morse 墙饱和于 sum D_k 而 lam*Q_k*Q_l 无界增长) 无下界,
    !>          因此 lam != 0 仅适用于阱底附近的局域态与小盒; 接近解离的
    !>          态或大盒计算必须 lam = 0 (或改用完整四阶力场, 12.2 路线图)。
    !>          接近解离的 Morse 束缚态尾部衰减率 kappa = sqrt(2(D-E)) 饱和,
    !>          盒截断误差不随格点数收敛, 只随盒宽收敛 (见测试注释)。
    subroutine poly_vib_morse_levels(n_modes, d_e, a_morse, lam, x_half, n_pts, &
                                     e_levels, stat)
        integer, intent(in)  :: n_modes                !< 模数 1..3
        real(dp), intent(in) :: d_e(n_modes)           !< 各模 Morse 解离能 (a.u.)
        real(dp), intent(in) :: a_morse(n_modes)       !< 各模 Morse 范围参数 (a.u.^-1)
        real(dp), intent(in) :: lam(n_modes, n_modes)  !< 双线性耦合 (上三角有效)
        real(dp), intent(in) :: x_half(n_modes)        !< 各模半盒宽 (a.u.)
        integer, intent(in)  :: n_pts                  !< 每模格点数
        real(dp), intent(out), allocatable :: e_levels(:)
        integer, intent(out) :: stat

        type(dvr_1d_t) :: d1, d2, d3
        real(dp), allocatable :: h_mat(:, :), vecs(:, :)
        real(dp) :: pot
        integer  :: n_tot, i1, j2, j3, ierr

        stat = -1
        if (n_modes < 1 .or. n_modes > 3) return
        if (n_pts < 5 .or. n_pts**n_modes > 3400) return
        if (any(d_e(1:n_modes) <= 0.0_dp)) return
        if (any(a_morse(1:n_modes) <= 0.0_dp)) return
        if (any(x_half(1:n_modes) <= 0.0_dp)) return

        call dvr_sinc_init(-x_half(1), x_half(1), n_pts, 1.0_dp, d1)
        if (n_modes >= 2) call dvr_sinc_init(-x_half(2), x_half(2), n_pts, 1.0_dp, d2)
        if (n_modes >= 3) call dvr_sinc_init(-x_half(3), x_half(3), n_pts, 1.0_dp, d3)

        n_tot = n_pts**n_modes
        allocate(h_mat(n_tot, n_tot), vecs(n_tot, n_tot))
        h_mat = 0.0_dp

        call add_mode_kinetic(n_modes, 1, n_pts, d1%t_mat, h_mat)
        if (n_modes >= 2) call add_mode_kinetic(n_modes, 2, n_pts, d2%t_mat, h_mat)
        if (n_modes >= 3) call add_mode_kinetic(n_modes, 3, n_pts, d3%t_mat, h_mat)

        ! 对角势: Morse 非谐 + 双线性耦合
        do j3 = 1, merge(n_pts, 1, n_modes >= 3)
            do j2 = 1, merge(n_pts, 1, n_modes >= 2)
                do i1 = 1, n_pts
                    pot = d_e(1) * (1.0_dp - exp(-a_morse(1) * d1%x(i1)))**2
                    if (n_modes >= 2) then
                        pot = pot + d_e(2) * (1.0_dp - exp(-a_morse(2) * d2%x(j2)))**2 &
                            + lam(1, 2) * d1%x(i1) * d2%x(j2)
                    end if
                    if (n_modes >= 3) then
                        pot = pot + d_e(3) * (1.0_dp - exp(-a_morse(3) * d3%x(j3)))**2 &
                            + lam(1, 3) * d1%x(i1) * d3%x(j3) + lam(2, 3) * d2%x(j2) * d3%x(j3)
                    end if
                    h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                        h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) + pot
                end do
            end do
        end do

        allocate(e_levels(n_tot))
        call diag_symmetric_matrix(n_tot, h_mat, e_levels, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine poly_vib_morse_levels

    !> \brief 简并模对 + Coriolis 振动角动量耦合的乘积 DVR 变分能级
    !>        (多原子拓展第五步, Watson z 轴 zeta 耦合可解核心)
    !> \details 两个简并横模 (线性分子弯曲型) 在固定 J_z = 0 扇区内化为
    !>          2D 谐振子 + 振动角动量耦合:
    !>          H = T1(+)T2 + 1/2 w^2 (Q1^2 + Q2^2) - g (Q1 p2 - Q2 p1),
    !>          L_z = Q1 p2 - Q2 p1 为振动角动量 (Watson pi_alpha 的 Q x p
    !>          结构), g = 2B*zeta 为有效 Coriolis 常数。
    !>          阶梯算符推导: 圆形量子 a_+/- 下 L_z = n_+ - n_-, 精确谱
    !>          E(v, l) = w (v+1) - g l,  l = v, v-2, ..., -v
    !>          (v=1: Pi 双重态 2w -/+ g; v=2: Sigma + Delta 三重态 3w, 3w -/+ 2g;
    !>           +/- l 简并在本阶保持, 真 l 型倍频需高阶项)。
    !> \note   数值要点: (1) DVR 动量矩阵由对易恒等式 p = i m [T, X] 构造
    !>          ([T,X] 为实反对称, p 纯虚厄米); (2) Coriolis 项使哈密顿量
    !>          成为复厄米 H_c = H_r + i M (M 实反对称), 经 2 倍维实对称化
    !>          [[H_r, -M], [M, H_r]] 用库内实 QL 求解器对角化, 谱等于
    !>          H_c 与 H_c* 之并 —— 每个物理能级精确出现两次 (测试断言);
    !>          (3) ZPE (l=0) 与 g 无关 —— Coriolis 耦合不移动非简并态。
    subroutine poly_vib_coriolis_levels(omega, g_coriolis, x_half, n_pts, &
                                        e_levels, stat)
        real(dp), intent(in)  :: omega           !< 简并对频率 (a.u.)
        real(dp), intent(in)  :: g_coriolis      !< 有效 Coriolis 常数 g = 2B*zeta (a.u.)
        real(dp), intent(in)  :: x_half          !< 半盒宽 (a.u., 两模相同)
        integer, intent(in)   :: n_pts           !< 每模格点数 (2*n_pts^2 <= 4000)
        real(dp), intent(out), allocatable :: e_levels(:)
        integer, intent(out)  :: stat

        type(dvr_1d_t) :: d1, d2
        real(dp), allocatable :: h_r(:, :), m_mat(:, :), h2(:, :), vecs(:, :)
        real(dp), allocatable :: dx1(:, :), dx2(:, :)
        real(dp) :: pot
        integer  :: n_tot, n_dbl, i1, i2, j1, j2, ierr, k

        stat = -1
        if (n_pts < 5 .or. 2 * n_pts**2 > 4000) return
        if (omega <= 0.0_dp .or. x_half <= 0.0_dp) return

        call dvr_sinc_init(-x_half, x_half, n_pts, 1.0_dp, d1)
        call dvr_sinc_init(-x_half, x_half, n_pts, 1.0_dp, d2)

        n_tot = n_pts * n_pts
        n_dbl = 2 * n_tot
        allocate(h_r(n_tot, n_tot), m_mat(n_tot, n_tot), vecs(n_dbl, n_dbl))
        h_r = 0.0_dp
        m_mat = 0.0_dp

        ! --- 实对称部分: Kronecker 和动能 + 对角谐振势 --------------------
        call add_mode_kinetic(2, 1, n_pts, d1%t_mat, h_r)
        call add_mode_kinetic(2, 2, n_pts, d2%t_mat, h_r)
        do i2 = 1, n_pts
            do i1 = 1, n_pts
                pot = 0.5_dp * omega**2 * (d1%x(i1)**2 + d2%x(i2)**2)
                k = basis_idx(i1, i2, 1, n_pts)
                h_r(k, k) = h_r(k, k) + pot
            end do
        end do

        ! --- 反厄米动量矩阵 D = [T, X] (实反对称); p = i m D (m = 1) -------
        allocate(dx1(n_pts, n_pts), dx2(n_pts, n_pts))
        do j1 = 1, n_pts
            do i1 = 1, n_pts
                dx1(i1, j1) = d1%t_mat(i1, j1) * d1%x(j1) - d1%x(i1) * d1%t_mat(i1, j1)
                dx2(i1, j1) = d2%t_mat(i1, j1) * d2%x(j1) - d2%x(i1) * d2%t_mat(i1, j1)
            end do
        end do

        ! --- Coriolis: H_cor = -g (Q1 p2 - Q2 p1) = i M,                   ---
        !     M = -g (Q1 (x) D2 - D1 (x) Q2)  (实反对称);
        !     Q1/Q2 在各自模指标上是对角的 (delta 限制)
        do i2 = 1, n_pts
            do i1 = 1, n_pts
                k = basis_idx(i1, i2, 1, n_pts)
                ! Q1 (x) D2: j1 = i1 (Q1 对角), j2 自由
                do j2 = 1, n_pts
                    m_mat(k, basis_idx(i1, j2, 1, n_pts)) = &
                        m_mat(k, basis_idx(i1, j2, 1, n_pts)) &
                        - g_coriolis * d1%x(i1) * dx2(i2, j2)
                end do
                ! D1 (x) Q2: j2 = i2 (Q2 对角), j1 自由
                do j1 = 1, n_pts
                    m_mat(k, basis_idx(j1, i2, 1, n_pts)) = &
                        m_mat(k, basis_idx(j1, i2, 1, n_pts)) &
                        + g_coriolis * dx1(i1, j1) * d2%x(i2)
                end do
            end do
        end do

        ! --- 2 倍维实对称化: [[H_r, -M], [M, H_r]] -------------------------
        allocate(h2(n_dbl, n_dbl))
        h2 = 0.0_dp
        h2(1:n_tot, 1:n_tot) = h_r
        h2(n_tot + 1:n_dbl, n_tot + 1:n_dbl) = h_r
        h2(1:n_tot, n_tot + 1:n_dbl) = -m_mat
        h2(n_tot + 1:n_dbl, 1:n_tot) = m_mat

        allocate(e_levels(n_dbl))
        call diag_symmetric_matrix(n_dbl, h2, e_levels, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine poly_vib_coriolis_levels

    !> \brief N 模 (N = 1..3) 一般力场 (谐振 + 双线性 + 三阶/四阶力常数)
    !>        的乘积 DVR 变分能级 (多原子拓展第六步)
    !> \details V(Q) = 1/2 sum w_k^2 Q_k^2 + sum_{k<l} lam_kl Q_k Q_l
    !>          + sum_{i<=j<=k} c_ijk Q_i Q_j Q_k
    !>          + sum_{i<=j<=k<=l} d_ijkl Q_i Q_j Q_k Q_l,
    !>          c/d 为展开式系数 (非 Taylor 1/3!、1/4! 约定), 乘积格点上对角。
    !>          动能复用 add_mode_kinetic。
    !> \note   验收: (1) 零张量与第三步谐振机制逐位一致; (2) 单模纯四阶的
    !>          首阶微扰 E_v = w(v+1/2) + 3 d (2v^2+2v+1)/(4 w^2) (小 d);
    !>          (3) 纯三阶首阶修正为零 (<v|Q^3|v> = 0, 对称性); (4) 可分离
    !>          多模 (非对角张量元为零) E(v1..vN) = sum_k E_k。
    !> \warning 力场是平衡位形附近的局域展开: 大振幅/接近解离时全局行为
    !>          由高阶项主导, 应改用 Morse (第四步) 或全局 PES。
    subroutine poly_vib_force_field_levels(n_modes, omegas, lam, cubic, quartic, &
                                           x_half, n_pts, e_levels, stat)
        integer, intent(in)  :: n_modes
        real(dp), intent(in) :: omegas(n_modes)
        real(dp), intent(in) :: lam(n_modes, n_modes)
        real(dp), intent(in) :: cubic(n_modes, n_modes, n_modes)
        real(dp), intent(in) :: quartic(n_modes, n_modes, n_modes, n_modes)
        real(dp), intent(in) :: x_half(n_modes)
        integer, intent(in)  :: n_pts
        real(dp), intent(out), allocatable :: e_levels(:)
        integer, intent(out) :: stat

        type(dvr_1d_t) :: d1, d2, d3
        real(dp), allocatable :: h_mat(:, :), vecs(:, :)
        real(dp) :: pot, qq(3)
        integer  :: n_tot, i1, j2, j3, i, j, k, l, ierr

        stat = -1
        if (n_modes < 1 .or. n_modes > 3) return
        if (n_pts < 5 .or. n_pts**n_modes > 3400) return
        if (any(omegas(1:n_modes) <= 0.0_dp)) return
        if (any(x_half(1:n_modes) <= 0.0_dp)) return

        call dvr_sinc_init(-x_half(1), x_half(1), n_pts, 1.0_dp, d1)
        if (n_modes >= 2) call dvr_sinc_init(-x_half(2), x_half(2), n_pts, 1.0_dp, d2)
        if (n_modes >= 3) call dvr_sinc_init(-x_half(3), x_half(3), n_pts, 1.0_dp, d3)

        n_tot = n_pts**n_modes
        allocate(h_mat(n_tot, n_tot), vecs(n_tot, n_tot))
        h_mat = 0.0_dp

        call add_mode_kinetic(n_modes, 1, n_pts, d1%t_mat, h_mat)
        if (n_modes >= 2) call add_mode_kinetic(n_modes, 2, n_pts, d2%t_mat, h_mat)
        if (n_modes >= 3) call add_mode_kinetic(n_modes, 3, n_pts, d3%t_mat, h_mat)

        do j3 = 1, merge(n_pts, 1, n_modes >= 3)
            do j2 = 1, merge(n_pts, 1, n_modes >= 2)
                do i1 = 1, n_pts
                    qq(1) = d1%x(i1)
                    if (n_modes >= 2) qq(2) = d2%x(j2)
                    if (n_modes >= 3) qq(3) = d3%x(j3)
                    pot = 0.5_dp * omegas(1)**2 * qq(1)**2
                    if (n_modes >= 2) then
                        pot = pot + 0.5_dp * omegas(2)**2 * qq(2)**2 &
                            + lam(1, 2) * qq(1) * qq(2)
                    end if
                    if (n_modes >= 3) then
                        pot = pot + 0.5_dp * omegas(3)**2 * qq(3)**2 &
                            + lam(1, 3) * qq(1) * qq(3) + lam(2, 3) * qq(2) * qq(3)
                    end if
                    do i = 1, n_modes
                        do j = i, n_modes
                            do k = j, n_modes
                                pot = pot + cubic(i, j, k) * qq(i) * qq(j) * qq(k)
                            end do
                        end do
                    end do
                    do i = 1, n_modes
                        do j = i, n_modes
                            do k = j, n_modes
                                do l = k, n_modes
                                    pot = pot + quartic(i, j, k, l) &
                                        * qq(i) * qq(j) * qq(k) * qq(l)
                                end do
                            end do
                        end do
                    end do
                    h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                        h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) + pot
                end do
            end do
        end do

        allocate(e_levels(n_tot))
        call diag_symmetric_matrix(n_tot, h_mat, e_levels, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine poly_vib_force_field_levels

end module mod_polyatomic_vib
