!> \brief 库仑三体系统、Hylleraas-Pekeris 坐标变换与两电子变分求解
!> \details 基于书中文献典藏第 42 节实现类氦两电子体系的高精度变分求解：
!>          - Hylleraas 显式关联坐标 (r1, r2, r12) 与 Pekeris 正交坐标 (u, v, w)
!>          - Rayleigh 梯度形式动能矩阵元（分部积分恒等式，仅需一阶导数）
!>          - 三维 Gauss-Laguerre 求积（张量积）计算全部矩阵元
!>          - Cholesky 对称化 + 对称本征求解器处理广义本征值问题 H c = E S c
!>          - 解析 Gamma 递推积分用于求积路径的独立交叉验证
!>          物理基准: He 1^1S 基态 E = -2.90372437... Hartree (Pekeris 1958)
!> \author LiHao
module mod_coulomb_threebody
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_constants, only: PI
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: perkeris_coordinate_transform
    public :: calc_perkeris_volume_element
    public :: calc_hylleraas_analytical_integral
    public :: gauss_laguerre_nodes
    public :: hylleraas_basis_indexing
    public :: solve_helium_ground_state_variational

    integer, parameter :: MAX_BASIS = 256
    integer, parameter :: MAX_QUAD  = 64

contains

    !> \brief 物理核间距 (r1, r2, r12) 与正交 Pekeris 坐标 (u, v, w) 双向解析变换
    !> \details u = r1 + r2 - r12,  v = r1 - r2 + r12,  w = -r1 + r2 + r12
    !>          r1 = (u+v)/2,     r2 = (u+w)/2,     r12 = (v+w)/2
    !>          三角形约束 |r1-r2| <= r12 <= r1+r2 解耦为 u, v, w >= 0
    pure subroutine perkeris_coordinate_transform(uvw, r_vals)
        real(dp), intent(in)  :: uvw(3)
        real(dp), intent(out) :: r_vals(3)

        r_vals(1) = 0.5_dp * (uvw(1) + uvw(2))   ! r1
        r_vals(2) = 0.5_dp * (uvw(1) + uvw(3))   ! r2
        r_vals(3) = 0.5_dp * (uvw(2) + uvw(3))   ! r12
    end subroutine perkeris_coordinate_transform

    !> \brief Pekeris 坐标变换 S 态约化雅可比体积元（不含角向常数）
    !> \details d^3 r1 d^3 r2 = (pi/16) (u+v)(u+w)(v+w) du dv dw；
    !>          常数因子 (pi/16) 在广义本征值问题 Hc = ESc 中严格相消，故可略去
    pure function calc_perkeris_volume_element(uvw) result(vol)
        real(dp), intent(in) :: uvw(3)
        real(dp) :: vol
        vol = (uvw(1) + uvw(2)) * (uvw(1) + uvw(3)) * (uvw(2) + uvw(3))
    end function calc_perkeris_volume_element

    !> \brief 解析 Hylleraas 单项式-指数积分 (Gamma 递推闭式)
    !> \details I(a,b,c; kappa) = int u^a v^b w^c e^{-kappa(u+v+w)} du dv dw
    !>          = a! b! c! / kappa^{a+b+c+3}，阶乘经对数累加避免溢出
    pure function calc_hylleraas_analytical_integral(a, b, c, kappa) result(val)
        integer, intent(in)  :: a, b, c
        real(dp), intent(in) :: kappa
        real(dp) :: val

        real(dp) :: log_fact
        integer  :: i, k

        log_fact = 0.0_dp
        do k = 0, 2
            do i = 2, max(a, b, c)
                if (k == 0 .and. i <= a) log_fact = log_fact + log(real(i, dp))
                if (k == 1 .and. i <= b) log_fact = log_fact + log(real(i, dp))
                if (k == 2 .and. i <= c) log_fact = log_fact + log(real(i, dp))
            end do
        end do
        val = exp(log_fact - real(a + b + c + 3, dp) * log(kappa))
    end function calc_hylleraas_analytical_integral

    !> \brief 生成 n 点标准 Gauss-Laguerre 求积节点与权重 (权函数 e^{-x})
    !> \details Golub-Welsch 方案: 由解析 Jacobi 矩阵特征分解直接得到正交多项式零点
    subroutine gauss_laguerre_nodes(n, x_nodes, w_nodes)
        integer, intent(in)   :: n
        real(dp), intent(out) :: x_nodes(n), w_nodes(n)

        real(dp) :: jac(n, n), eigvec(n, n)
        integer  :: i, stat

        jac = 0.0_dp
        do i = 1, n
            jac(i, i) = real(2 * i - 1, dp)
        end do
        do i = 1, n - 1
            jac(i, i + 1) = -real(i, dp)
            jac(i + 1, i) = -real(i, dp)
        end do

        call diag_symmetric_matrix(n, jac, x_nodes, eigvec, stat)
        do i = 1, n
            w_nodes(i) = eigvec(1, i)**2   ! mu0 = int e^{-x} dx = 1
        end do
    end subroutine gauss_laguerre_nodes

    !> \brief Hylleraas 基组索引: 枚举 |{l+m+n <= n_max}| 个基函数的 (l, m, n) 三元组
    pure subroutine hylleraas_basis_indexing(n_max, n_basis, idx_l, idx_m, idx_n)
        integer, intent(in)  :: n_max
        integer, intent(out) :: n_basis
        integer, intent(out) :: idx_l(MAX_BASIS), idx_m(MAX_BASIS), idx_n(MAX_BASIS)

        integer :: l, m, n, cnt

        cnt = 0
        do l = 0, n_max
            do m = 0, n_max - l
                do n = 0, n_max - l - m
                    if (cnt < MAX_BASIS) then
                        cnt = cnt + 1
                        idx_l(cnt) = l
                        idx_m(cnt) = m
                        idx_n(cnt) = n
                    end if
                end do
            end do
        end do
        n_basis = cnt
    end subroutine hylleraas_basis_indexing

    !> \brief 类氦两电子体系基态能量的 Hylleraas-Pekeris 变分求解
    !> \details 基组: phi_{lmn}(u,v,w) = u^l v^m w^n exp(-alpha(2u + v + w)/2)
    !>          其中 exp(-alpha(2u+v+w)/2) = exp(-alpha(r1+r2))。
    !>          动能采用 Rayleigh 梯度形式（分部积分恒等式，仅需一阶导数）:
    !>            T_ij = 1/2 int [grad_1 phi_i . grad_1 phi_j + grad_2 phi_i . grad_2 phi_j] dtau
    !>          S 态球坐标梯度展开:
    !>            grad_1 f = (df/dr1) r1_hat + (df/dr12) u12_hat,  其中
    !>            df/dr1  = f_u + f_v - f_w,   df/dr2 = f_u - f_v + f_w,   df/dr12 = -f_u + f_v + f_w
    !>            r1_hat . u12_hat = (r1^2 - r2^2 + r12^2) / (2 r1 r12)
    !>            r2_hat . u12_hat = (r2^2 - r1^2 + r12^2) / (2 r2 r12)
    !>          全部矩阵元由 3D Gauss-Laguerre 张量积求积计算，广义本征值问题
    !>          H c = E S c 经 Cholesky 分解对称化后调用对称本征求解器。
    subroutine solve_helium_ground_state_variational(n_max, z_charge, alpha, e_ground, stat)
        integer, intent(in)  :: n_max
        real(dp), intent(in) :: z_charge, alpha
        real(dp), intent(out):: e_ground
        integer, intent(out), optional :: stat

        integer  :: n_basis, i, j, k, stat_loc
        integer  :: idx_l(MAX_BASIS), idx_m(MAX_BASIS), idx_n(MAX_BASIS)
        integer  :: n_quad, ix, iy, iz
        real(dp), allocatable :: xq(:), wq(:)
        real(dp), allocatable :: hx(:, :), smat(:, :)
        real(dp), allocatable :: a_mat(:, :), eigs(:), eigvec(:, :)
        real(dp), allocatable :: chol(:, :), linv(:, :), tmp(:, :)
        real(dp), allocatable :: chi(:), du_fac(:), dv_fac(:), dw_fac(:)
        real(dp) :: uval, vval, wval, r1, r2, r12
        real(dp) :: fi, fj, dfi_u, dfi_v, dfi_w, dfj_u, dfj_v, dfj_w
        real(dp) :: g1i_r, g1j_r, g1i_u, g1j_u, g2i_r, g2j_r
        real(dp) :: a1, a2, w3d, acc_t, acc_s, acc_v

        if (present(stat)) stat = 0
        if (n_max < 0 .or. z_charge <= 0.0_dp .or. alpha <= 0.0_dp) then
            if (present(stat)) stat = -1
            e_ground = 0.0_dp
            return
        end if

        n_quad = min(MAX_QUAD, max(24, 6 * n_max))
        allocate(xq(n_quad), wq(n_quad))
        call gauss_laguerre_nodes(n_quad, xq, wq)
        call hylleraas_basis_indexing(n_max, n_basis, idx_l, idx_m, idx_n)
        allocate(hx(n_basis, n_basis), smat(n_basis, n_basis))
        allocate(a_mat(n_basis, n_basis), eigs(n_basis), eigvec(n_basis, n_basis))
        allocate(chol(n_basis, n_basis), linv(n_basis, n_basis), tmp(n_basis, n_basis))
        allocate(chi(n_basis), du_fac(n_basis), dv_fac(n_basis), dw_fac(n_basis))

        hx = 0.0_dp
        smat = 0.0_dp

        ! 坐标代换: x_u = 2*alpha*u, x_v = alpha*v, x_w = alpha*w —— 恰好使
        ! phi_i * phi_j 的总指数 e^{-2x_u} e^{-x_v} e^{-x_w} 与 GL 权重
        ! e^{-x_u} e^{-x_v} e^{-x_w} 逐轴完全吸收, 重叠积分对多项式严格精确
        do iz = 1, n_quad
            wval = xq(iz) / alpha
            do iy = 1, n_quad
                vval = xq(iy) / alpha
                do ix = 1, n_quad
                    uval = xq(ix) / (2.0_dp * alpha)

                    r1  = 0.5_dp * (uval + vval)
                    r2  = 0.5_dp * (uval + wval)
                    r12 = 0.5_dp * (vval + wval)

                    ! 每轴物理指数 e^{-2 alpha u}, e^{-alpha v}, e^{-alpha w} 恰好由
                    ! GL 权重 e^{-x_u}, e^{-x_v}, e^{-x_w} 完全提供（坐标代换见上），
                    ! 故基函数只以纯多项式 chi 进入被积函数；
                    ! 物理一阶导数 = (对数导数因子) * phi，因子不含指数
                    do k = 1, n_basis
                        chi(k) = uval**idx_l(k) * vval**idx_m(k) * wval**idx_n(k)
                        du_fac(k) = real(idx_l(k), dp) / uval - alpha
                        dv_fac(k) = real(idx_m(k), dp) / vval - 0.5_dp * alpha
                        dw_fac(k) = real(idx_n(k), dp) / wval - 0.5_dp * alpha
                    end do

                    w3d = (uval + vval) * (uval + wval) * (vval + wval) &
                          * wq(ix) * wq(iy) * wq(iz) / (2.0_dp * alpha**3)

                    a1 = (r1**2 - r2**2 + r12**2) / (2.0_dp * r1 * r12)
                    a2 = (r2**2 - r1**2 + r12**2) / (2.0_dp * r2 * r12)

                    do j = 1, n_basis
                        fj = chi(j)
                        dfj_u = du_fac(j); dfj_v = dv_fac(j); dfj_w = dw_fac(j)
                        do i = 1, n_basis
                            fi = chi(i)
                            dfi_u = du_fac(i); dfi_v = dv_fac(i); dfi_w = dw_fac(i)

                            acc_s = fi * fj

                            ! 物理坐标一阶导数因子（链式法则，逆映射为单位系数）
                            g1i_r = dfi_u + dfi_v - dfi_w
                            g1j_r = dfj_u + dfj_v - dfj_w
                            g1i_u = -dfi_u + dfi_v + dfi_w
                            g1j_u = -dfj_u + dfj_v + dfj_w
                            g2i_r = dfi_u - dfi_v + dfi_w
                            g2j_r = dfj_u - dfj_v + dfj_w

                            ! Rayleigh 梯度形式动能 (df/dr12 对两电子为同一函数导数, 复用 g1i_u/g1j_u)
                            ! 正交基完备性: A^2 f_u g_u + (1 - A^2) f_u g_u 严格相消为 f_u g_u
                            acc_t = g1i_r * g1j_r + a1 * (g1i_r * g1j_u + g1i_u * g1j_r) &
                                    + g1i_u * g1j_u &
                                    + g2i_r * g2j_r + a2 * (g2i_r * g1j_u + g1i_u * g2j_r) &
                                    + g1i_u * g1j_u

                            ! 动能被积函数 = 因子乘积 * chi_i * chi_j (指数由 GL 权重提供)
                            acc_t = acc_t * fi * fj
                            acc_v = -(z_charge / r1 + z_charge / r2 - 1.0_dp / r12) * fi * fj

                            hx(i, j) = hx(i, j) + w3d * (0.5_dp * acc_t + acc_v)
                            smat(i, j) = smat(i, j) + w3d * acc_s
                        end do
                    end do
                end do
            end do
        end do

        ! --- Cholesky 分解 S = L L^T ---
        chol = 0.0_dp
        do i = 1, n_basis
            acc_s = smat(i, i)
            do k = 1, i - 1
                acc_s = acc_s - chol(i, k)**2
            end do
            if (acc_s <= 0.0_dp) then
                if (present(stat)) stat = -2
                e_ground = 0.0_dp
                return
            end if
            chol(i, i) = sqrt(acc_s)
            do j = i + 1, n_basis
                acc_s = smat(j, i)
                do k = 1, i - 1
                    acc_s = acc_s - chol(j, k) * chol(i, k)
                end do
                chol(j, i) = acc_s / chol(i, i)
            end do
        end do

        ! --- L^{-1} 计算 (下三角回代) ---
        linv = 0.0_dp
        do i = 1, n_basis
            linv(i, i) = 1.0_dp / chol(i, i)
            do j = 1, i - 1
                acc_s = 0.0_dp
                do k = j, i - 1
                    acc_s = acc_s + chol(i, k) * linv(k, j)
                end do
                linv(i, j) = -acc_s / chol(i, i)
            end do
        end do

        ! --- A = L^{-1} H L^{-T}，对称化后调用对称本征求解器 ---
        tmp = matmul(hx, transpose(linv))
        a_mat = matmul(linv, tmp)
        a_mat = 0.5_dp * (a_mat + transpose(a_mat))

        call diag_symmetric_matrix(n_basis, a_mat, eigs, eigvec, stat_loc)
        if (stat_loc /= 0) then
            if (present(stat)) stat = -3
            e_ground = 0.0_dp
            return
        end if

        e_ground = eigs(1)
    end subroutine solve_helium_ground_state_variational

end module mod_coulomb_threebody
