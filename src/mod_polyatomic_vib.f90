!> \brief 多原子分子 N 模简正振动: 乘积 DVR 变分与双线性耦合 (多原子拓展第三步)
!> \details 质量权重简正坐标下的多原子振动哈密顿量
!>          H = sum_k [ T_k (+) I ] + V(Q),
!>          V(Q) = sum_k (1/2) w_k^2 Q_k^2 + sum_{k<l} lam_kl Q_k Q_l,
!>          T_k 为各模的 Sinc-DVR 动能 (复用 mod_dvr_grid), 势能与双线性耦合
!>          在乘积格点上是对角的, 全哈密顿量经 Kronecker 积组装后用库内 QL
!>          对称求解器对角化。双模双线性耦合存在解析简正频率 (动力学矩阵
!>          本征值), 为变分结果提供机器精度的验收基准。
!> \note   本模块为多原子体系的入门步: 简正坐标 + 谐振 + 双线性力常数耦合。
!>          Morse 型非谐与一般力常数张量属后续增量 (见 12.2 路线图)。
!> \author LiHao
module mod_polyatomic_vib
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_dvr_grid, only: dvr_1d_t, dvr_sinc_init
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: coupled2_normal_freqs
    public :: poly_vib_product_levels
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
        call dvr_sinc_init(-x_half(2), x_half(2), n_pts, 1.0_dp, d2)
        call dvr_sinc_init(-x_half(3), x_half(3), n_pts, 1.0_dp, d3)

        n_tot = n_pts**n_modes
        allocate(h_mat(n_tot, n_tot), vecs(n_tot, n_tot))
        h_mat = 0.0_dp

        select case (n_modes)
        case (1)
            h_mat = d1%t_mat
            do i1 = 1, n_pts
                h_mat(i1, i1) = h_mat(i1, i1) &
                    + 0.5_dp * omegas(1)**2 * d1%x(i1)**2
            end do

        case (2)
            ! T1 (+) I
            do j2 = 1, n_pts
                do i1 = 1, n_pts
                    do i1p = 1, n_pts
                        h_mat(basis_idx(i1p, j2, 1, n_pts), basis_idx(i1, j2, 1, n_pts)) = &
                            h_mat(basis_idx(i1p, j2, 1, n_pts), basis_idx(i1, j2, 1, n_pts)) &
                            + d1%t_mat(i1p, i1)
                    end do
                end do
            end do
            ! I (+) T2
            do j2 = 1, n_pts
                do j2p = 1, n_pts
                    do i1 = 1, n_pts
                        h_mat(basis_idx(i1, j2p, 1, n_pts), basis_idx(i1, j2, 1, n_pts)) = &
                            h_mat(basis_idx(i1, j2p, 1, n_pts), basis_idx(i1, j2, 1, n_pts)) &
                            + d2%t_mat(j2p, j2)
                    end do
                end do
            end do
            ! 对角势 + 双线性耦合
            do j2 = 1, n_pts
                do i1 = 1, n_pts
                    pot = 0.5_dp * omegas(1)**2 * d1%x(i1)**2 &
                        + 0.5_dp * omegas(2)**2 * d2%x(j2)**2 &
                        + lam(1, 2) * d1%x(i1) * d2%x(j2)
                    h_mat(basis_idx(i1, j2, 1, n_pts), basis_idx(i1, j2, 1, n_pts)) = &
                        h_mat(basis_idx(i1, j2, 1, n_pts), basis_idx(i1, j2, 1, n_pts)) + pot
                end do
            end do

        case (3)
            ! T1 (+) I (+) I
            do j3 = 1, n_pts
                do j2 = 1, n_pts
                    do i1 = 1, n_pts
                        do i1p = 1, n_pts
                            h_mat(basis_idx(i1p, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                                h_mat(basis_idx(i1p, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) &
                                + d1%t_mat(i1p, i1)
                        end do
                    end do
                end do
            end do
            ! I (+) T2 (+) I
            do j3 = 1, n_pts
                do j2 = 1, n_pts
                    do j2p = 1, n_pts
                        do i1 = 1, n_pts
                            h_mat(basis_idx(i1, j2p, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                                h_mat(basis_idx(i1, j2p, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) &
                                + d2%t_mat(j2p, j2)
                        end do
                    end do
                end do
            end do
            ! I (+) I (+) T3
            do j3 = 1, n_pts
                do j3p = 1, n_pts
                    do j2 = 1, n_pts
                        do i1 = 1, n_pts
                            h_mat(basis_idx(i1, j2, j3p, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                                h_mat(basis_idx(i1, j2, j3p, n_pts), basis_idx(i1, j2, j3, n_pts)) &
                                + d3%t_mat(j3p, j3)
                        end do
                    end do
                end do
            end do
            ! 对角势 + 双线性耦合
            do j3 = 1, n_pts
                do j2 = 1, n_pts
                    do i1 = 1, n_pts
                        pot = 0.5_dp * omegas(1)**2 * d1%x(i1)**2 &
                            + 0.5_dp * omegas(2)**2 * d2%x(j2)**2 &
                            + 0.5_dp * omegas(3)**2 * d3%x(j3)**2 &
                            + lam(1, 2) * d1%x(i1) * d2%x(j2) &
                            + lam(1, 3) * d1%x(i1) * d3%x(j3) &
                            + lam(2, 3) * d2%x(j2) * d3%x(j3)
                        h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) = &
                            h_mat(basis_idx(i1, j2, j3, n_pts), basis_idx(i1, j2, j3, n_pts)) + pot
                    end do
                end do
            end do
        end select

        allocate(e_levels(n_tot))
        call diag_symmetric_matrix(n_tot, h_mat, e_levels, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine poly_vib_product_levels

end module mod_polyatomic_vib
