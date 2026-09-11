!> \brief 离散变量表象（DVR）与空间网格核心算法模块
!> \details 包含 Colbert-Miller Sinc-DVR、Gauss-Legendre DVR 及 Fourier Grid Hamiltonian (FGH) 束缚态本征求解器。
!> \author LiHao
module mod_dvr_grid
    use mod_constants, only: dp, PI
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: dvr_1d_t
    public :: dvr_legendre_t
    public :: dvr_sinc_init
    public :: dvr_legendre_init
    public :: fgh_solve_bound_states

    !> \brief 一维 Sinc-DVR 网格对象
    type :: dvr_1d_t
        integer :: n_points = 0
        real(dp) :: x_min = 0.0_dp
        real(dp) :: x_max = 0.0_dp
        real(dp) :: dx = 0.0_dp
        real(dp) :: mass = 1.0_dp
        real(dp), allocatable :: x(:)          !< 离散坐标点向量 (n_points)
        real(dp), allocatable :: t_mat(:, :)   !< 动能矩阵 T_ij (n_points, n_points)
    end type dvr_1d_t

    !> \brief Gauss-Legendre DVR 角向网格对象
    type :: dvr_legendre_t
        integer :: n_points = 0
        real(dp), allocatable :: x(:)          !< Legendre 零点 cos(theta)_i
        real(dp), allocatable :: weights(:)    !< Gauss 权重 w_i
        real(dp), allocatable :: theta(:)      !< 对应极角 theta_i (rad)
        real(dp), allocatable :: j2_mat(:, :)  !< 角动量 J^2 算符矩阵
    end type dvr_legendre_t

contains

    !> \brief 初始化 Colbert-Miller Sinc-DVR（正弦离散变量表象）
    !> \param[in] x_min 坐标下界
    !> \param[in] x_max 坐标上界
    !> \param[in] n_pts 网格点数
    !> \param[in] mass 体系约化质量 (a.u.)
    !> \param[out] dvr 初始化的 Sinc-DVR 结构体
    subroutine dvr_sinc_init(x_min, x_max, n_pts, mass, dvr)
        real(dp), intent(in) :: x_min, x_max
        integer, intent(in) :: n_pts
        real(dp), intent(in) :: mass
        type(dvr_1d_t), intent(out) :: dvr

        integer :: i, j
        real(dp) :: factor, diff_sq

        dvr%n_points = n_pts
        dvr%x_min = x_min
        dvr%x_max = x_max
        dvr%dx = (x_max - x_min) / real(n_pts + 1, dp)
        dvr%mass = mass

        if (allocated(dvr%x)) deallocate(dvr%x)
        if (allocated(dvr%t_mat)) deallocate(dvr%t_mat)
        allocate(dvr%x(n_pts))
        allocate(dvr%t_mat(n_pts, n_pts))

        ! 生成等间距坐标
        do i = 1, n_pts
            dvr%x(i) = x_min + real(i, dp) * dvr%dx
        end do

        ! Colbert-Miller 解析动能矩阵元 T_ij = (hbar^2 / 2m*dx^2) * [...]
        factor = 1.0_dp / (2.0_dp * mass * dvr%dx**2)
        do i = 1, n_pts
            do j = 1, n_pts
                if (i == j) then
                    dvr%t_mat(i, j) = factor * (PI**2 / 3.0_dp)
                else
                    diff_sq = real((i - j)**2, dp)
                    if (mod(i - j, 2) == 0) then
                        dvr%t_mat(i, j) = factor * (2.0_dp / diff_sq)
                    else
                        dvr%t_mat(i, j) = -factor * (2.0_dp / diff_sq)
                    end if
                end if
            end do
        end do
    end subroutine dvr_sinc_init

    !> \brief 初始化 Gauss-Legendre DVR 角向网格
    subroutine dvr_legendre_init(n_pts, dvr)
        integer, intent(in) :: n_pts
        type(dvr_legendre_t), intent(out) :: dvr

        integer :: i, j, m
        real(dp) :: z, z1, p1, p2, p3, pp

        dvr%n_points = n_pts
        if (allocated(dvr%x)) deallocate(dvr%x)
        if (allocated(dvr%weights)) deallocate(dvr%weights)
        if (allocated(dvr%theta)) deallocate(dvr%theta)
        if (allocated(dvr%j2_mat)) deallocate(dvr%j2_mat)
        allocate(dvr%x(n_pts))
        allocate(dvr%weights(n_pts))
        allocate(dvr%theta(n_pts))
        allocate(dvr%j2_mat(n_pts, n_pts))

        ! Newton-Raphson 法求 Legendre 多项式零点与权重
        m = (n_pts + 1) / 2
        do i = 1, m
            z = cos(PI * (real(i, dp) - 0.25_dp) / (real(n_pts, dp) + 0.5_dp))
            find_root: do
                p1 = 1.0_dp
                p2 = 0.0_dp
                do j = 1, n_pts
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0_dp * real(j, dp) - 1.0_dp) * z * p2 - (real(j - 1, dp)) * p3) / real(j, dp)
                end do
                pp = real(n_pts, dp) * (z * p1 - p2) / (z * z - 1.0_dp)
                z1 = z
                z = z1 - p1 / pp
                if (abs(z - z1) < 1.0e-15_dp) exit find_root
            end do find_root

            dvr%x(i) = -z
            dvr%x(n_pts + 1 - i) = z
            dvr%weights(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp)
            dvr%weights(n_pts + 1 - i) = dvr%weights(i)
        end do

        do i = 1, n_pts
            dvr%theta(i) = acos(max(-1.0_dp, min(1.0_dp, dvr%x(i))))
        end do

        ! 构造角动量动能矩阵 J^2
        do i = 1, n_pts
            do j = 1, n_pts
                if (i == j) then
                    dvr%j2_mat(i, j) = real(n_pts * (n_pts + 1), dp) / 3.0_dp - &
                                       (1.0_dp - dvr%x(i)**2) / (1.0_dp - dvr%x(i)**2 + 1.0e-30_dp)
                else
                    dvr%j2_mat(i, j) = 2.0_dp / (dvr%x(i) - dvr%x(j))**2
                end if
            end do
        end do
    end subroutine dvr_legendre_init

    !> \brief 使用 Fourier Grid Hamiltonian (FGH) 求解任意分子势能面束缚态
    !> \param[in] dvr Sinc-DVR 结构体
    !> \param[in] v_pot 各网格点上的势能值 V(x_i) (a.u.)
    !> \param[out] eig_vals 本征能级数组 (n_points)
    !> \param[out] eig_vecs 本征波函数数组 (n_points, n_points)
    !> \param[out] stat 求解状态码
    subroutine fgh_solve_bound_states(dvr, v_pot, eig_vals, eig_vecs, stat)
        type(dvr_1d_t), intent(in) :: dvr
        real(dp), intent(in) :: v_pot(:)
        real(dp), intent(out) :: eig_vals(:)
        real(dp), intent(out) :: eig_vecs(:, :)
        integer, intent(out) :: stat

        integer :: n, i
        real(dp), allocatable :: h_mat(:, :)

        n = dvr%n_points
        allocate(h_mat(n, n))

        ! H_ij = T_ij + V(x_i) * delta_ij
        h_mat = dvr%t_mat
        do i = 1, n
            h_mat(i, i) = h_mat(i, i) + v_pot(i)
        end do

        ! 对角化实对称哈密顿量
        call diag_symmetric_matrix(n, h_mat, eig_vals, eig_vecs, stat)

        ! 波函数归一化由积分权 sqrt(dx) 保证: psi(x_i) = eig_vecs(i, v) / sqrt(dx)
        deallocate(h_mat)
    end subroutine fgh_solve_bound_states

end module mod_dvr_grid
