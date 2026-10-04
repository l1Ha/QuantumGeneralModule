!> \brief 和积哈密顿量模块 (Sum-of-Products, SOP)
!> \details 基于书中第18.1节与第18.8节规范实现高维和积算符表示与矩阵自由作用：
!>          H = sum_alpha c_alpha prod_kappa h_alpha^(kappa)
!>          避免显式构造全张量积基稠密矩阵 (N_tot = prod n_kappa)，
!>          矩阵-向量作用计算复杂度为 O(M D n^(D+1))。
!> \author LiHao
module mod_sop_hamiltonian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: sop_term_1d_t
    public :: sop_hamiltonian_3d_t
    public :: init_sop_hamiltonian_3d
    public :: destroy_sop_hamiltonian_3d
    public :: sop_apply_3d

    !> \brief 一维局域算符因子 (矩阵形式 n x n)
    type :: sop_term_1d_t
        real(dp), allocatable :: mat(:, :)
    end type sop_term_1d_t

    !> \brief 三维和积哈密顿量结构体 (D = 3)
    type :: sop_hamiltonian_3d_t
        integer :: n1 = 0, n2 = 0, n3 = 0     !< 各维度网格点数
        integer :: n_terms = 0                !< 和积项数 M
        real(dp), allocatable :: coeffs(:)     !< 各项系数 c_alpha
        type(sop_term_1d_t), allocatable :: h1(:) !< 维度 1 局域算符 (大小 M)
        type(sop_term_1d_t), allocatable :: h2(:) !< 维度 2 局域算符 (大小 M)
        type(sop_term_1d_t), allocatable :: h3(:) !< 维度 3 局域算符 (大小 M)
    end type sop_hamiltonian_3d_t

contains

    !> \brief 初始化 3D SOP 哈密顿量
    subroutine init_sop_hamiltonian_3d(sop, n1, n2, n3, n_terms)
        type(sop_hamiltonian_3d_t), intent(out) :: sop
        integer, intent(in)                     :: n1, n2, n3, n_terms
        integer :: m

        sop%n1 = n1; sop%n2 = n2; sop%n3 = n3
        sop%n_terms = n_terms
        allocate(sop%coeffs(n_terms))
        allocate(sop%h1(n_terms), sop%h2(n_terms), sop%h3(n_terms))
        sop%coeffs = 0.0_dp

        do m = 1, n_terms
            allocate(sop%h1(m)%mat(n1, n1))
            allocate(sop%h2(m)%mat(n2, n2))
            allocate(sop%h3(m)%mat(n3, n3))
            sop%h1(m)%mat = 0.0_dp
            sop%h2(m)%mat = 0.0_dp
            sop%h3(m)%mat = 0.0_dp
        end do
    end subroutine init_sop_hamiltonian_3d

    !> \brief 销毁并释放 3D SOP 哈密顿量
    subroutine destroy_sop_hamiltonian_3d(sop)
        type(sop_hamiltonian_3d_t), intent(inout) :: sop
        integer :: m

        if (allocated(sop%coeffs)) deallocate(sop%coeffs)
        if (allocated(sop%h1)) then
            do m = 1, sop%n_terms
                if (allocated(sop%h1(m)%mat)) deallocate(sop%h1(m)%mat)
                if (allocated(sop%h2(m)%mat)) deallocate(sop%h2(m)%mat)
                if (allocated(sop%h3(m)%mat)) deallocate(sop%h3(m)%mat)
            end do
            deallocate(sop%h1, sop%h2, sop%h3)
        end if
        sop%n1 = 0; sop%n2 = 0; sop%n3 = 0
        sop%n_terms = 0
    end subroutine destroy_sop_hamiltonian_3d

    !> \brief 矩阵自由作用: out_vec = H * in_vec
    !> \details 采用逐维张量缩并，不显式展开全 Kronecker 矩阵
    subroutine sop_apply_3d(sop, in_vec, out_vec)
        type(sop_hamiltonian_3d_t), intent(in) :: sop
        real(dp), intent(in)                   :: in_vec(sop%n1, sop%n2, sop%n3)
        real(dp), intent(out)                  :: out_vec(sop%n1, sop%n2, sop%n3)

        real(dp) :: tmp1(sop%n1, sop%n2, sop%n3)
        real(dp) :: tmp2(sop%n1, sop%n2, sop%n3)
        real(dp) :: tmp3(sop%n1, sop%n2, sop%n3)
        integer  :: m, i1, i2, i3, k
        real(dp) :: c

        out_vec = 0.0_dp

        do m = 1, sop%n_terms
            c = sop%coeffs(m)
            if (abs(c) < 1.0e-30_dp) cycle

            ! 1. 作用维度 1: tmp1(i1, i2, i3) = sum_k h1(i1, k) * in_vec(k, i2, i3)
            tmp1 = 0.0_dp
            do i3 = 1, sop%n3
                do i2 = 1, sop%n2
                    do i1 = 1, sop%n1
                        do k = 1, sop%n1
                            tmp1(i1, i2, i3) = tmp1(i1, i2, i3) + sop%h1(m)%mat(i1, k) * in_vec(k, i2, i3)
                        end do
                    end do
                end do
            end do

            ! 2. 作用维度 2: tmp2(i1, i2, i3) = sum_k h2(i2, k) * tmp1(i1, k, i3)
            tmp2 = 0.0_dp
            do i3 = 1, sop%n3
                do i2 = 1, sop%n2
                    do k = 1, sop%n2
                        do i1 = 1, sop%n1
                            tmp2(i1, i2, i3) = tmp2(i1, i2, i3) + sop%h2(m)%mat(i2, k) * tmp1(i1, k, i3)
                        end do
                    end do
                end do
            end do

            ! 3. 作用维度 3: tmp3(i1, i2, i3) = sum_k h3(i3, k) * tmp2(i1, i2, k)
            tmp3 = 0.0_dp
            do i3 = 1, sop%n3
                do k = 1, sop%n3
                    do i2 = 1, sop%n2
                        do i1 = 1, sop%n1
                            tmp3(i1, i2, i3) = tmp3(i1, i2, i3) + sop%h3(m)%mat(i3, k) * tmp2(i1, i2, k)
                        end do
                    end do
                end do
            end do

            out_vec = out_vec + c * tmp3
        end do
    end subroutine sop_apply_3d

end module mod_sop_hamiltonian
