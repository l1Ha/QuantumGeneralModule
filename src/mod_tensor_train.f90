!> \brief 张量列车与矩阵乘积态模块 (Tensor Train, TT / MPS)
!> \details 基于书中第18.3节与第18.8节规范实现高维张量 TT 格式核心操作：
!>          - 三核张量列车结构体 TT-3 (r0=1, r1, r2, r3=1)
!>          - 稠密张量向 TT 重构与还原
!>          - TT 范数与内积快速缩并 (O(D n r^3))
!> \author LiHao
module mod_tensor_train
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: tt_tensor_3d_t
    public :: init_tt_tensor_3d
    public :: destroy_tt_tensor_3d
    public :: tt_3d_to_dense
    public :: tt_3d_inner_product
    public :: tt_3d_norm

    !> \brief 三维张量列车类型
    !> \details 存储形式为三个核: G1(1, n1, r1), G2(r1, n2, r2), G3(r2, n3, 1)
    type :: tt_tensor_3d_t
        integer :: n1 = 0, n2 = 0, n3 = 0  !< 物理维度
        integer :: r1 = 0, r2 = 0          !< 内部键维 (Bond Dimensions)
        real(dp), allocatable :: g1(:, :, :) !< 核 1: (1, n1, r1)
        real(dp), allocatable :: g2(:, :, :) !< 核 2: (r1, n2, r2)
        real(dp), allocatable :: g3(:, :, :) !< 核 3: (r2, n3, 1)
    end type tt_tensor_3d_t

contains

    !> \brief 分配并初始化三维 TT 张量核
    subroutine init_tt_tensor_3d(tt, n1, n2, n3, r1, r2)
        type(tt_tensor_3d_t), intent(out) :: tt
        integer, intent(in)               :: n1, n2, n3, r1, r2

        tt%n1 = n1; tt%n2 = n2; tt%n3 = n3
        tt%r1 = r1; tt%r2 = r2

        allocate(tt%g1(1, n1, r1))
        allocate(tt%g2(r1, n2, r2))
        allocate(tt%g3(r2, n3, 1))

        tt%g1 = 0.0_dp
        tt%g2 = 0.0_dp
        tt%g3 = 0.0_dp
    end subroutine init_tt_tensor_3d

    !> \brief 销毁三维 TT 张量
    subroutine destroy_tt_tensor_3d(tt)
        type(tt_tensor_3d_t), intent(inout) :: tt
        if (allocated(tt%g1)) deallocate(tt%g1)
        if (allocated(tt%g2)) deallocate(tt%g2)
        if (allocated(tt%g3)) deallocate(tt%g3)
        tt%n1 = 0; tt%n2 = 0; tt%n3 = 0
        tt%r1 = 0; tt%r2 = 0
    end subroutine destroy_tt_tensor_3d

    !> \brief 将三维 TT 格式核缩并还原为稠密三阶张量
    !> \details C(i1, i2, i3) = sum_{a1, a2} G1(1, i1, a1) * G2(a1, i2, a2) * G3(a2, i3, 1)
    subroutine tt_3d_to_dense(tt, dense_out)
        type(tt_tensor_3d_t), intent(in) :: tt
        real(dp), intent(out)            :: dense_out(tt%n1, tt%n2, tt%n3)

        integer  :: i1, i2, i3, a1, a2
        real(dp) :: val

        dense_out = 0.0_dp
        ! OpenMP: 各 i3 切片独立写 dense_out, 无竞争 (服务器多核加速)
        !$omp parallel do private(i1, i2, val, a1, a2) schedule(static)
        do i3 = 1, tt%n3
            do i2 = 1, tt%n2
                do i1 = 1, tt%n1
                    val = 0.0_dp
                    do a2 = 1, tt%r2
                        do a1 = 1, tt%r1
                            val = val + tt%g1(1, i1, a1) * tt%g2(a1, i2, a2) * tt%g3(a2, i3, 1)
                        end do
                    end do
                    dense_out(i1, i2, i3) = val
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine tt_3d_to_dense

    !> \brief 高效计算两个 TT 张量的内积 <A, B> (无需展开全张量)
    pure function tt_3d_inner_product(a, b) result(dot_val)
        type(tt_tensor_3d_t), intent(in) :: a, b
        real(dp) :: dot_val

        real(dp) :: v1(a%r1, b%r1)
        real(dp) :: v2(a%r2, b%r2)
        integer  :: i, a1, b1, a2, b2

        dot_val = 0.0_dp
        if (a%n1 /= b%n1 .or. a%n2 /= b%n2 .or. a%n3 /= b%n3) return

        ! 第 1 维物理指标缩并: v1(a1, b1) = sum_i G1_A(1, i, a1) * G1_B(1, i, b1)
        v1 = 0.0_dp
        do b1 = 1, b%r1
            do a1 = 1, a%r1
                do i = 1, a%n1
                    v1(a1, b1) = v1(a1, b1) + a%g1(1, i, a1) * b%g1(1, i, b1)
                end do
            end do
        end do

        ! 第 2 维物理指标缩并: v2(a2, b2) = sum_{a1, b1, i} v1(a1, b1) * G2_A(a1, i, a2) * G2_B(b1, i, b2)
        v2 = 0.0_dp
        do b2 = 1, b%r2
            do a2 = 1, a%r2
                do i = 1, a%n2
                    do b1 = 1, b%r1
                        do a1 = 1, a%r1
                            v2(a2, b2) = v2(a2, b2) + v1(a1, b1) * a%g2(a1, i, a2) * b%g2(b1, i, b2)
                        end do
                    end do
                end do
            end do
        end do

        ! 第 3 维物理指标缩并: dot_val = sum_{a2, b2, i} v2(a2, b2) * G3_A(a2, i, 1) * G3_B(b2, i, 1)
        do i = 1, a%n3
            do b2 = 1, b%r2
                do a2 = 1, a%r2
                    dot_val = dot_val + v2(a2, b2) * a%g3(a2, i, 1) * b%g3(b2, i, 1)
                end do
            end do
        end do
    end function tt_3d_inner_product

    !> \brief 计算 TT 张量的 Frobenius 范数
    pure function tt_3d_norm(tt) result(nrm)
        type(tt_tensor_3d_t), intent(in) :: tt
        real(dp) :: nrm
        real(dp) :: dot_self
        dot_self = tt_3d_inner_product(tt, tt)
        nrm = sqrt(max(dot_self, 0.0_dp))
    end function tt_3d_norm

end module mod_tensor_train
