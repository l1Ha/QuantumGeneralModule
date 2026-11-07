!> \brief 多原子分子非对称陀螺转动哈密顿量 (Watson 刚转子极限)
!> \details 多原子体系拓展的第一步 (书中 3.4 节 Watson 哈密顿量的可解核心):
!>          刚性非对称陀螺 H = A Ja^2 + B Jb^2 + C Jc^2 (A >= B >= C) 在
!>          对称陀螺基 |J K M> (K 为 a 轴投影) 中的 Wang 矩阵构造与对角化。
!>          约定: H = (B+C)/2 * J^2 + (A - (B+C)/2) * Ja^2
!>                    + (B-C)/4 * (J_+^2 + J_-^2)
!>          极限自检: B=C 时退化为长对称陀螺 E = B J(J+1) + (A-B)K^2;
!>          A=B 时退化为扁对称陀螺 E = A J(J+1) - (A-C) K_c^2;
!>          A=B=C 球陀螺: E = A J(J+1), (2J+1) 重简并。
!> \author LiHao
module mod_watson_rotor
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    !> \brief 转动常数换算因子: B[cm^-1] = ROT_FACTOR / I[amu*angstrom^2]
    !> \details B = hbar^2 / (2 I) / (h c), 以 amu*angstrom^2 为惯量单位时
    !>          hbar^2/(2*amu*angstrom^2*hc) = 16.857631 cm^-1 (CODATA 2018)。
    real(dp), parameter, public :: ROT_FACTOR_AMU_A2_CM1 = 16.857631_dp

    public :: asymmetric_top_levels
    public :: asymmetric_top_levels_cd
    public :: inertia_constants_from_geometry
    public :: wang_matrix_element

contains

    !> \brief Wang 非对称陀螺矩阵元: 对角 (K) 与 Delta K = +2 耦合 (K -> K+2)
    pure function wang_matrix_element(a_const, b_const, c_const, j_qnum, k_bra, k_ket) result(h_val)
        real(dp), intent(in) :: a_const, b_const, c_const
        integer, intent(in)  :: j_qnum, k_bra, k_ket
        real(dp) :: h_val
        real(dp) :: jj, k_d, t1, t2
        integer :: k1, k2

        jj = real(j_qnum * (j_qnum + 1), dp)
        k1 = min(k_bra, k_ket)
        k2 = max(k_bra, k_ket)

        if (k2 - k1 == 0) then
            ! 对角元: (B+C)/2 * J(J+1) + (A - (B+C)/2) * K^2
            h_val = 0.5_dp * (b_const + c_const) * jj &
                    + (a_const - 0.5_dp * (b_const + c_const)) * real(k1 * k1, dp)
        else if (k2 - k1 == 2) then
            ! Delta K = +2: (B-C)/4 * sqrt([J(J+1)-K1(K1+1)] [J(J+1)-(K1+1)(K1+2)])
            k_d = real(k1 * (k1 + 1), dp)
            t1 = jj - k_d
            k_d = real((k1 + 1) * (k1 + 2), dp)
            t2 = jj - k_d
            h_val = 0.25_dp * (b_const - c_const) * sqrt(max(t1 * t2, 0.0_dp))
        else
            h_val = 0.0_dp
        end if
    end function wang_matrix_element

    !> \brief 给定 J 的非对称陀螺能级 (2J+1 条, 升序, cm^-1)
    !> \details 构造 (2J+1)x(2J+1) Wang 矩阵 (K = -J..J) 并用库内 QL 对称求解器
    !>          对角化。A >= B >= C 为惯例排序; 不满足时按值直接使用 (结果等价)。
    subroutine asymmetric_top_levels(a_const, b_const, c_const, j_qnum, energies, stat)
        real(dp), intent(in)  :: a_const, b_const, c_const
        integer, intent(in)   :: j_qnum
        real(dp), intent(out), allocatable :: energies(:)
        integer, intent(out)  :: stat

        integer  :: n_dim, k, kp, ierr
        real(dp), allocatable :: h_mat(:, :), vecs(:, :)

        if (j_qnum < 0 .or. a_const <= 0.0_dp .or. b_const <= 0.0_dp .or. c_const <= 0.0_dp) then
            stat = -1
            allocate(energies(0))
            return
        end if

        n_dim = 2 * j_qnum + 1
        allocate(h_mat(n_dim, n_dim), vecs(n_dim, n_dim), energies(n_dim))
        h_mat = 0.0_dp
        do k = -j_qnum, j_qnum
            do kp = k, j_qnum
                h_mat(kp + j_qnum + 1, k + j_qnum + 1) = &
                    wang_matrix_element(a_const, b_const, c_const, j_qnum, k, kp)
                h_mat(k + j_qnum + 1, kp + j_qnum + 1) = &
                    h_mat(kp + j_qnum + 1, k + j_qnum + 1)
            end do
        end do

        call diag_symmetric_matrix(n_dim, h_mat, energies, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine asymmetric_top_levels

    !> \brief 非对称陀螺能级 + 四阶离心畸变 (Watson A reduction 全项)
    !> \details 在刚转子 Wang 矩阵上叠加 A reduction 四次项:
    !>          H_cd = -D_J J^4 - D_JK J^2 Ja^2 - D_K Ja^4
    !>                 + d1 J^2 (J_+^2 + J_-^2) + d2 (J_+^4 + J_-^4)
    !>          (标准光谱学约定: D_J > 0 降低能级, E = B J(J+1) - D_J J^2(J+1)^2;
    !>           对角贡献: -D_J [J(J+1)]^2 - D_JK J(J+1) K^2 - D_K K^4)。
    !>          可选 d1/d2 (缺省为零) 为 Delta K = +/-2/4 的非对角项, 矩阵元
    !>          由角动量梯级给出:
    !>            <K+2|d1 J^2 J_+^2|K> = d1 [J(J+1)]^2
    !>                * sqrt([J(J+1)-K(K+1)][J(J+1)-(K+1)(K+2)])
    !>            <K+4|d2 J_+^4|K> = d2 * sqrt(prod_{m=0..3} [J(J+1)-(K+m)(K+m+1)])
    !>          (d1 与 Wang Delta K = 2 项共享同一梯级几何因子)。
    !> \note   迹定理: d1/d2 项纯非对角, Tr = 0, 故 sum(energies) 与 d1/d2
    !>          无关 (精确, 非微扰) —— Watson A/S reduction 在求和规则下
    !>          解耦, 测试组 7 断言。SPFIT 的 (2K+/-3) 对称化形式与本原基
    !>          梯级形式在对称化波函数基下等价, 符号约定对照属外部基准。
    !>          极限: D_JK = D_K = 0 且 A = B = C 时退化为
    !>          E = [A - D_J J(J+1)] J(J+1) (球陀螺离心畸变)。
    subroutine asymmetric_top_levels_cd(a_const, b_const, c_const, &
                                        d_j, d_jk, d_k, j_qnum, energies, stat, &
                                        d1_in, d2_in)
        real(dp), intent(in)  :: a_const, b_const, c_const
        real(dp), intent(in)  :: d_j, d_jk, d_k
        integer, intent(in)   :: j_qnum
        real(dp), intent(out), allocatable :: energies(:)
        real(dp), intent(in), optional :: d1_in, d2_in
        integer, intent(out)  :: stat

        integer  :: n_dim, k, kp, m, ierr
        real(dp) :: jj, kk, d1_v, d2_v, prod4
        real(dp), allocatable :: h_mat(:, :), vecs(:, :)

        if (j_qnum < 0 .or. a_const <= 0.0_dp .or. b_const <= 0.0_dp .or. c_const <= 0.0_dp) then
            stat = -1
            allocate(energies(0))
            return
        end if

        d1_v = 0.0_dp
        d2_v = 0.0_dp
        if (present(d1_in)) d1_v = d1_in
        if (present(d2_in)) d2_v = d2_in

        n_dim = 2 * j_qnum + 1
        allocate(h_mat(n_dim, n_dim), vecs(n_dim, n_dim), energies(n_dim))
        jj = real(j_qnum * (j_qnum + 1), dp)
        h_mat = 0.0_dp
        do k = -j_qnum, j_qnum
            kk = real(k * k, dp)
            h_mat(k + j_qnum + 1, k + j_qnum + 1) = &
                0.5_dp * (b_const + c_const) * jj &
                + (a_const - 0.5_dp * (b_const + c_const)) * kk &
                - d_j * jj * jj - d_jk * jj * kk - d_k * kk * kk
            do kp = k + 2, j_qnum, 2
                h_mat(kp + j_qnum + 1, k + j_qnum + 1) = &
                    wang_matrix_element(a_const, b_const, c_const, j_qnum, k, kp)
                h_mat(k + j_qnum + 1, kp + j_qnum + 1) = &
                    h_mat(kp + j_qnum + 1, k + j_qnum + 1)
            end do
            ! d1: Delta K = +/-2, 矩阵元 = d1 [J(J+1)]^2 * <K+2|J_+^2|K>
            if (d1_v /= 0.0_dp) then
                if (k + 2 <= j_qnum) then
                    kp = k + 2
                    prod4 = (jj - real(k * (k + 1), dp)) &
                            * (jj - real((k + 1) * (k + 2), dp))
                    h_mat(kp + j_qnum + 1, k + j_qnum + 1) = &
                        h_mat(kp + j_qnum + 1, k + j_qnum + 1) &
                        + d1_v * jj * jj * sqrt(max(prod4, 0.0_dp))
                    h_mat(k + j_qnum + 1, kp + j_qnum + 1) = &
                        h_mat(kp + j_qnum + 1, k + j_qnum + 1)
                end if
            end if
            ! d2: Delta K = +/-4, 矩阵元 = d2 * sqrt(prod_{m=0..3} [jj-(K+m)(K+m+1)])
            if (d2_v /= 0.0_dp) then
                if (k + 4 <= j_qnum) then
                    kp = k + 4
                    prod4 = 1.0_dp
                    do m = 0, 3
                        prod4 = prod4 * (jj - real((k + m) * (k + m + 1), dp))
                    end do
                    h_mat(kp + j_qnum + 1, k + j_qnum + 1) = &
                        h_mat(kp + j_qnum + 1, k + j_qnum + 1) &
                        + d2_v * sqrt(max(prod4, 0.0_dp))
                    h_mat(k + j_qnum + 1, kp + j_qnum + 1) = &
                        h_mat(kp + j_qnum + 1, k + j_qnum + 1)
                end if
            end if
        end do

        call diag_symmetric_matrix(n_dim, h_mat, energies, vecs, ierr)
        stat = merge(0, -2, ierr == 0)
    end subroutine asymmetric_top_levels_cd

    !> \brief 由核坐标与质量计算主惯量与转动常数 (多原子第一步: 几何 -> A,B,C)
    !> \details 以质心为原点的惯量张量 I_ab = sum_i m_i (|r_i|^2 delta_ab - r_a r_b),
    !>          对角化得主惯量矩 (升序 I_a <= I_b <= I_c), 转动常数 = ROT_FACTOR/I。
    !>          平面分子满足垂直轴定理 I_c = I_a + I_b (调用方可用于自检)。
    subroutine inertia_constants_from_geometry(positions, masses, a_const, b_const, c_const, &
                                               i_princ, com, stat)
        real(dp), intent(in)  :: positions(:, :)   !< (3, N) 坐标 (任意单位, 建议埃)
        real(dp), intent(in)  :: masses(:)         !< (N) 质量 (amu)
        real(dp), intent(out) :: a_const, b_const, c_const  !< 转动常数 (cm^-1)
        real(dp), intent(out) :: i_princ(3)        !< 主惯量矩 (amu*angstrom^2, 升序)
        real(dp), intent(out) :: com(3)            !< 质心坐标
        integer, intent(out)  :: stat

        integer  :: i, n_atoms, ierr
        real(dp) :: total_mass, r2
        real(dp) :: inertia(3, 3), d(3), vecs(3, 3)

        n_atoms = size(masses)
        stat = -1
        if (n_atoms < 2 .or. size(positions, 1) /= 3 .or. size(positions, 2) /= n_atoms) return
        if (any(masses <= 0.0_dp)) return

        total_mass = sum(masses)
        com = 0.0_dp
        do i = 1, n_atoms
            com = com + masses(i) * positions(:, i)
        end do
        com = com / total_mass

        inertia = 0.0_dp
        do i = 1, n_atoms
            r2 = sum((positions(:, i) - com)**2)
            inertia(1, 1) = inertia(1, 1) + masses(i) * (r2 - (positions(1, i) - com(1))**2)
            inertia(2, 2) = inertia(2, 2) + masses(i) * (r2 - (positions(2, i) - com(2))**2)
            inertia(3, 3) = inertia(3, 3) + masses(i) * (r2 - (positions(3, i) - com(3))**2)
            inertia(1, 2) = inertia(1, 2) - masses(i) * (positions(1, i) - com(1)) * (positions(2, i) - com(2))
            inertia(1, 3) = inertia(1, 3) - masses(i) * (positions(1, i) - com(1)) * (positions(3, i) - com(3))
            inertia(2, 3) = inertia(2, 3) - masses(i) * (positions(2, i) - com(2)) * (positions(3, i) - com(3))
        end do
        inertia(2, 1) = inertia(1, 2)
        inertia(3, 1) = inertia(1, 3)
        inertia(3, 2) = inertia(2, 3)

        call diag_symmetric_matrix(3, inertia, d, vecs, ierr)
        if (ierr /= 0) return

        i_princ = d
        if (minval(d) <= 0.0_dp) return
        a_const = ROT_FACTOR_AMU_A2_CM1 / d(1)
        b_const = ROT_FACTOR_AMU_A2_CM1 / d(2)
        c_const = ROT_FACTOR_AMU_A2_CM1 / d(3)
        stat = 0
    end subroutine inertia_constants_from_geometry

end module mod_watson_rotor
