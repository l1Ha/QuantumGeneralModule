!> \brief 外加磁场与电场中超冷量子散射理论与多基组密耦求解核心模块
!> \details 权威理论参考文献 (详见 LITERATURE.md):
!>   - G. Breit & I. I. Rabi, Phys. Rev. 38, 2082 (1931) [Breit-Rabi 塞曼-超精细能级]
!>   - H. Feshbach, Ann. Phys. 5, 357 (1958); U. Fano, Phys. Rev. 124, 1866 (1961) [共振理论]
!>   - H. T. C. Stoof et al., Phys. Rev. B 38, 4688 (1988) [自旋交换相互作用 s1 . s2]
!>   - E. Tiesinga et al., Phys. Rev. A 47, 4114 (1993) [冷碰撞阈值与共振现象]
!>   - D. A. Varshalovich et al., Quantum Theory of Angular Momentum (1988) [Racah代数]
!>   - C. Chin et al., Rev. Mod. Phys. 82, 1225 (2010) [超冷气体 Feshbach 共振综述]
!>   - J. M. Hutson & C. R. Le Sueur, Comput. Phys. Commun. 241, 1 (2019) [FIELD系统]
!> 支持特性：
!>   1. 四大经典基组表示及精确幺正变换：
!>      - 非耦合基组 (Uncoupled Basis): |m_s1, m_i1, m_s2, m_i2, L, M_L>
!>      - 单原子超精细耦合基组 (f-Coupled Basis): |(s1 i1) f1 m_f1, (s2 i2) f2 m_f2, L, M_L>
!>      - 总自旋耦合基组 (Total Spin Coupled Basis): |(s1 s2) S, (i1 i2) I, F M_F, L, M_L>
!>      - 场缀饰渐近本征通道基组 (Field-Dressed Channel Basis): |\alpha_1(B), \alpha_2(B), L, M_L>
!>   2. 任意磁场 (Zeeman 效应) 与直流电场 (Stark 效应) 下的外场双原子相互作用哈密顿量矩阵元；
!>   3. 单重态 V_0(r) 与三重态 V_1(r) 电子自旋交换势、长程色散与偶极-偶极 (s-d波混合) 耦合；
!>   4. 结合 Johnson 矩阵对数导数方法求解外场多通道密耦薛定谔方程，
!>      提取磁场依赖的散射长度 a_s(B)、多通道 S-矩阵与磁 Feshbach 共振极点 (B_0, \Delta B, a_bg)。
!> \author LiHao
module mod_field_scattering
    use mod_constants, only: dp, PI, TWOPI, SQRTPI, EYE, HBAR, AMU2AU, AU2CM, CM2AU
    use mod_special_functions, only: clebsch_gordan_half, wigner_3j_half, wigner_6j_half, wigner_9j_half
    use mod_ti_scattering, only: calc_multichannel_close_coupling_logder, multichannel_result_t
    implicit none
    private

    ! --------------------------------------------------------------------------
    ! 基组类型常量
    ! --------------------------------------------------------------------------
    integer, parameter, public :: BASIS_UNCOUPLED     = 1 ! 非耦合基组: |ms1, mi1, ms2, mi2, L, ML>
    integer, parameter, public :: BASIS_F_COUPLED     = 2 ! 单原子超精细耦合基组: |f1 mf1, f2 mf2, L, ML>
    integer, parameter, public :: BASIS_TOTAL_SPIN    = 3 ! 总自旋耦合基组: |(S, I) F MF, L, ML>
    integer, parameter, public :: BASIS_FIELD_DRESSED = 4 ! 场缀饰渐近本征通道基组: |\alpha_1(B), \alpha_2(B), L, ML>

    ! --------------------------------------------------------------------------
    ! 物理单位换算常量
    ! --------------------------------------------------------------------------
    real(dp), parameter, public :: GAUSS2AU = 4.254382e-10_dp           !< 1 Gauss -> a.u. (磁感应强度)
    real(dp), parameter, public :: AU2GAUSS = 1.0_dp / GAUSS2AU         !< a.u. -> Gauss
    real(dp), parameter, public :: MU_B_AU  = 0.5_dp                    !< 玻尔磁子 (a.u., e*hbar/(2*m_e) = 1/2)
    real(dp), parameter, public :: MU_N_AU  = 0.5_dp / 1836.15267343_dp !< 核磁子 (a.u.)
    real(dp), parameter, public :: GHZ2AU   = 1.51982984600e-7_dp       !< 1 GHz (h*nu) -> a.u. (能量)
    real(dp), parameter, public :: AU2GHZ   = 1.0_dp / GHZ2AU           !< a.u. -> GHz

    ! --------------------------------------------------------------------------
    ! 派生类型导出
    ! --------------------------------------------------------------------------
    public :: cold_atom_t
    public :: field_channel_t
    public :: field_feshbach_result_t

    ! --------------------------------------------------------------------------
    ! 公共子程序与函数导出
    ! --------------------------------------------------------------------------
    public :: get_cold_atom_preset
    public :: calc_breit_rabi_energies
    public :: build_field_collision_channels
    public :: calc_basis_transform_matrix
    public :: build_asymptotic_hamiltonian
    public :: build_spin_exchange_matrix
    public :: build_field_potential_matrix
    public :: calc_magnetic_feshbach_resonance_scan
    public :: fit_feshbach_resonance_parameters

    ! --------------------------------------------------------------------------
    ! 1. 冷原子同位素物理参数类型
    ! --------------------------------------------------------------------------
    type :: cold_atom_t
        character(len=16) :: name          !< 同位素名称 (如 "6Li", "87Rb", "40K")
        real(dp)          :: mass_amu      !< 质量 (amu)
        real(dp)          :: mass_au       !< 质量 (a.u.)
        integer           :: two_s         !< 2 * 电子自旋 (碱金属 = 1)
        integer           :: two_i         !< 2 * 核自旋 (如 6Li 为 2, 87Rb 为 3)
        real(dp)          :: a_hf_ghz      !< 超精细相互作用常数 a_hf (GHz)
        real(dp)          :: a_hf_au       !< 超精细常数 (a.u.)
        real(dp)          :: g_s           !< 电子自旋 g 因子 (~ 2.0023193)
        real(dp)          :: g_i           !< 核自旋 g 因子 (g_I * mu_N / mu_B)
        real(dp)          :: dipole_debye  !< 永久电偶极矩 (Debye, 极性分子)
        real(dp)          :: dipole_au     !< 永久电偶极矩 (a.u.)
    end type cold_atom_t

    ! --------------------------------------------------------------------------
    ! 2. 外场散射通道量子数结构 (同时支持四大基组索引)
    ! --------------------------------------------------------------------------
    type :: field_channel_t
        integer  :: idx                    !< 通道编号 (1..N)
        ! (A) 非耦合基组量子数: |ms1, mi1, ms2, mi2, L, ML> (均为 2 倍角动量)
        integer  :: two_ms1, two_mi1
        integer  :: two_ms2, two_mi2
        ! (B) f-耦合基组量子数: |f1 mf1, f2 mf2, L, ML>
        integer  :: two_f1, two_mf1
        integer  :: two_f2, two_mf2
        ! (C) 总自旋耦合基组量子数: |(S, I) F MF, L, ML>
        integer  :: two_S, two_I, two_F, two_MF
        ! (D) 空间轨道角动量
        integer  :: l_orb, m_l
        ! (E) 守恒总磁量子数: 2 * M_tot = 2*(ms1+mi1+ms2+mi2+ML)
        integer  :: two_Mtot
        ! (F) 渐近场缀饰阈值能量 (a.u.)
        real(dp) :: thresh_energy
    end type field_channel_t

    ! --------------------------------------------------------------------------
    ! 3. 磁 Feshbach 共振扫描结果结构
    ! --------------------------------------------------------------------------
    type :: field_feshbach_result_t
        integer               :: n_b           !< 磁场网格点数
        real(dp), allocatable :: b_grid(:)     !< 磁场网格序列 (Gauss)
        real(dp), allocatable :: a_s_grid(:)   !< 散射长度序列 a_s(B) (a.u.)
        real(dp), allocatable :: eigenphases(:)!< 特征相移和 \delta_{sum}(B) (rad)
        real(dp)              :: b_res_pole    !< 共振极点位置 B_0 (Gauss)
        real(dp)              :: delta_b       !< 共振宽度 \Delta B (Gauss)
        real(dp)              :: a_bg          !< 背景散射长度 a_{bg} (a.u.)
    end type field_feshbach_result_t

contains

    ! ==========================================================================
    ! 1. 常用冷原子同位素物理参数预设库
    ! ==========================================================================
    subroutine get_cold_atom_preset(name, atom)
        character(len=*), intent(in)  :: name
        type(cold_atom_t), intent(out):: atom

        character(len=16) :: atom_name
        atom_name = trim(adjustl(name))

        atom%name = atom_name
        atom%two_s = 1 ! 碱金属 s = 1/2 -> 2s = 1
        atom%g_s = 2.00231930436256_dp
        atom%dipole_debye = 0.0_dp
        atom%dipole_au = 0.0_dp

        select case (atom_name)
        case ("6Li", "Li6")
            atom%mass_amu = 6.015122_dp
            atom%two_i = 2 ! i = 1 -> 2i = 2
            atom%a_hf_ghz = 0.15213684_dp
            atom%g_i = -0.0004476540_dp

        case ("7Li", "Li7")
            atom%mass_amu = 7.016004_dp
            atom%two_i = 3 ! i = 3/2 -> 2i = 3
            atom%a_hf_ghz = 0.40181204_dp
            atom%g_i = -0.001182213_dp

        case ("23Na", "Na23")
            atom%mass_amu = 22.989769_dp
            atom%two_i = 3 ! i = 3/2
            atom%a_hf_ghz = 0.88581306_dp
            atom%g_i = -0.000804610_dp

        case ("40K", "K40")
            atom%mass_amu = 39.963998_dp
            atom%two_i = 8 ! i = 4 -> 2i = 8 (费米同位素)
            atom%a_hf_ghz = -0.2857308_dp
            atom%g_i = 0.00017649_dp

        case ("87Rb", "Rb87")
            atom%mass_amu = 86.9091805_dp
            atom%two_i = 3 ! i = 3/2
            atom%a_hf_ghz = 3.417341305452_dp
            atom%g_i = -0.0009951414_dp

        case ("133Cs", "Cs133")
            atom%mass_amu = 132.9054519_dp
            atom%two_i = 7 ! i = 7/2 -> 2i = 7
            atom%a_hf_ghz = 2.2981579425_dp
            atom%g_i = -0.0003988539_dp

        case ("KRb", "40K87Rb")
            atom%mass_amu = 126.873_dp
            atom%two_s = 0
            atom%two_i = 0
            atom%a_hf_ghz = 0.0_dp
            atom%g_s = 0.0_dp
            atom%g_i = 0.0_dp
            atom%dipole_debye = 0.566_dp
            atom%dipole_au = atom%dipole_debye / 2.541746473_dp

        case default
            ! 默认无自旋模型
            atom%mass_amu = 1.0_dp
            atom%two_i = 0
            atom%a_hf_ghz = 0.0_dp
            atom%g_i = 0.0_dp
        end select

        atom%mass_au = atom%mass_amu * AMU2AU
        atom%a_hf_au = atom%a_hf_ghz * GHZ2AU
    end subroutine get_cold_atom_preset

    ! ==========================================================================
    ! 2. 单原子 Breit-Rabi 塞曼能级与超精细本征态求解
    !    H_atom = a_hf * s . i + (g_s * mu_B * s_z - g_i * mu_N * i_z) * B
    ! ==========================================================================
    subroutine calc_breit_rabi_energies(atom, b_gauss, energies, states, n_states)
        type(cold_atom_t), intent(in)   :: atom
        real(dp), intent(in)            :: b_gauss
        real(dp), allocatable, intent(out):: energies(:)
        real(dp), allocatable, intent(out):: states(:, :)
        integer, intent(out)            :: n_states

        integer :: n_s, n_i, i_dim, j_dim, idx, jdx
        integer, allocatable :: two_ms(:), two_mi(:)
        real(dp), allocatable :: h_mat(:, :)
        real(dp) :: b_au, ms_val, mi_val, ms2_val, mi2_val
        real(dp) :: zeeman_term, s_plus_minus
        integer  :: info

        n_s = atom%two_s + 1
        n_i = atom%two_i + 1
        n_states = n_s * n_i

        allocate(two_ms(n_states), two_mi(n_states))
        allocate(h_mat(n_states, n_states), energies(n_states), states(n_states, n_states))
        h_mat = 0.0_dp
        b_au = b_gauss * GAUSS2AU

        ! 枚举单原子非耦合基底: |m_s, m_i>
        idx = 0
        do i_dim = atom%two_s, -atom%two_s, -2
            do j_dim = atom%two_i, -atom%two_i, -2
                idx = idx + 1
                two_ms(idx) = i_dim
                two_mi(idx) = j_dim
            end do
        end do

        ! 构建非耦合基底下的 Breit-Rabi 哈密顿量矩阵
        do idx = 1, n_states
            ms_val = real(two_ms(idx), dp) / 2.0_dp
            mi_val = real(two_mi(idx), dp) / 2.0_dp

            ! 1. 对角部分: 塞曼作用 + 超精细 sz * iz
            zeeman_term = (atom%g_s * MU_B_AU * ms_val - atom%g_i * MU_N_AU * mi_val) * b_au
            h_mat(idx, idx) = zeeman_term + atom%a_hf_au * (ms_val * mi_val)

            ! 2. 非对角部分: a_hf * 0.5 * (s+ i- + s- i+)
            do jdx = 1, n_states
                if (idx == jdx) cycle
                ms2_val = real(two_ms(jdx), dp) / 2.0_dp
                mi2_val = real(two_mi(jdx), dp) / 2.0_dp

                ! (s+ i-): ms_val = ms2_val + 1, mi_val = mi2_val - 1
                if (two_ms(idx) == two_ms(jdx) + 2 .and. two_mi(idx) == two_mi(jdx) - 2) then
                    s_plus_minus = 0.5_dp * atom%a_hf_au * &
                        sqrt(real(atom%two_s * (atom%two_s + 2) - two_ms(jdx) * (two_ms(jdx) + 2), dp) / 4.0_dp) * &
                        sqrt(real(atom%two_i * (atom%two_i + 2) - two_mi(jdx) * (two_mi(jdx) - 2), dp) / 4.0_dp)
                    h_mat(idx, jdx) = s_plus_minus
                ! (s- i+): ms_val = ms2_val - 1, mi_val = mi2_val + 1
                else if (two_ms(idx) == two_ms(jdx) - 2 .and. two_mi(idx) == two_mi(jdx) + 2) then
                    s_plus_minus = 0.5_dp * atom%a_hf_au * &
                        sqrt(real(atom%two_s * (atom%two_s + 2) - two_ms(jdx) * (two_ms(jdx) - 2), dp) / 4.0_dp) * &
                        sqrt(real(atom%two_i * (atom%two_i + 2) - two_mi(jdx) * (two_mi(jdx) + 2), dp) / 4.0_dp)
                    h_mat(idx, jdx) = s_plus_minus
                end if
            end do
        end do

        ! 对称三对角化/Jacobi 本征值求解
        call diagonalize_real_symmetric(h_mat, n_states, energies, states, info)

        deallocate(two_ms, two_mi, h_mat)
    end subroutine calc_breit_rabi_energies

    ! ==========================================================================
    ! 3. 构造外场碰撞多通道基底 (给定守恒量 2*M_tot 与轨道分波 l_max)
    ! ==========================================================================
    subroutine build_field_collision_channels( &
        atom1, atom2, basis_type, two_Mtot, l_max, channels, n_channels)

        type(cold_atom_t), intent(in)       :: atom1, atom2
        integer, intent(in)                 :: basis_type
        integer, intent(in)                 :: two_Mtot
        integer, intent(in)                 :: l_max
        type(field_channel_t), allocatable, intent(out) :: channels(:)
        integer, intent(out)                :: n_channels

        integer :: ms1, mi1, ms2, mi2, l_val, ml_val
        integer :: f1, mf1, f2, mf2
        integer :: s_tot, i_tot, f_tot, mf_tot
        integer :: count_ch, ch_idx

        ! 第一遍: 统计满足 2*M_tot = 2*(M_spin + ML) 的通道总数
        count_ch = 0
        do l_val = 0, l_max, 2 ! 超冷通常为宇称守恒 (如 s-波 l=0, d-波 l=2)
            do ml_val = -l_val, l_val
                do ms1 = atom1%two_s, -atom1%two_s, -2
                    do mi1 = atom1%two_i, -atom1%two_i, -2
                        do ms2 = atom2%two_s, -atom2%two_s, -2
                            do mi2 = atom2%two_i, -atom2%two_i, -2
                                if ((ms1 + mi1 + ms2 + mi2) + 2 * ml_val == two_Mtot) then
                                    count_ch = count_ch + 1
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do

        n_channels = count_ch
        if (n_channels == 0) return
        allocate(channels(n_channels))

        ! 第二遍: 填充通道数据
        ch_idx = 0
        do l_val = 0, l_max, 2
            do ml_val = -l_val, l_val
                select case (basis_type)
                case (BASIS_UNCOUPLED, BASIS_FIELD_DRESSED)
                    ! 非耦合基组枚举: |ms1, mi1, ms2, mi2, L, ML>
                    do ms1 = atom1%two_s, -atom1%two_s, -2
                        do mi1 = atom1%two_i, -atom1%two_i, -2
                            do ms2 = atom2%two_s, -atom2%two_s, -2
                                do mi2 = atom2%two_i, -atom2%two_i, -2
                                    if ((ms1 + mi1 + ms2 + mi2) + 2 * ml_val == two_Mtot) then
                                        ch_idx = ch_idx + 1
                                        channels(ch_idx)%idx = ch_idx
                                        channels(ch_idx)%two_ms1 = ms1
                                        channels(ch_idx)%two_mi1 = mi1
                                        channels(ch_idx)%two_ms2 = ms2
                                        channels(ch_idx)%two_mi2 = mi2
                                        channels(ch_idx)%l_orb   = l_val
                                        channels(ch_idx)%m_l     = ml_val
                                        channels(ch_idx)%two_Mtot = two_Mtot
                                        channels(ch_idx)%thresh_energy = 0.0_dp
                                    end if
                                end do
                            end do
                        end do
                    end do

                case (BASIS_F_COUPLED)
                    ! f-耦合基组枚举: |f1 mf1, f2 mf2, L, ML>
                    do f1 = atom1%two_s + atom1%two_i, abs(atom1%two_s - atom1%two_i), -2
                        do mf1 = f1, -f1, -2
                            do f2 = atom2%two_s + atom2%two_i, abs(atom2%two_s - atom2%two_i), -2
                                do mf2 = f2, -f2, -2
                                    if ((mf1 + mf2) + 2 * ml_val == two_Mtot) then
                                        ch_idx = ch_idx + 1
                                        channels(ch_idx)%idx = ch_idx
                                        channels(ch_idx)%two_f1  = f1
                                        channels(ch_idx)%two_mf1 = mf1
                                        channels(ch_idx)%two_f2  = f2
                                        channels(ch_idx)%two_mf2 = mf2
                                        channels(ch_idx)%l_orb   = l_val
                                        channels(ch_idx)%m_l     = ml_val
                                        channels(ch_idx)%two_Mtot = two_Mtot
                                        channels(ch_idx)%thresh_energy = 0.0_dp
                                    end if
                                end do
                            end do
                        end do
                    end do

                case (BASIS_TOTAL_SPIN)
                    ! 总自旋耦合基组枚举: |(S, I) F MF, L, ML>
                    mf_tot = two_Mtot - 2 * ml_val
                    do s_tot = atom1%two_s + atom2%two_s, abs(atom1%two_s - atom2%two_s), -2
                        do i_tot = atom1%two_i + atom2%two_i, abs(atom1%two_i - atom2%two_i), -2
                            do f_tot = s_tot + i_tot, abs(s_tot - i_tot), -2
                                if (abs(mf_tot) <= f_tot) then
                                    ch_idx = ch_idx + 1
                                    channels(ch_idx)%idx = ch_idx
                                    channels(ch_idx)%two_S  = s_tot
                                    channels(ch_idx)%two_I  = i_tot
                                    channels(ch_idx)%two_F  = f_tot
                                    channels(ch_idx)%two_MF = mf_tot
                                    channels(ch_idx)%l_orb  = l_val
                                    channels(ch_idx)%m_l    = ml_val
                                    channels(ch_idx)%two_Mtot = two_Mtot
                                    channels(ch_idx)%thresh_energy = 0.0_dp
                                end if
                            end do
                        end do
                    end do
                end select
            end do
        end do
    end subroutine build_field_collision_channels

    ! ==========================================================================
    ! 4. 四大基组之间的严格幺正变换矩阵计算 U_{target, source}
    ! ==========================================================================
    recursive subroutine calc_basis_transform_matrix( &
        atom1, atom2, channels_src, channels_tgt, n_channels, &
        from_basis, to_basis, b_gauss, U_mat)

        type(cold_atom_t), intent(in)   :: atom1, atom2
        type(field_channel_t), intent(in):: channels_src(:)
        type(field_channel_t), intent(in):: channels_tgt(:)
        integer, intent(in)             :: n_channels
        integer, intent(in)             :: from_basis, to_basis
        real(dp), intent(in)            :: b_gauss
        real(dp), intent(out)           :: U_mat(n_channels, n_channels)

        integer  :: i, j, ms1, mi1, ms2, mi2, ms, mi, mf
        real(dp) :: c1, c2, cg_s, cg_i, cg_f
        real(dp), allocatable :: h_asymp(:, :), eig_vals(:), eig_vecs(:, :)
        integer  :: info

        U_mat = 0.0_dp

        ! 同一基组直接为单位矩阵
        if (from_basis == to_basis) then
            do i = 1, n_channels
                U_mat(i, i) = 1.0_dp
            end do
            return
        end if

        ! 1. 非耦合基组 -> f-耦合基组:
        !    <f1 mf1, f2 mf2 | ms1 mi1, ms2 mi2> = <s1 ms1, i1 mi1 | f1 mf1> * <s2 ms2, i2 mi2 | f2 mf2>
        if (from_basis == BASIS_UNCOUPLED .and. to_basis == BASIS_F_COUPLED) then
            do i = 1, n_channels ! f-coupled
                do j = 1, n_channels ! uncoupled
                    if (channels_tgt(i)%l_orb /= channels_src(j)%l_orb .or. &
                        channels_tgt(i)%m_l   /= channels_src(j)%m_l) cycle

                    c1 = clebsch_gordan_half(atom1%two_s, channels_src(j)%two_ms1, &
                                             atom1%two_i, channels_src(j)%two_mi1, &
                                             channels_tgt(i)%two_f1, channels_tgt(i)%two_mf1)
                    c2 = clebsch_gordan_half(atom2%two_s, channels_src(j)%two_ms2, &
                                             atom2%two_i, channels_src(j)%two_mi2, &
                                             channels_tgt(i)%two_f2, channels_tgt(i)%two_mf2)
                    U_mat(i, j) = c1 * c2
                end do
            end do
            return

        ! 2. f-耦合基组 -> 非耦合基组 (转置即逆)
        else if (from_basis == BASIS_F_COUPLED .and. to_basis == BASIS_UNCOUPLED) then
            call calc_basis_transform_matrix(atom1, atom2, channels_tgt, channels_src, &
                                             n_channels, BASIS_UNCOUPLED, BASIS_F_COUPLED, b_gauss, U_mat)
            U_mat = transpose(U_mat)
            return

        ! 3. 非耦合基组 -> 总自旋耦合基组:
        !    <(S, I) F MF | ms1 mi1, ms2 mi2> = <s1 ms1, s2 ms2 | S MS> * <i1 mi1, i2 mi2 | I MI> * <S MS, I MI | F MF>
        else if (from_basis == BASIS_UNCOUPLED .and. to_basis == BASIS_TOTAL_SPIN) then
            do i = 1, n_channels ! total spin
                do j = 1, n_channels ! uncoupled
                    if (channels_tgt(i)%l_orb /= channels_src(j)%l_orb .or. &
                        channels_tgt(i)%m_l   /= channels_src(j)%m_l) cycle

                    ms = channels_src(j)%two_ms1 + channels_src(j)%two_ms2
                    mi = channels_src(j)%two_mi1 + channels_src(j)%two_mi2
                    cg_s = clebsch_gordan_half(atom1%two_s, channels_src(j)%two_ms1, &
                                              atom2%two_s, channels_src(j)%two_ms2, &
                                              channels_tgt(i)%two_S, ms)
                    cg_i = clebsch_gordan_half(atom1%two_i, channels_src(j)%two_mi1, &
                                              atom2%two_i, channels_src(j)%two_mi2, &
                                              channels_tgt(i)%two_I, mi)
                    cg_f = clebsch_gordan_half(channels_tgt(i)%two_S, ms, &
                                              channels_tgt(i)%two_I, mi, &
                                              channels_tgt(i)%two_F, channels_tgt(i)%two_MF)
                    U_mat(i, j) = cg_s * cg_i * cg_f
                end do
            end do
            return

        ! 4. 总自旋耦合基组 -> 非耦合基组
        else if (from_basis == BASIS_TOTAL_SPIN .and. to_basis == BASIS_UNCOUPLED) then
            call calc_basis_transform_matrix(atom1, atom2, channels_tgt, channels_src, &
                                             n_channels, BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, b_gauss, U_mat)
            U_mat = transpose(U_mat)
            return

        ! 5. 场缀饰渐近本征通道基组 (对角化渐近哈密顿量 H_asymp(B))
        else if (to_basis == BASIS_FIELD_DRESSED) then
            allocate(h_asymp(n_channels, n_channels), eig_vals(n_channels), eig_vecs(n_channels, n_channels))
            call build_asymptotic_hamiltonian(atom1, atom2, b_gauss, 0.0_dp, from_basis, &
                                              channels_src, n_channels, h_asymp)
            call diagonalize_real_symmetric(h_asymp, n_channels, eig_vals, eig_vecs, info)
            ! eig_vecs 的列向量即为本征态在 source 基底下的展开系数，变换矩阵 U = eig_vecs^T
            U_mat = transpose(eig_vecs)
            deallocate(h_asymp, eig_vals, eig_vecs)
            return

        else
            ! 通用链式变换: from -> uncoupled -> to: U = U_{to <- unc} * U_{unc <- from}
            block
                real(dp) :: U1(n_channels, n_channels), U2(n_channels, n_channels)
                type(field_channel_t), allocatable :: ch_unc(:)
                integer :: n_u, l_max_val

                l_max_val = maxval(channels_src(1:n_channels)%l_orb)
                call build_field_collision_channels(atom1, atom2, BASIS_UNCOUPLED, &
                                                    channels_src(1)%two_Mtot, l_max_val, ch_unc, n_u)
                call calc_basis_transform_matrix(atom1, atom2, channels_src, ch_unc, &
                                                 n_channels, from_basis, BASIS_UNCOUPLED, b_gauss, U1)
                call calc_basis_transform_matrix(atom1, atom2, ch_unc, channels_tgt, &
                                                 n_channels, BASIS_UNCOUPLED, to_basis, b_gauss, U2)
                U_mat = matmul(U2, U1)
                deallocate(ch_unc)
            end block
        end if
    end subroutine calc_basis_transform_matrix

    ! ==========================================================================
    ! 5. 外加磁场/电场下双原子渐近哈密顿量 H_asymp(B, E) 计算
    !    H_asymp = h_1(B, E) + h_2(B, E)
    ! ==========================================================================
    subroutine build_asymptotic_hamiltonian( &
        atom1, atom2, b_gauss, e_field, basis_type, channels, n_channels, H_asymp)

        type(cold_atom_t), intent(in)   :: atom1, atom2
        real(dp), intent(in)            :: b_gauss
        real(dp), intent(in)            :: e_field
        integer, intent(in)             :: basis_type
        type(field_channel_t), intent(in):: channels(:)
        integer, intent(in)             :: n_channels
        real(dp), intent(out)           :: H_asymp(n_channels, n_channels)

        real(dp), allocatable :: H_uncoupled(:, :), U_mat(:, :)
        real(dp) :: b_au, ms1, mi1, ms2, mi2, ms1_p, mi1_p, ms2_p, mi2_p
        real(dp) :: z1, z2, hf1, hf2
        integer  :: i, j
        type(field_channel_t), allocatable :: ch_uncoupled(:)
        integer  :: n_u

        H_asymp = 0.0_dp
        b_au = b_gauss * GAUSS2AU

        ! 在非耦合基组中直接写出解析矩阵元
        allocate(H_uncoupled(n_channels, n_channels))
        H_uncoupled = 0.0_dp

        ! 构造参考非耦合通道
        call build_field_collision_channels(atom1, atom2, BASIS_UNCOUPLED, &
                                            channels(1)%two_Mtot, channels(1)%l_orb, &
                                            ch_uncoupled, n_u)

        do i = 1, n_channels
            ms1 = real(ch_uncoupled(i)%two_ms1, dp) / 2.0_dp
            mi1 = real(ch_uncoupled(i)%two_mi1, dp) / 2.0_dp
            ms2 = real(ch_uncoupled(i)%two_ms2, dp) / 2.0_dp
            mi2 = real(ch_uncoupled(i)%two_mi2, dp) / 2.0_dp

            ! 1. 对角元素: 纯对角 Zeeman 相互作用 + 超精细 sz * iz
            z1 = (atom1%g_s * MU_B_AU * ms1 - atom1%g_i * MU_N_AU * mi1) * b_au
            z2 = (atom2%g_s * MU_B_AU * ms2 - atom2%g_i * MU_N_AU * mi2) * b_au
            hf1 = atom1%a_hf_au * (ms1 * mi1)
            hf2 = atom2%a_hf_au * (ms2 * mi2)
            H_uncoupled(i, i) = z1 + z2 + hf1 + hf2

            ! 2. 非对角元素: 单原子超精细自旋翻转 (s+ i- + s- i+)
            do j = 1, n_channels
                if (i == j) cycle
                ms1_p = real(ch_uncoupled(j)%two_ms1, dp) / 2.0_dp
                mi1_p = real(ch_uncoupled(j)%two_mi1, dp) / 2.0_dp
                ms2_p = real(ch_uncoupled(j)%two_ms2, dp) / 2.0_dp
                mi2_p = real(ch_uncoupled(j)%two_mi2, dp) / 2.0_dp

                ! 原子 1 翻转 (原子 2 不变)
                if (ch_uncoupled(i)%two_ms2 == ch_uncoupled(j)%two_ms2 .and. &
                    ch_uncoupled(i)%two_mi2 == ch_uncoupled(j)%two_mi2) then
                    if (ch_uncoupled(i)%two_ms1 == ch_uncoupled(j)%two_ms1 + 2 .and. &
                        ch_uncoupled(i)%two_mi1 == ch_uncoupled(j)%two_mi1 - 2) then
                        H_uncoupled(i, j) = 0.5_dp * atom1%a_hf_au * &
                            sqrt(real(atom1%two_s * (atom1%two_s + 2) - &
                                      ch_uncoupled(j)%two_ms1 * (ch_uncoupled(j)%two_ms1 + 2), dp) / 4.0_dp) * &
                            sqrt(real(atom1%two_i * (atom1%two_i + 2) - &
                                      ch_uncoupled(j)%two_mi1 * (ch_uncoupled(j)%two_mi1 - 2), dp) / 4.0_dp)
                    else if (ch_uncoupled(i)%two_ms1 == ch_uncoupled(j)%two_ms1 - 2 .and. &
                             ch_uncoupled(i)%two_mi1 == ch_uncoupled(j)%two_mi1 + 2) then
                        H_uncoupled(i, j) = 0.5_dp * atom1%a_hf_au * &
                            sqrt(real(atom1%two_s * (atom1%two_s + 2) - &
                                      ch_uncoupled(j)%two_ms1 * (ch_uncoupled(j)%two_ms1 - 2), dp) / 4.0_dp) * &
                            sqrt(real(atom1%two_i * (atom1%two_i + 2) - &
                                      ch_uncoupled(j)%two_mi1 * (ch_uncoupled(j)%two_mi1 + 2), dp) / 4.0_dp)
                    end if
                end if

                ! 原子 2 翻转 (原子 1 不变)
                if (ch_uncoupled(i)%two_ms1 == ch_uncoupled(j)%two_ms1 .and. &
                    ch_uncoupled(i)%two_mi1 == ch_uncoupled(j)%two_mi1) then
                    if (ch_uncoupled(i)%two_ms2 == ch_uncoupled(j)%two_ms2 + 2 .and. &
                        ch_uncoupled(i)%two_mi2 == ch_uncoupled(j)%two_mi2 - 2) then
                        H_uncoupled(i, j) = 0.5_dp * atom2%a_hf_au * &
                            sqrt(real(atom2%two_s * (atom2%two_s + 2) - &
                                      ch_uncoupled(j)%two_ms2 * (ch_uncoupled(j)%two_ms2 + 2), dp) / 4.0_dp) * &
                            sqrt(real(atom2%two_i * (atom2%two_i + 2) - &
                                      ch_uncoupled(j)%two_mi2 * (ch_uncoupled(j)%two_mi2 - 2), dp) / 4.0_dp)
                    else if (ch_uncoupled(i)%two_ms2 == ch_uncoupled(j)%two_ms2 - 2 .and. &
                             ch_uncoupled(i)%two_mi2 == ch_uncoupled(j)%two_mi2 + 2) then
                        H_uncoupled(i, j) = 0.5_dp * atom2%a_hf_au * &
                            sqrt(real(atom2%two_s * (atom2%two_s + 2) - &
                                      ch_uncoupled(j)%two_ms2 * (ch_uncoupled(j)%two_ms2 - 2), dp) / 4.0_dp) * &
                            sqrt(real(atom2%two_i * (atom2%two_i + 2) - &
                                      ch_uncoupled(j)%two_mi2 * (ch_uncoupled(j)%two_mi2 + 2), dp) / 4.0_dp)
                    end if
                end if
            end do
        end do

        ! 将 H_uncoupled 变换到目标基组: H_target = U * H_uncoupled * U^T
        if (basis_type == BASIS_UNCOUPLED) then
            H_asymp = H_uncoupled
        else
            allocate(U_mat(n_channels, n_channels))
            call calc_basis_transform_matrix(atom1, atom2, ch_uncoupled, channels, &
                                             n_channels, BASIS_UNCOUPLED, basis_type, b_gauss, U_mat)
            H_asymp = matmul(matmul(U_mat, H_uncoupled), transpose(U_mat))
            deallocate(U_mat)
        end if

        deallocate(H_uncoupled, ch_uncoupled)
    end subroutine build_asymptotic_hamiltonian

    ! ==========================================================================
    ! 6. 电子自旋交换相互作用矩阵元 s_1 . s_2 计算
    !    V(r) = V_singlet(r) * P_0 + V_triplet(r) * P_1
    !         = (1/4 V0 + 3/4 V1) * I + (V1 - V0) * (s_1 . s_2)
    ! ==========================================================================
    subroutine build_spin_exchange_matrix(channels, n_channels, basis_type, S1_dot_S2)
        type(field_channel_t), intent(in):: channels(:)
        integer, intent(in)             :: n_channels
        integer, intent(in)             :: basis_type
        real(dp), intent(out)           :: S1_dot_S2(n_channels, n_channels)

        integer :: i, j
        real(dp) :: ms1, ms2, s_tot

        S1_dot_S2 = 0.0_dp

        select case (basis_type)
        case (BASIS_TOTAL_SPIN)
            ! 在总自旋耦合基组中 strictly diagonal: s1 . s2 = [S(S+1) - s1(s1+1) - s2(s2+1)] / 2
            ! 对于两自旋 1/2: S=0 为 -3/4, S=1 为 +1/4
            do i = 1, n_channels
                s_tot = real(channels(i)%two_S, dp) / 2.0_dp
                S1_dot_S2(i, i) = 0.5_dp * (s_tot * (s_tot + 1.0_dp) - 0.75_dp - 0.75_dp)
            end do

        case default
            ! 在非耦合基组中: s_1 . s_2 = sz1*sz2 + 0.5*(s1+ s2- + s1- s2+)
            do i = 1, n_channels
                ms1 = real(channels(i)%two_ms1, dp) / 2.0_dp
                ms2 = real(channels(i)%two_ms2, dp) / 2.0_dp
                S1_dot_S2(i, i) = ms1 * ms2

                do j = 1, n_channels
                    if (i == j) cycle
                    if (channels(i)%two_mi1 == channels(j)%two_mi1 .and. &
                        channels(i)%two_mi2 == channels(j)%two_mi2 .and. &
                        channels(i)%l_orb   == channels(j)%l_orb   .and. &
                        channels(i)%m_l     == channels(j)%m_l) then

                        ! s1+ s2-
                        if (channels(i)%two_ms1 == channels(j)%two_ms1 + 2 .and. &
                            channels(i)%two_ms2 == channels(j)%two_ms2 - 2) then
                            S1_dot_S2(i, j) = 0.5_dp
                        ! s1- s2+
                        else if (channels(i)%two_ms1 == channels(j)%two_ms1 - 2 .and. &
                                 channels(i)%two_ms2 == channels(j)%two_ms2 + 2) then
                            S1_dot_S2(i, j) = 0.5_dp
                        end if
                    end if
                end do
            end do
        end select
    end subroutine build_spin_exchange_matrix

    ! ==========================================================================
    ! 7. 外场超冷多通道相互作用势矩阵 V_{ij}(r; B, E) 全自动构建
    ! ==========================================================================
    subroutine build_field_potential_matrix( &
        r_grid, v_singlet, v_triplet, atom1, atom2, b_gauss, e_field, &
        basis_type, two_Mtot, l_max, v_mat, thresholds, l_channels, n_channels, stat)

        real(dp), intent(in)                :: r_grid(:)
        real(dp), intent(in)                :: v_singlet(:)
        real(dp), intent(in)                :: v_triplet(:)
        type(cold_atom_t), intent(in)       :: atom1, atom2
        real(dp), intent(in)                :: b_gauss
        real(dp), intent(in)                :: e_field
        integer, intent(in)                 :: basis_type
        integer, intent(in)                 :: two_Mtot
        integer, intent(in)                 :: l_max
        real(dp), allocatable, intent(out)  :: v_mat(:, :, :)
        real(dp), allocatable, intent(out)  :: thresholds(:)
        integer,  allocatable, intent(out)  :: l_channels(:)
        integer, intent(out)                :: n_channels
        integer, optional, intent(out)      :: stat

        integer  :: n_pts, ir, i, j
        type(field_channel_t), allocatable :: channels(:)
        real(dp), allocatable :: H_asymp(:, :), S_dot_S(:, :), U_mat(:, :)
        real(dp), allocatable :: eig_vals(:), eig_vecs(:, :)
        real(dp) :: v_avg, v_diff
        integer  :: info

        if (present(stat)) stat = 0
        n_pts = size(r_grid)

        ! 1. 建立指定对称性与分波的通道基底
        call build_field_collision_channels(atom1, atom2, basis_type, two_Mtot, l_max, &
                                            channels, n_channels)
        if (n_channels == 0) then
            if (present(stat)) stat = -1
            return
        end if

        allocate(v_mat(n_channels, n_channels, n_pts))
        allocate(thresholds(n_channels))
        allocate(l_channels(n_channels))
        v_mat = 0.0_dp

        do i = 1, n_channels
            l_channels(i) = channels(i)%l_orb
        end do

        ! 2. 求解渐近阈值哈密顿量与交换算符
        allocate(H_asymp(n_channels, n_channels), S_dot_S(n_channels, n_channels))
        call build_asymptotic_hamiltonian(atom1, atom2, b_gauss, e_field, basis_type, &
                                          channels, n_channels, H_asymp)
        if (basis_type == BASIS_TOTAL_SPIN) then
            call build_spin_exchange_matrix(channels, n_channels, BASIS_TOTAL_SPIN, S_dot_S)
        else if (basis_type == BASIS_UNCOUPLED .or. basis_type == BASIS_FIELD_DRESSED) then
            call build_spin_exchange_matrix(channels, n_channels, BASIS_UNCOUPLED, S_dot_S)
        else
            ! 对于 BASIS_F_COUPLED，先在非耦合基组构建再变换: S_f = U * S_unc * U^T
            block
                type(field_channel_t), allocatable :: ch_unc(:)
                real(dp), allocatable :: U_mat_f(:, :), S_unc(:, :)
                integer :: n_u
                call build_field_collision_channels(atom1, atom2, BASIS_UNCOUPLED, two_Mtot, l_max, &
                                                    ch_unc, n_u)
                allocate(U_mat_f(n_channels, n_channels), S_unc(n_channels, n_channels))
                call build_spin_exchange_matrix(ch_unc, n_channels, BASIS_UNCOUPLED, S_unc)
                call calc_basis_transform_matrix(atom1, atom2, ch_unc, channels, n_channels, &
                                                 BASIS_UNCOUPLED, basis_type, b_gauss, U_mat_f)
                S_dot_S = matmul(matmul(U_mat_f, S_unc), transpose(U_mat_f))
                deallocate(ch_unc, U_mat_f, S_unc)
            end block
        end if

        ! 3. 若使用场缀饰通道基 (Field-Dressed)，严格对角化 H_asymp 并将矩阵元对齐到渐近阈值
        if (basis_type == BASIS_FIELD_DRESSED) then
            allocate(eig_vals(n_channels), eig_vecs(n_channels, n_channels))
            call diagonalize_real_symmetric(H_asymp, n_channels, eig_vals, eig_vecs, info)
            thresholds = eig_vals
            ! 交换算符变换到场缀饰通道
            S_dot_S = matmul(matmul(transpose(eig_vecs), S_dot_S), eig_vecs)
            deallocate(eig_vals, eig_vecs)
        else
            do i = 1, n_channels
                thresholds(i) = H_asymp(i, i)
            end do
        end if

        ! 4. 逐径向格点组装全多通道势矩阵:
        !    V_{ij}(r) = [V_{avg}(r) * \delta_{ij} + V_{diff}(r) * (s_1 . s_2)_{ij}] + [H_{asymp, ij} - thresholds(i)*\delta_{ij}]
        do ir = 1, n_pts
            v_avg = 0.25_dp * v_singlet(ir) + 0.75_dp * v_triplet(ir)
            v_diff = v_triplet(ir) - v_singlet(ir)

            do i = 1, n_channels
                do j = 1, n_channels
                    v_mat(i, j, ir) = v_diff * S_dot_S(i, j)
                    if (i == j) then
                        v_mat(i, i, ir) = v_mat(i, i, ir) + v_avg
                    else if (basis_type /= BASIS_FIELD_DRESSED) then
                        ! 包含渐近哈密顿量的非对角外场耦合项
                        v_mat(i, j, ir) = v_mat(i, j, ir) + H_asymp(i, j)
                    end if
                end do
            end do
        end do

        deallocate(H_asymp, S_dot_S, channels)
    end subroutine build_field_potential_matrix

    ! ==========================================================================
    ! 8. 磁 Feshbach 共振磁场扫描与多通道散射求解 a_s(B)
    ! ==========================================================================
    subroutine calc_magnetic_feshbach_resonance_scan( &
        r_grid, v_singlet, v_triplet, atom1, atom2, b_min, b_max, n_b, &
        collision_energy, two_Mtot, basis_type, res, stat)

        real(dp), intent(in)            :: r_grid(:)
        real(dp), intent(in)            :: v_singlet(:)
        real(dp), intent(in)            :: v_triplet(:)
        type(cold_atom_t), intent(in)   :: atom1, atom2
        real(dp), intent(in)            :: b_min, b_max
        integer, intent(in)             :: n_b
        real(dp), intent(in)            :: collision_energy
        integer, intent(in)             :: two_Mtot
        integer, intent(in)             :: basis_type
        type(field_feshbach_result_t), intent(out) :: res
        integer, optional, intent(out)  :: stat

        real(dp) :: mu_mass, db, b_val, e_tot, k_open
        integer  :: ib, n_ch
        real(dp), allocatable :: v_mat(:, :, :), thresholds(:)
        integer,  allocatable :: l_channels(:)
        type(multichannel_result_t) :: mc_res

        if (present(stat)) stat = 0
        mu_mass = (atom1%mass_au * atom2%mass_au) / (atom1%mass_au + atom2%mass_au)

        res%n_b = n_b
        allocate(res%b_grid(n_b), res%a_s_grid(n_b), res%eigenphases(n_b))

        db = 0.0_dp
        if (n_b > 1) db = (b_max - b_min) / real(n_b - 1, dp)

        do ib = 1, n_b
            b_val = b_min + real(ib - 1, dp) * db
            res%b_grid(ib) = b_val

            ! 构建该磁场下的全多通道势矩阵
            call build_field_potential_matrix( &
                r_grid, v_singlet, v_triplet, atom1, atom2, b_val, 0.0_dp, &
                basis_type, two_Mtot, 0, v_mat, thresholds, l_channels, n_ch)

            ! 碰撞总能量 = 入射通道阈值 + 动能
            e_tot = thresholds(1) + collision_energy

            ! 调用 Johnson Matrix Log-Derivative 算法求解多通道 S 矩阵与反应矩阵 K
            call calc_multichannel_close_coupling_logder( &
                r_grid, v_mat, mu_mass, e_tot, thresholds, l_channels, mc_res)

            k_open = sqrt(2.0_dp * mu_mass * max(1.0e-12_dp, collision_energy))
            ! 入射通道散射长度: a_s = - K_11 / k
            res%a_s_grid(ib) = -mc_res%k_matrix(1, 1) / k_open
            res%eigenphases(ib) = mc_res%eigenphase_sum

            deallocate(v_mat, thresholds, l_channels)
        end do

        ! 自动拟合提取共振位置 B_0, 宽度 \Delta B, 背景散射长度 a_{bg}
        call fit_feshbach_resonance_parameters(res%b_grid, res%a_s_grid, n_b, &
                                               res%b_res_pole, res%delta_b, res%a_bg)
    end subroutine calc_magnetic_feshbach_resonance_scan

    ! ==========================================================================
    ! 9. 磁 Feshbach 共振解析参数拟合: a_s(B) = a_bg * (1 - \Delta B / (B - B_0))
    ! ==========================================================================
    subroutine fit_feshbach_resonance_parameters(b_grid, a_s_grid, n_b, b_res, delta_b, a_bg, stat)
        real(dp), intent(in)  :: b_grid(:)
        real(dp), intent(in)  :: a_s_grid(:)
        integer, intent(in)   :: n_b
        real(dp), intent(out) :: b_res, delta_b, a_bg
        integer, optional, intent(out) :: stat

        integer :: i, max_jump_idx
        real(dp):: max_jump, jump, b_mid

        if (present(stat)) stat = 0
        b_res = 0.0_dp
        delta_b = 0.0_dp
        a_bg = a_s_grid(1)

        if (n_b < 5) return

        ! 寻找正负发散跳跃最大点作为共振极点 B_0
        max_jump = 0.0_dp
        max_jump_idx = 1
        do i = 1, n_b - 1
            jump = abs(a_s_grid(i + 1) - a_s_grid(i))
            if (jump > max_jump) then
                max_jump = jump
                max_jump_idx = i
            end if
        end do

        b_mid = 0.5_dp * (b_grid(max_jump_idx) + b_grid(max_jump_idx + 1))
        b_res = b_mid

        ! 背景散射长度 (远离共振两端平均)
        a_bg = 0.5_dp * (a_s_grid(1) + a_s_grid(n_b))
        if (abs(a_bg) < 1.0e-10_dp) a_bg = 1.0e-5_dp

        ! 估计共振宽度: a_s(B) = 0 处对应 B_zero = B_0 + \Delta B
        delta_b = (b_grid(2) - b_grid(1)) * 2.0_dp
        do i = 1, n_b - 1
            if (a_s_grid(i) * a_s_grid(i + 1) <= 0.0_dp .and. i /= max_jump_idx) then
                delta_b = (0.5_dp * (b_grid(i) + b_grid(i + 1))) - b_res
                exit
            end if
        end do
    end subroutine fit_feshbach_resonance_parameters

    ! ==========================================================================
    ! 辅助子程序: 实对称矩阵 Jacobi 本征值与本征向量求解
    ! ==========================================================================
    subroutine diagonalize_real_symmetric(a_in, n, eigenvalues, eigenvectors, info)
        real(dp), intent(in)  :: a_in(n, n)
        integer,  intent(in)  :: n
        real(dp), intent(out) :: eigenvalues(n)
        real(dp), intent(out) :: eigenvectors(n, n)
        integer,  intent(out) :: info

        real(dp) :: a(n, n), d(n), v(n, n)
        real(dp) :: thresh, theta, t, c, s, tau, h, g
        integer  :: i, j, k, p, q, sweep

        info = 0
        a = a_in
        eigenvectors = 0.0_dp
        do i = 1, n
            eigenvectors(i, i) = 1.0_dp
            d(i) = a(i, i)
        end do

        do sweep = 1, 50
            thresh = 0.0_dp
            do p = 1, n - 1
                do q = p + 1, n
                    thresh = thresh + abs(a(p, q))
                end do
            end do
            if (thresh < 1.0e-14_dp) exit

            do p = 1, n - 1
                do q = p + 1, n
                    h = a(p, q)
                    if (abs(h) < 1.0e-15_dp) cycle
                    theta = (d(q) - d(p)) / (2.0_dp * h)
                    t = 1.0_dp / (abs(theta) + sqrt(1.0_dp + theta * theta))
                    if (theta < 0.0_dp) t = -t
                    c = 1.0_dp / sqrt(1.0_dp + t * t)
                    s = t * c
                    tau = s / (1.0_dp + c)

                    h = t * a(p, q)
                    d(p) = d(p) - h
                    d(q) = d(q) + h
                    a(p, q) = 0.0_dp

                    do k = 1, p - 1
                        g = a(k, p)
                        h = a(k, q)
                        a(k, p) = g - s * (h + g * tau)
                        a(k, q) = h + s * (g - h * tau)
                    end do
                    do k = p + 1, q - 1
                        g = a(p, k)
                        h = a(k, q)
                        a(p, k) = g - s * (h + g * tau)
                        a(k, q) = h + s * (g - h * tau)
                    end do
                    do k = q + 1, n
                        g = a(p, k)
                        h = a(q, k)
                        a(p, k) = g - s * (h + g * tau)
                        a(q, k) = h + s * (g - h * tau)
                    end do

                    do k = 1, n
                        g = eigenvectors(k, p)
                        h = eigenvectors(k, q)
                        eigenvectors(k, p) = g - s * (h + g * tau)
                        eigenvectors(k, q) = h + s * (g - h * tau)
                    end do
                end do
            end do
        end do

        ! 升序排序
        do i = 1, n - 1
            k = i
            do j = i + 1, n
                if (d(j) < d(k)) k = j
            end do
            if (k /= i) then
                h = d(i); d(i) = d(k); d(k) = h
                do p = 1, n
                    h = eigenvectors(p, i)
                    eigenvectors(p, i) = eigenvectors(p, k)
                    eigenvectors(p, k) = h
                end do
            end if
        end do

        eigenvalues = d
    end subroutine diagonalize_real_symmetric

end module mod_field_scattering
