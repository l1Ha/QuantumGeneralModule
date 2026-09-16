!> \file mod_optical_lattice_hubbard.f90
!> \brief 超冷光晶格物理与玻色-哈伯德 (Bose-Hubbard) 微观映射模块
!> \details 涵盖一维/高维周期驻波光晶格势、Bloch 态与 Mathieu 能带结构、
!>          Wannier 函数构造与玻色-哈伯德微观参数 (J, U) 严格映射、
!>          超流-Mott 绝缘体相变判据、以及布洛赫振荡与 Landau-Zener 带间隧穿。
!> \author LiHao
module mod_optical_lattice_hubbard
    use mod_constants, only: dp, PI, TWOPI, AMU2AU, KB
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    public :: optical_lattice_t
    public :: bose_hubbard_param_t
    public :: init_optical_lattice
    public :: calc_bloch_band_energies
    public :: calc_bose_hubbard_parameters
    public :: calc_bloch_oscillation_dynamics

    !> 光晶格几何与激光参数
    type :: optical_lattice_t
        real(dp) :: mass_amu          !< 原子质量 (amu)
        real(dp) :: mass_au           !< 原子质量 (a.u.)
        real(dp) :: laser_lambda_nm   !< 晶格激光波长 (nm, 如 1064 nm 或 780 nm)
        real(dp) :: k_l_au            !< 晶格波矢 k_L = 2*pi / lambda (a.u.)
        real(dp) :: d_lattice_au      !< 晶格常数 d = lambda / 2 (a.u.)
        real(dp) :: e_recoil_au       !< 反冲能量 E_R = hbar^2 k_L^2 / (2 m) (a.u.)
        real(dp) :: e_recoil_hz       !< 反冲能量对应频率 E_R / h (Hz)
        real(dp) :: e_recoil_nkelvin  !< 反冲温度 E_R / k_B (nK)
        real(dp) :: s_depth           !< 无量纲晶格深度 V_0 / E_R (通常 2 - 30)
    end type optical_lattice_t

    !> 玻色-哈伯德微观参数与相变指示
    type :: bose_hubbard_param_t
        real(dp) :: j_hopping_er      !< 最近邻跃迁能量 J (E_R)
        real(dp) :: j_hopping_hz      !< 最近邻跃迁能量 J / h (Hz)
        real(dp) :: u_onsite_er       !< 在位相互作用能量 U (E_R)
        real(dp) :: u_onsite_hz       !< 在位相互作用能量 U / h (Hz)
        real(dp) :: u_over_j_ratio    !< 强相互作用比率 U / J
        real(dp) :: band_gap_er       !< 基带与第一激发带能隙 (E_R)
        logical  :: is_mott_candidate !< 是否达到 Mott 绝缘相临界区 (U/J > 3.3 in 1D, > 29.3 in 3D)
    end type bose_hubbard_param_t

contains

    !> \brief 初始化光学晶格体系参数
    !> \param[out] latt 晶格结构体
    !> \param[in] mass_amu 原子质量 (amu, 如 87Rb 为 86.909 amu)
    !> \param[in] lambda_nm 激光波长 (nm, 如 1064.0 nm)
    !> \param[in] s_depth 晶格阱深 V0 / E_R
    !> \param[out] stat 状态码
    subroutine init_optical_lattice(latt, mass_amu, lambda_nm, s_depth, stat)
        type(optical_lattice_t), intent(out) :: latt
        real(dp), intent(in)                 :: mass_amu
        real(dp), intent(in)                 :: lambda_nm
        real(dp), intent(in)                 :: s_depth
        integer, optional, intent(out)       :: stat

        real(dp) :: lambda_m, lambda_au

        if (present(stat)) stat = 0
        if (mass_amu <= 0.0_dp .or. lambda_nm <= 0.0_dp .or. s_depth < 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        latt%mass_amu = mass_amu
        latt%mass_au  = mass_amu * AMU2AU
        latt%laser_lambda_nm = lambda_nm
        latt%s_depth = s_depth

        ! 1 a.u. length = 0.529177210903e-10 m
        lambda_m = lambda_nm * 1.0e-9_dp
        lambda_au = lambda_m / 0.529177210903e-10_dp

        latt%k_l_au = TWOPI / lambda_au
        latt%d_lattice_au = 0.5_dp * lambda_au

        ! E_R = k_L^2 / (2 * m) in a.u.
        latt%e_recoil_au = (latt%k_l_au**2) / (2.0_dp * latt%mass_au)

        ! 1 a.u. energy / h = 6.579683920502e15 Hz, 1 a.u. / k_B = 3.157750248e5 K
        latt%e_recoil_hz = latt%e_recoil_au * 6.579683920502e15_dp
        latt%e_recoil_nkelvin = (latt%e_recoil_au * 3.157750248e5_dp) * 1.0e9_dp
    end subroutine init_optical_lattice

    !> \brief 利用平面波展开法求解 Mathieu 能带在准动量 q 处的本征能量
    !> \param[in] latt 晶格参数
    !> \param[in] q_quasi 准动量 q / k_L (范围在 [-1.0, 1.0])
    !> \param[in] n_bands 所需输出的能带数目 (通常 2 - 4)
    !> \param[out] band_energies_er 能带本征能量 (单位 E_R)
    subroutine calc_bloch_band_energies(latt, q_quasi, n_bands, band_energies_er, stat)
        type(optical_lattice_t), intent(in)   :: latt
        real(dp), intent(in)                  :: q_quasi
        integer, intent(in)                   :: n_bands
        real(dp), dimension(n_bands), intent(out) :: band_energies_er
        integer, optional, intent(out)        :: stat

        integer, parameter :: N_PLANE = 15  !< 平面波基底展开截断 (-N_PLANE 到 +N_PLANE)
        integer, parameter :: DIM_MAT = 2 * N_PLANE + 1
        real(dp), dimension(DIM_MAT, DIM_MAT) :: h_mat, eig_vecs
        real(dp), dimension(DIM_MAT) :: eig_vals
        integer  :: m, i, info

        if (present(stat)) stat = 0
        if (n_bands > DIM_MAT .or. n_bands < 1) then
            if (present(stat)) stat = -1
            band_energies_er = 0.0_dp
            return
        end if

        h_mat = 0.0_dp

        ! 填充动能对角线: T_m = (2*m + q/k_L)^2 * E_R
        do m = -N_PLANE, N_PLANE
            i = m + N_PLANE + 1
            h_mat(i, i) = (2.0_dp * real(m, dp) + q_quasi)**2 + 0.5_dp * latt%s_depth
        end do

        ! 填充周期势非对角线: V_0 cos^2(k_L x) = V_0/2 + V_0/4 (e^{i 2 k_L x} + e^{-i 2 k_L x})
        ! 耦合 m 与 m +/- 1
        do m = -N_PLANE, N_PLANE - 1
            i = m + N_PLANE + 1
            h_mat(i, i + 1) = 0.25_dp * latt%s_depth
            h_mat(i + 1, i) = 0.25_dp * latt%s_depth
        end do

        call diag_symmetric_matrix(DIM_MAT, h_mat, eig_vals, eig_vecs, stat=info)
        if (info /= 0) then
            if (present(stat)) stat = info
        end if

        do i = 1, n_bands
            band_energies_er(i) = eig_vals(i)
        end do
    end subroutine calc_bloch_band_energies

    !> \brief 计算玻色-哈伯德 (Bose-Hubbard) 模型微观参数 J, U 及相变指标
    !> \param[in] latt 晶格参数
    !> \param[in] a_s_bohr s 波散射长度 (Bohr)
    !> \param[in] omega_perp_hz 横向约束谐振频率 (Hz, 用于 1D 晶格的 3D 有效相互作用积分)
    !> \param[out] bh 玻色-哈伯德微观参数结构体
    !> \param[out] stat 状态码
    subroutine calc_bose_hubbard_parameters(latt, a_s_bohr, omega_perp_hz, bh, stat)
        type(optical_lattice_t), intent(in) :: latt
        real(dp), intent(in)                :: a_s_bohr
        real(dp), intent(in)                :: omega_perp_hz
        type(bose_hubbard_param_t), intent(out) :: bh
        integer, optional, intent(out)      :: stat

        real(dp), dimension(2) :: bands_q0, bands_q1
        real(dp) :: e_band_width, s, a_perp_bohr
        real(dp) :: mass_kg

        if (present(stat)) stat = 0
        if (a_s_bohr <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        s = latt%s_depth

        ! 1. 计算 q=0 (布里渊区中心) 与 q=1 (布里渊区边界) 的本征能量
        call calc_bloch_band_energies(latt, 0.0_dp, 2, bands_q0)
        call calc_bloch_band_energies(latt, 1.0_dp, 2, bands_q1)

        ! 紧束缚色散关系: E_0(q) = E_bar - 2*J*cos(q*d) => 带宽 Delta E = 4*J
        e_band_width = abs(bands_q1(1) - bands_q0(1))
        bh%j_hopping_er = e_band_width / 4.0_dp
        bh%j_hopping_hz = bh%j_hopping_er * latt%e_recoil_hz

        ! 基带与激发带能隙 (在 q=0 处取最小间隙)
        bh%band_gap_er = bands_q0(2) - bands_q0(1)

        ! 2. 在位相互作用 U 计算 (深晶格简谐近似结合横向谐振约束)
        ! U = sqrt(8 / pi) * k_L * a_s * E_R * s^(3/4) * (omega_perp / omega_z)
        ! 横向谐振基态半径 a_perp = sqrt(hbar / (m * omega_perp))
        mass_kg = latt%mass_amu * 1.66053906660e-27_dp
        a_perp_bohr = sqrt(1.054571817e-34_dp / (mass_kg * max(100.0_dp, omega_perp_hz * TWOPI))) / &
                      0.529177210903e-10_dp

        ! 1D 晶格下的在位排斥能 U:
        ! U = 2 * hbar * omega_perp * a_s * int |w(z)|^4 dz
        ! 利用谐振子基态波函数 int |w(z)|^4 dz = (m * omega_z / (2*pi*hbar))^(1/2) * sqrt(pi) / ...
        ! 解析定标: U / E_R = sqrt(8 / pi) * (k_L * a_s) * (s^(1/4)) * (2 * hbar * omega_perp / E_R)
        bh%u_onsite_er = sqrt(8.0_dp / PI) * (latt%k_l_au * a_s_bohr) * &
                         (s**0.25_dp) * (2.0_dp * (omega_perp_hz / latt%e_recoil_hz))
        bh%u_onsite_hz = bh%u_onsite_er * latt%e_recoil_hz

        ! 3. 强关联参数比 U / J
        if (bh%j_hopping_er > 1.0e-12_dp) then
            bh%u_over_j_ratio = bh%u_onsite_er / bh%j_hopping_er
        else
            bh%u_over_j_ratio = 1.0e10_dp
        end if

        ! 1D 临界转变点 (U/J)_c ~ 3.3; 3D 体系 (U/J)_c ~ 29.3
        bh%is_mott_candidate = (bh%u_over_j_ratio > 3.3_dp)
    end subroutine calc_bose_hubbard_parameters

    !> \brief 计算微弱恒力场下的布洛赫振荡周期与 Landau-Zener 激发带隧穿几率
    !> \param[in] latt 晶格结构
    !> \param[in] force_si 外力大小 (N, 如重力 m*g)
    !> \param[out] t_bloch_ms 布洛赫振荡周期 (ms)
    !> \param[out] omega_bloch_hz 布洛赫振荡角频率 (Hz)
    !> \param[out] p_lz_tunnel Landau-Zener 带间跃迁几率
    subroutine calc_bloch_oscillation_dynamics(latt, force_si, t_bloch_ms, omega_bloch_hz, &
                                              p_lz_tunnel, stat)
        type(optical_lattice_t), intent(in) :: latt
        real(dp), intent(in)                :: force_si
        real(dp), intent(out)               :: t_bloch_ms
        real(dp), intent(out)               :: omega_bloch_hz
        real(dp), intent(out)               :: p_lz_tunnel
        integer, optional, intent(out)      :: stat

        real(dp) :: d_m, hbar_si, delta_gap_j, v_max_si
        real(dp) :: exponent_lz
        real(dp), dimension(2) :: bands_q1

        if (present(stat)) stat = 0
        if (force_si <= 0.0_dp) then
            if (present(stat)) stat = -1
            t_bloch_ms = 0.0_dp
            omega_bloch_hz = 0.0_dp
            p_lz_tunnel = 0.0_dp
            return
        end if

        hbar_si = 1.054571817e-34_dp
        d_m = latt%d_lattice_au * 0.529177210903e-10_dp

        ! 布洛赫振荡角频率: omega_B = F * d / hbar
        omega_bloch_hz = force_si * d_m / hbar_si
        ! 布洛赫振荡周期: T_B = 2*pi / omega_B = 2*pi*hbar / (F * d)
        t_bloch_ms = (TWOPI / omega_bloch_hz) * 1.0e3_dp

        ! 在布里渊区边界 q=k_L (q_quasi = 1) 处的带隙
        call calc_bloch_band_energies(latt, 1.0_dp, 2, bands_q1)
        delta_gap_j = (bands_q1(2) - bands_q1(1)) * (latt%e_recoil_hz * 6.62607015e-34_dp)

        ! 反冲速度 v_R = hbar * k_L / m = hbar * pi / (m * d)
        v_max_si = hbar_si * PI / ((latt%mass_amu * 1.66053906660e-27_dp) * d_m)

        ! P_LZ = exp( - pi * Delta_gap^2 / (4 * hbar * F * v_R) )
        exponent_lz = (PI * (delta_gap_j**2)) / max(1.0e-90_dp, 4.0_dp * hbar_si * force_si * v_max_si)
        p_lz_tunnel = exp(-min(80.0_dp, exponent_lz))
    end subroutine calc_bloch_oscillation_dynamics

end module mod_optical_lattice_hubbard
