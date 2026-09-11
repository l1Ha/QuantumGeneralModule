!> \brief 强场原子物理模型势与强场参数模块
!> \details 提供单活性电子（SAE）软核库仑势及其一阶导数、常见稀有气体原子模型参数、
!>          Keldysh 参数、有质动力能、高次谐波截断能量定律以及 ADK 隧穿电离率计算。
!> \author LiHao
module mod_coulomb_atomic
    use mod_constants, only: dp, PI, TWOPI, EV2AU, AU2EV
    implicit none
    private

    public :: atom_config_t
    public :: get_atom_config
    public :: soft_core_coulomb_potential
    public :: soft_core_coulomb_derivative
    public :: keldysh_parameter
    public :: ponderomotive_energy
    public :: hhg_cutoff_energy
    public :: quiver_radius
    public :: adk_ionization_rate

    !> \brief 原子/分子单活性电子物理参数结构体
    type :: atom_config_t
        character(len=16) :: name = "H"
        real(dp) :: ip_au = 0.5_dp           !< 电离能 (Hartree)
        real(dp) :: soft_core_a = 1.0_dp     !< 1D 软核屏蔽参数 (a.u.)
        real(dp) :: z_eff = 1.0_dp           !< 有效残余核电荷
        integer  :: l_quantum = 0            !< 轨道角动量 l
        integer  :: m_quantum = 0            !< 磁量子数 m
    end type atom_config_t

contains

    !> \brief 获取常见原子的标准物理参数配置
    !> \param[in] atom_name 原子标识，如 'H', 'He', 'Ne', 'Ar', 'Kr', 'Xe'
    !> \param[out] cfg 原子参数结构体
    subroutine get_atom_config(atom_name, cfg)
        character(len=*), intent(in) :: atom_name
        type(atom_config_t), intent(out) :: cfg
        character(len=16) :: name_clean

        name_clean = trim(adjustl(atom_name))
        cfg%name = name_clean
        cfg%z_eff = 1.0_dp
        cfg%l_quantum = 0
        cfg%m_quantum = 0

        select case(name_clean)
        case('H', 'h', 'Hydrogen')
            cfg%ip_au = 0.5_dp                   ! 13.606 eV
            cfg%soft_core_a = 1.0_dp             ! 1D 基态 E0 = -0.5 a.u.
        case('He', 'he', 'Helium')
            cfg%ip_au = 0.90357_dp               ! 24.587 eV
            cfg%soft_core_a = 0.751_dp
        case('Ne', 'ne', 'Neon')
            cfg%ip_au = 0.79248_dp               ! 21.565 eV
            cfg%soft_core_a = 0.885_dp
            cfg%l_quantum = 1                    ! 2p
        case('Ar', 'ar', 'Argon')
            cfg%ip_au = 0.57915_dp               ! 15.760 eV
            cfg%soft_core_a = 1.340_dp
            cfg%l_quantum = 1                    ! 3p
        case('Kr', 'kr', 'Krypton')
            cfg%ip_au = 0.51446_dp               ! 13.999 eV
            cfg%soft_core_a = 1.480_dp
            cfg%l_quantum = 1                    ! 4p
        case('Xe', 'xe', 'Xenon')
            cfg%ip_au = 0.44577_dp               ! 12.130 eV
            cfg%soft_core_a = 1.620_dp
            cfg%l_quantum = 1                    ! 5p
        case default
            cfg%ip_au = 0.5_dp
            cfg%soft_core_a = 1.0_dp
        end select
    end subroutine get_atom_config

    !> \brief 一维软核库仑势: V(x) = -Z_eff / sqrt(x^2 + a^2)
    pure function soft_core_coulomb_potential(x, soft_a, z_eff) result(v)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: soft_a
        real(dp), intent(in), optional :: z_eff
        real(dp) :: v, z_val

        z_val = 1.0_dp
        if (present(z_eff)) z_val = z_eff
        v = -z_val / sqrt(x**2 + soft_a**2)
    end function soft_core_coulomb_potential

    !> \brief 软核库仑势一阶导数: dV/dx = Z_eff * x / (x^2 + a^2)^(3/2)
    !> \details 用于 Ehrenfest 偶极加速度中受力计算: F(x) = -dV/dx
    pure function soft_core_coulomb_derivative(x, soft_a, z_eff) result(dvdx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: soft_a
        real(dp), intent(in), optional :: z_eff
        real(dp) :: dvdx, z_val

        z_val = 1.0_dp
        if (present(z_eff)) z_val = z_eff
        dvdx = z_val * x / ((x**2 + soft_a**2)**(1.5_dp))
    end function soft_core_coulomb_derivative

    !> \brief 计算 Keldysh 参数: gamma = omega * sqrt(2 * Ip) / E_0
    pure function keldysh_parameter(omega, e_peak, ip_au) result(gamma)
        real(dp), intent(in) :: omega, e_peak, ip_au
        real(dp) :: gamma

        if (e_peak <= 1.0e-15_dp) then
            gamma = 1.0e10_dp
        else
            gamma = omega * sqrt(2.0_dp * ip_au) / e_peak
        end if
    end function keldysh_parameter

    !> \brief 计算激光场中有质动力能 (Ponderomotive Energy): U_p = E_0^2 / (4 * omega^2)
    pure function ponderomotive_energy(e_peak, omega) result(up)
        real(dp), intent(in) :: e_peak, omega
        real(dp) :: up

        if (omega <= 1.0e-15_dp) then
            up = 0.0_dp
        else
            up = (e_peak**2) / (4.0_dp * omega**2)
        end if
    end function ponderomotive_energy

    !> \brief 计算高次谐波半经典截断截止能量: E_cutoff = I_p + 3.17 * U_p
    pure function hhg_cutoff_energy(ip_au, e_peak, omega) result(e_cut)
        real(dp), intent(in) :: ip_au, e_peak, omega
        real(dp) :: e_cut

        e_cut = ip_au + 3.1725955_dp * ponderomotive_energy(e_peak, omega)
    end function hhg_cutoff_energy

    !> \brief 计算自由电子在激光场中的晃动振幅 (Quiver radius): alpha_0 = E_0 / omega^2
    pure function quiver_radius(e_peak, omega) result(alpha_0)
        real(dp), intent(in) :: e_peak, omega
        real(dp) :: alpha_0

        if (omega <= 1.0e-15_dp) then
            alpha_0 = 0.0_dp
        else
            alpha_0 = e_peak / (omega**2)
        end if
    end function quiver_radius

    !> \brief ADK (Ammosov-Delone-Krainov) 静态/准静态隧穿电离率
    !> \param[in] e_field 瞬时电场绝对值 |E(t)| (a.u.)
    !> \param[in] ip_au 电离势 (a.u.)
    !> \param[in] z_eff 有效电荷
    !> \param[in] l 轨道角动量
    !> \param[in] m 磁量子数
    pure function adk_ionization_rate(e_field, ip_au, z_eff, l, m) result(w_adk)
        real(dp), intent(in) :: e_field, ip_au
        real(dp), intent(in), optional :: z_eff
        integer, intent(in), optional :: l, m
        real(dp) :: w_adk

        real(dp) :: z, ef, n_star, c_nl_sq, factor_f, e_crit
        integer :: l_val, m_val

        z = 1.0_dp
        l_val = 0
        m_val = 0
        if (present(z_eff)) z = z_eff
        if (present(l)) l_val = l
        if (present(m)) m_val = abs(m)

        ef = abs(e_field)
        if (ef <= 1.0e-6_dp .or. ip_au <= 1.0e-6_dp) then
            w_adk = 0.0_dp
            return
        end if

        n_star = z / sqrt(2.0_dp * ip_au)

        ! C_nl^2 系数经验公式
        c_nl_sq = (2.0_dp**(2.0_dp * n_star)) / (n_star * gamma_approx(n_star + 1.0_dp) * gamma_approx(n_star))
        factor_f = real(2 * l_val + 1, dp) * real(factorial_int(l_val + m_val), dp) / &
                   real(2**m_val * factorial_int(m_val) * factorial_int(l_val - m_val), dp)

        e_crit = 2.0_dp * (2.0_dp * ip_au)**(1.5_dp)
        w_adk = c_nl_sq * factor_f * ip_au * &
                ((e_crit / ef)**(2.0_dp * n_star - real(m_val, dp) - 1.0_dp)) * &
                exp(-e_crit / (3.0_dp * ef))
    end function adk_ionization_rate

    !> \brief Lanczos 伽马函数近似 Gamma(x)
    pure function gamma_approx(x) result(g)
        real(dp), intent(in) :: x
        real(dp) :: g
        real(dp) :: t, s
        real(dp), parameter :: c(6) = [ &
            76.18009172947146_dp, -86.50532032941677_dp, &
            24.01409824083091_dp, -1.231739572450155_dp, &
            0.1208650973866179e-2_dp, -0.5395239384953e-5_dp ]
        integer :: j

        if (x <= 0.0_dp) then
            g = 1.0_dp
            return
        end if
        t = x + 4.5_dp
        s = 1.000000000190015_dp
        do j = 1, 6
            s = s + c(j) / (x + real(j - 1, dp))
        end do
        g = exp((x - 0.5_dp) * log(t) - t + 0.918938533204673_dp) * s
    end function gamma_approx

    !> \brief 整数阶乘辅助函数
    pure function factorial_int(n) result(res)
        integer, intent(in) :: n
        integer :: res, i
        res = 1
        if (n <= 1) return
        do i = 2, n
            res = res * i
        end do
    end function factorial_int

end module mod_coulomb_atomic
