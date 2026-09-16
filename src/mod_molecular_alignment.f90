!> \file mod_molecular_alignment.f90
!> \brief 强激光诱导分子定向、取向与超转子动力学模块
!> \details 涵盖刚体转子各向异性极化势相互作用、非绝热激光脉冲引起的波包演化、
!>          无场转动复苏序参量 <cos^2 theta>(t) 与双色场取向 <cos theta>(t)、
!>          以及光学离心机 (Optical Centrifuge) 加速至分子超转子 (Superrotor) 的极化动力学。
!> \author LiHao
module mod_molecular_alignment
    use mod_constants, only: dp, PI, TWOPI, EYE, KB
    implicit none
    private

    public :: rotor_molecule_t
    public :: init_rotor_molecule
    public :: calc_cos2_matrix_elements
    public :: calc_cos_matrix_elements
    public :: simulate_laser_induced_alignment
    public :: calc_optical_centrifuge_kick
    public :: calc_superrotor_dissociation

    ! 换算常数
    real(dp), parameter :: CM1_TO_AU = 4.5563352529120e-6_dp
    real(dp), parameter :: FS_TO_AU  = 41.341374575751_dp
    real(dp), parameter :: PS_TO_AU  = 41341.374575751_dp
    real(dp), parameter :: DEBYE_TO_AU = 0.393430307_dp
    real(dp), parameter :: WCM2_TO_AU = 1.0_dp / 3.5094452e16_dp

    !> 刚体线性分子转子参数
    type :: rotor_molecule_t
        character(len=16) :: name
        real(dp) :: b_rot_cm1         !< 转动常数 B_e (cm^-1)
        real(dp) :: b_rot_au          !< 转动常数 B_e (a.u.)
        real(dp) :: delta_alpha_au    !< 极化率各向异性 Delta alpha = alpha_par - alpha_perp (a.u.)
        real(dp) :: dipole_debye      !< 永久偶极矩 (Debye)
        real(dp) :: dipole_au         !< 永久偶极矩 (a.u.)
        real(dp) :: t_rev_ps          !< 经典转动全复苏周期 T_rev = 1/(2*B*c) (ps)
    end type rotor_molecule_t

contains

    !> \brief 初始化线性转子分子参数
    !> \param[out] mol 分子结构体
    !> \param[in] name 分子标识符 (如 "N2", "O2", "CO", "I2")
    !> \param[in] b_rot_cm1 转动常数 (cm^-1, 如 N2 ~ 1.998 cm^-1, I2 ~ 0.037 cm^-1)
    !> \param[in] delta_alpha_au 极化率各向异性 (a.u., 如 N2 ~ 6.7 a.u., I2 ~ 46 a.u.)
    !> \param[in] dipole_debye 永久偶极矩 (Debye, 如 N2 为 0, CO ~ 0.112 D)
    !> \param[out] stat 状态码
    subroutine init_rotor_molecule(mol, name, b_rot_cm1, delta_alpha_au, dipole_debye, stat)
        type(rotor_molecule_t), intent(out) :: mol
        character(len=*), intent(in)        :: name
        real(dp), intent(in)                :: b_rot_cm1
        real(dp), intent(in)                :: delta_alpha_au
        real(dp), intent(in)                :: dipole_debye
        integer, optional, intent(out)      :: stat

        if (present(stat)) stat = 0
        if (b_rot_cm1 <= 0.0_dp .or. delta_alpha_au <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        mol%name = trim(name)
        mol%b_rot_cm1 = b_rot_cm1
        mol%b_rot_au = b_rot_cm1 * CM1_TO_AU
        mol%delta_alpha_au = delta_alpha_au
        mol%dipole_debye = dipole_debye
        mol%dipole_au = dipole_debye * DEBYE_TO_AU

        ! T_rev = 1 / (2 * B * c): 1 cm^-1 对应光周期 33.3564 ps => T_rev = 33.3564 / (2 * B)
        mol%t_rev_ps = 33.3564095198152_dp / (2.0_dp * b_rot_cm1)
    end subroutine init_rotor_molecule

    !> \brief 计算 cos^2(theta) 在 |J, M> 球谐基下的对角与 off-diagonal (Delta J = 2) 矩阵元
    subroutine calc_cos2_matrix_elements(j_max, m_proj, cos2_diag, cos2_off2)
        integer, intent(in)                    :: j_max, m_proj
        real(dp), dimension(0:j_max), intent(out)   :: cos2_diag
        real(dp), dimension(0:j_max-2), intent(out) :: cos2_off2

        integer  :: j
        real(dp) :: j_d, m_d

        cos2_diag = 0.0_dp
        cos2_off2 = 0.0_dp
        m_d = real(abs(m_proj), dp)

        do j = abs(m_proj), j_max
            j_d = real(j, dp)
            ! <J, M | cos^2 theta | J, M> = 1/3 + 2/3 * (J(J+1) - 3*M^2) / ((2J-1)*(2J+3))
            if (j == 0) then
                cos2_diag(j) = 1.0_dp / 3.0_dp
            else
                cos2_diag(j) = 1.0_dp / 3.0_dp + &
                    (2.0_dp / 3.0_dp) * (j_d * (j_d + 1.0_dp) - 3.0_dp * (m_d**2)) / &
                    max(1.0e-6_dp, (2.0_dp * j_d - 1.0_dp) * (2.0_dp * j_d + 3.0_dp))
            end if
        end do

        do j = abs(m_proj), j_max - 2
            j_d = real(j, dp)
            ! <J+2, M | cos^2 theta | J, M>
            cos2_off2(j) = sqrt(max(0.0_dp, &
                (j_d - m_d + 1.0_dp) * (j_d - m_d + 2.0_dp) * &
                (j_d + m_d + 1.0_dp) * (j_d + m_d + 2.0_dp)) / &
                max(1.0e-6_dp, (2.0_dp * j_d + 1.0_dp) * ((2.0_dp * j_d + 3.0_dp)**2) * (2.0_dp * j_d + 5.0_dp)))
        end do
    end subroutine calc_cos2_matrix_elements

    !> \brief 计算 cos(theta) 在 |J, M> 球谐基下的 off-diagonal (Delta J = 1) 矩阵元
    subroutine calc_cos_matrix_elements(j_max, m_proj, cos_off1)
        integer, intent(in)                    :: j_max, m_proj
        real(dp), dimension(0:j_max-1), intent(out) :: cos_off1

        integer  :: j
        real(dp) :: j_d, m_d

        cos_off1 = 0.0_dp
        m_d = real(abs(m_proj), dp)

        do j = abs(m_proj), j_max - 1
            j_d = real(j, dp)
            ! <J+1, M | cos theta | J, M> = sqrt( ((J+1)^2 - M^2) / ((2J+1)*(2J+3)) )
            cos_off1(j) = sqrt(max(0.0_dp, ((j_d + 1.0_dp)**2 - m_d**2)) / &
                               max(1.0e-6_dp, (2.0_dp * j_d + 1.0_dp) * (2.0_dp * j_d + 3.0_dp)))
        end do
    end subroutine calc_cos_matrix_elements

    !> \brief 模拟强脉冲诱导的分子转动无场排列复苏动态 <cos^2 theta>(t)
    !> \param[in] mol 分子转子参数
    !> \param[in] laser_i0_wcm2 激光峰值光强 (W/cm^2, 如 1e13 - 5e13)
    !> \param[in] pulse_fwhm_fs 脉冲半高全宽 (fs, 如 50 - 150 fs)
    !> \param[in] temp_k 气体热系综温度 (K, 0K 表示纯基态 J=0)
    !> \param[in] j_max 截断最大角动量量子数 (建议 >= 16)
    !> \param[in] n_time 时间采样步数
    !> \param[in] t_span_ps 总时间跨度 (ps, 建议约为 1.2 * T_rev)
    !> \param[out] t_grid_ps 输出时间轴 (ps)
    !> \param[out] cos2_trace 输出排列度演化曲线 <cos^2 theta>(t)
    !> \param[out] stat 状态码
    subroutine simulate_laser_induced_alignment(mol, laser_i0_wcm2, pulse_fwhm_fs, temp_k, &
                                               j_max, n_time, t_span_ps, t_grid_ps, cos2_trace, stat)
        type(rotor_molecule_t), intent(in)   :: mol
        real(dp), intent(in)                 :: laser_i0_wcm2
        real(dp), intent(in)                 :: pulse_fwhm_fs
        real(dp), intent(in)                 :: temp_k
        integer, intent(in)                  :: j_max, n_time
        real(dp), intent(in)                 :: t_span_ps
        real(dp), dimension(n_time), intent(out) :: t_grid_ps
        real(dp), dimension(n_time), intent(out) :: cos2_trace
        integer, optional, intent(out)       :: stat

        integer  :: it, j, j0, n_j
        real(dp) :: dt_ps, dt_au, t_ps, t_au, sigma_t_au
        real(dp) :: e_field_peak, e_t, v_polar_prefac
        real(dp) :: bolt_weight, z_part, e_rot_k
        real(dp), dimension(0:j_max)   :: cos2_diag, e_free
        real(dp), dimension(0:j_max-2) :: cos2_off2
        complex(dp), dimension(0:j_max) :: psi, psi_new
        real(dp), dimension(n_time) :: cos2_single_j

        if (present(stat)) stat = 0
        if (n_time < 10 .or. j_max < 4 .or. t_span_ps <= 0.0_dp) then
            if (present(stat)) stat = -1
            return
        end if

        n_j = j_max + 1
        dt_ps = t_span_ps / real(n_time - 1, dp)
        dt_au = dt_ps * PS_TO_AU
        sigma_t_au = (pulse_fwhm_fs / 2.354820045_dp) * FS_TO_AU

        ! 激光电场峰值与极化势相互作用强度 V = -1/4 * Delta_alpha * E^2(t) * cos^2(theta)
        e_field_peak = sqrt(laser_i0_wcm2 * WCM2_TO_AU)
        v_polar_prefac = 0.25_dp * mol%delta_alpha_au * (e_field_peak**2)

        do j = 0, j_max
            e_free(j) = mol%b_rot_au * real(j * (j + 1), dp)
        end do

        ! 计算 M = 0 矩阵元
        call calc_cos2_matrix_elements(j_max, 0, cos2_diag, cos2_off2)

        do it = 1, n_time
            t_grid_ps(it) = real(it - 1, dp) * dt_ps
        end do
        cos2_trace = 0.0_dp
        z_part = 0.0_dp

        ! 遍历初始态 J0 (热系综玻尔兹曼加权)
        do j0 = 0, min(8, j_max - 2)
            e_rot_k = (e_free(j0) / CM1_TO_AU) * 1.438777_dp  ! K
            if (temp_k > 0.1_dp) then
                bolt_weight = real(2 * j0 + 1, dp) * exp(-e_rot_k / temp_k)
            else
                if (j0 == 0) then
                    bolt_weight = 1.0_dp
                else
                    bolt_weight = 0.0_dp
                end if
            end if
            if (bolt_weight < 1.0e-5_dp) cycle

            z_part = z_part + bolt_weight

            ! 初始化态矢量
            psi = (0.0_dp, 0.0_dp)
            psi(j0) = (1.0_dp, 0.0_dp)

            do it = 1, n_time
                t_ps = t_grid_ps(it)
                t_au = (t_ps - 0.5_dp) * PS_TO_AU  ! 脉冲中心位于 0.5 ps

                ! 高斯光强包络
                if (abs(t_au) < 5.0_dp * sigma_t_au) then
                    e_t = exp(- (t_au / sigma_t_au)**2)
                else
                    e_t = 0.0_dp
                end if

                ! 演化一步: H = H_0 + V(t)
                psi_new = psi
                do j = 0, j_max
                    psi_new(j) = psi(j) * exp(-EYE * (e_free(j) - v_polar_prefac * e_t * cos2_diag(j)) * dt_au)
                end do
                do j = 2, j_max
                    psi_new(j) = psi_new(j) + &
                        EYE * dt_au * (v_polar_prefac * e_t * cos2_off2(j - 2)) * psi(j - 2)
                end do
                do j = 0, j_max - 2
                    psi_new(j) = psi_new(j) + &
                        EYE * dt_au * (v_polar_prefac * e_t * cos2_off2(j)) * psi(j + 2)
                end do

                ! 归一化
                psi = psi_new / sqrt(max(1.0e-30_dp, sum(abs(psi_new)**2)))

                ! 计算当前期望值 <cos^2 theta> = sum_j |c_j|^2 * diag + 2 * sum_j Re(c_j^* c_{j+2}) * off2
                cos2_single_j(it) = sum((abs(psi)**2) * cos2_diag)
                do j = 0, j_max - 2
                    cos2_single_j(it) = cos2_single_j(it) + &
                        2.0_dp * real(conjg(psi(j + 2)) * psi(j), dp) * cos2_off2(j)
                end do
            end do

            cos2_trace = cos2_trace + bolt_weight * cos2_single_j
        end do

        if (z_part > 1.0e-12_dp) then
            cos2_trace = cos2_trace / z_part
        else
            cos2_trace = 1.0_dp / 3.0_dp
        end if
    end subroutine simulate_laser_induced_alignment

    !> \brief 计算光学离心机 (Optical Centrifuge) 加速分子的末态角动量
    !> \details 电场偏振以恒定角加速度旋转: theta_E(t) = 0.5 * beta * t^2
    !>          分子偶极被绝热锁相在旋转场中，末态转动量子数 J_max ~ (I_mol * beta * t_pulse) / hbar
    !> \param[in] mol 分子结构
    !> \param[in] chirp_beta_thz2 离心机角加速度 beta (THz^2, 通常 0.1 - 1.0 THz/ps)
    !> \param[in] pulse_duration_ps 离心机脉宽 (ps, 通常 50 - 150 ps)
    !> \param[out] j_superrotor 加速达到的分子超转子量子数 J
    subroutine calc_optical_centrifuge_kick(mol, chirp_beta_thz2, pulse_duration_ps, j_superrotor, stat)
        type(rotor_molecule_t), intent(in) :: mol
        real(dp), intent(in)               :: chirp_beta_thz2
        real(dp), intent(in)               :: pulse_duration_ps
        integer, intent(out)               :: j_superrotor
        integer, optional, intent(out)     :: stat

        real(dp) :: omega_rot_final, b_ang_freq

        if (present(stat)) stat = 0
        if (chirp_beta_thz2 <= 0.0_dp .or. pulse_duration_ps <= 0.0_dp) then
            if (present(stat)) stat = -1
            j_superrotor = 0
            return
        end if

        ! 末态旋转角频率: omega_final = beta * T_pulse (rad/ps -> rad/s)
        omega_rot_final = (chirp_beta_thz2 * pulse_duration_ps) * 1.0e12_dp

        ! 分子转动频率与 J 的对应关系: omega = 2 * B_e * c * 2*pi * J
        ! => J = omega / (4 * pi * B_e * c)
        b_ang_freq = 4.0_dp * PI * (mol%b_rot_cm1 * 2.99792458e10_dp)
        j_superrotor = nint(omega_rot_final / b_ang_freq)
        j_superrotor = max(0, j_superrotor)
    end subroutine calc_optical_centrifuge_kick

    !> \brief 评估极端超转子 (Superrotor) 的离心势垒与自解离判据
    !> \details 当极速旋转的离心势垒能量 E_rot(J) = B * J * (J + 1) 超过分子解离能 D_e 时，发生离心破键
    pure subroutine calc_superrotor_dissociation(mol, j_rot, d_e_ev, is_dissociated, rot_energy_ev)
        type(rotor_molecule_t), intent(in) :: mol
        integer, intent(in)                :: j_rot
        real(dp), intent(in)               :: d_e_ev
        logical, intent(out)               :: is_dissociated
        real(dp), intent(out)              :: rot_energy_ev

        real(dp) :: b_ev

        b_ev = mol%b_rot_cm1 * 1.239841984e-4_dp
        rot_energy_ev = b_ev * real(j_rot * (j_rot + 1), dp)

        is_dissociated = (rot_energy_ev >= d_e_ev)
    end subroutine calc_superrotor_dissociation

end module mod_molecular_alignment
