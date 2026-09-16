!> \file mod_resonant_xray_scattering.f90
!> \brief 共振非弹性 X 射线散射 (RIXS) 与内壳层光谱动力学模块
!> \details 基于 Kramers-Heisenberg 二阶微扰公式，求解 X 射线吸收谱 (XAS)、发射谱 (XES)
!>          以及 2D 能量损失共振非弹性散射截面，涵盖核心空穴 Auger 寿命展宽、
!>          拉曼/荧光过渡区行为及电子-声子耦合 Huang-Rhys 振动 Franck-Condon 级数。
!> \author LiHao
module mod_resonant_xray_scattering
    use mod_constants, only: dp, PI
    implicit none
    private

    public :: rixs_system_t
    public :: init_rixs_system
    public :: calc_xas_cross_section
    public :: calc_kramers_heisenberg_cross_section
    public :: calc_rixs_2d_map
    public :: calc_huang_rhys_vibrational_rixs

    !> RIXS 多能级系统配置结构体
    type :: rixs_system_t
        real(dp) :: e_initial                        !< 初态基态能量 (eV)
        integer  :: n_intermediate                   !< 中间核心激发态数目
        real(dp), allocatable :: e_intermediate(:)   !< 中间态能量 E_m (eV)
        real(dp), allocatable :: gamma_core(:)       !< 中间态核心空穴寿命展宽 HWHM (eV)
        real(dp), allocatable :: d_in(:)             !< 入射偶极跃迁矩阵元 <m|D_1|i>
        integer  :: n_final                          !< 终态低能激发态数目
        real(dp), allocatable :: e_final(:)          !< 终态能量 E_f (eV)
        real(dp), allocatable :: gamma_final(:)      !< 终态本征寿命/仪器展宽 HWHM (eV)
        real(dp), allocatable :: d_out(:, :)         !< 出射偶极跃迁矩阵元 <f|D_2^dagger|m> (n_final, n_intermediate)
    end type rixs_system_t

contains

    !> \brief 初始化 RIXS 能级系统与偶极跃迁参数
    !> \param[out] sys RIXS 系统对象
    !> \param[in] e_init 初态能量 (eV)
    !> \param[in] e_inter 中间核心能级数组 (eV)
    !> \param[in] gamma_core 核心空穴展宽数组 (eV)
    !> \param[in] d_in 入射跃迁矩数组
    !> \param[in] e_fin 终态能级数组 (eV)
    !> \param[in] gamma_fin 终态展宽数组 (eV)
    !> \param[in] d_out 出射跃迁矩矩阵 (n_final x n_inter)
    !> \param[out] stat 状态码 (0 成功)
    subroutine init_rixs_system(sys, e_init, e_inter, gamma_core, d_in, &
                                e_fin, gamma_fin, d_out, stat)
        type(rixs_system_t), intent(out) :: sys
        real(dp), intent(in)             :: e_init
        real(dp), intent(in)             :: e_inter(:), gamma_core(:), d_in(:)
        real(dp), intent(in)             :: e_fin(:), gamma_fin(:), d_out(:, :)
        integer, optional, intent(out)   :: stat

        integer :: ni, nf

        if (present(stat)) stat = 0
        ni = size(e_inter)
        nf = size(e_fin)

        if (ni < 1 .or. nf < 1 .or. size(gamma_core) /= ni .or. &
            size(d_in) /= ni .or. size(gamma_fin) /= nf .or. &
            size(d_out, 1) /= nf .or. size(d_out, 2) /= ni) then
            if (present(stat)) stat = -1
            return
        end if

        sys%e_initial = e_init
        sys%n_intermediate = ni
        sys%n_final = nf

        allocate(sys%e_intermediate(ni))
        allocate(sys%gamma_core(ni))
        allocate(sys%d_in(ni))
        allocate(sys%e_final(nf))
        allocate(sys%gamma_final(nf))
        allocate(sys%d_out(nf, ni))

        sys%e_intermediate = e_inter
        sys%gamma_core = gamma_core
        sys%d_in = d_in
        sys%e_final = e_fin
        sys%gamma_final = gamma_fin
        sys%d_out = d_out
    end subroutine init_rixs_system

    !> \brief 计算一阶 X 射线吸收截面 XAS (Optical theorem)
    !> \param[in] sys RIXS 系统对象
    !> \param[in] omega_in_ev 入射光子能量 (eV)
    !> \return 吸收截面强度 (Lorentzian 展宽)
    pure function calc_xas_cross_section(sys, omega_in_ev) result(sigma_xas)
        type(rixs_system_t), intent(in) :: sys
        real(dp), intent(in)            :: omega_in_ev
        real(dp)                        :: sigma_xas

        integer  :: m
        real(dp) :: delta_e, gamma_m

        sigma_xas = 0.0_dp
        if (.not. allocated(sys%e_intermediate)) return

        do m = 1, sys%n_intermediate
            delta_e = (sys%e_initial + omega_in_ev) - sys%e_intermediate(m)
            gamma_m = max(1.0e-5_dp, sys%gamma_core(m))
            ! 洛伦兹吸收线型: |<m|D|i>|^2 * (Gamma / 2pi) / (Delta_E^2 + Gamma^2)
            sigma_xas = sigma_xas + (sys%d_in(m)**2) * (gamma_m / PI) / (delta_e**2 + gamma_m**2)
        end do
    end function calc_xas_cross_section

    !> \brief 计算给定入射能量与能量损失处的 Kramers-Heisenberg 二阶散射截面
    !> \details d^2 sigma / (dOmega d(omega_loss)) ~ sum_f |sum_m <f|D_2|m><m|D_1|i> / (E_i + w1 - E_m + i*Gamma_m)|^2
    !>          * L(omega_loss - (E_f - E_i), gamma_f)
    !> \param[in] sys RIXS 系统对象
    !> \param[in] omega_in_ev 入射光子能量 w1 (eV)
    !> \param[in] omega_loss_ev 能量损失 Omega = w1 - w2 (eV)
    !> \return 二阶非弹性散射强度 (a.u.)
    pure function calc_kramers_heisenberg_cross_section(sys, omega_in_ev, omega_loss_ev) result(sigma_rixs)
        type(rixs_system_t), intent(in) :: sys
        real(dp), intent(in)            :: omega_in_ev, omega_loss_ev
        real(dp)                        :: sigma_rixs

        integer  :: f, m
        real(dp) :: amp_re, amp_im, denom, denom_re, denom_im
        real(dp) :: loss_diff, gamma_f, lineshape

        sigma_rixs = 0.0_dp
        if (.not. allocated(sys%e_intermediate) .or. .not. allocated(sys%e_final)) return

        ! 对每个终态 |f> 进行非相干求和
        do f = 1, sys%n_final
            amp_re = 0.0_dp
            amp_im = 0.0_dp

            ! 对中间态 |m> 进行相干量子干涉叠加
            do m = 1, sys%n_intermediate
                denom_re = (sys%e_initial + omega_in_ev) - sys%e_intermediate(m)
                denom_im = max(1.0e-5_dp, sys%gamma_core(m))
                denom = denom_re**2 + denom_im**2

                ! 跃迁偶极乘积 M_fm = <f|D_2|m> * <m|D_1|i>
                ! (M_fm) / (denom_re + i * denom_im) = M_fm * (denom_re - i * denom_im) / denom
                amp_re = amp_re + (sys%d_out(f, m) * sys%d_in(m)) * denom_re / denom
                amp_im = amp_im - (sys%d_out(f, m) * sys%d_in(m)) * denom_im / denom
            end do

            ! 终态能量损失的洛伦兹展宽 lineshape
            loss_diff = omega_loss_ev - (sys%e_final(f) - sys%e_initial)
            gamma_f = max(1.0e-5_dp, sys%gamma_final(f))
            lineshape = (gamma_f / PI) / (loss_diff**2 + gamma_f**2)

            ! 截面累加: |Amplitude|^2 * Lineshape
            sigma_rixs = sigma_rixs + (amp_re**2 + amp_im**2) * lineshape
        end do
    end function calc_kramers_heisenberg_cross_section

    !> \brief 生成 2D RIXS 散射谱强度矩阵
    !> \param[in] sys RIXS 系统对象
    !> \param[in] n_in 入射能量网格点数
    !> \param[in] omega_in_grid 入射能量网格 (eV)
    !> \param[in] n_loss 能量损失网格点数
    !> \param[in] omega_loss_grid 能量损失网格 (eV)
    !> \param[out] rixs_map 2D 强度谱矩阵 (n_loss x n_in)
    subroutine calc_rixs_2d_map(sys, n_in, omega_in_grid, n_loss, omega_loss_grid, rixs_map)
        type(rixs_system_t), intent(in) :: sys
        integer, intent(in)             :: n_in, n_loss
        real(dp), intent(in)            :: omega_in_grid(n_in)
        real(dp), intent(in)            :: omega_loss_grid(n_loss)
        real(dp), intent(out)           :: rixs_map(n_loss, n_in)

        integer :: i_in, i_loss

        do i_in = 1, n_in
            do i_loss = 1, n_loss
                rixs_map(i_loss, i_in) = calc_kramers_heisenberg_cross_section(sys, &
                    omega_in_grid(i_in), omega_loss_grid(i_loss))
            end do
        end do
    end subroutine calc_rixs_2d_map

    !> \brief 计算电子-声子耦合 Huang-Rhys 振动 Franck-Condon RIXS 级数强度
    !> \details 基于 Ament-van den Brink 模型 (RMP 83, 705, 2011):
    !>          计算振动级数能量损失 n * omega_0 处的相对峰值强度:
    !>          A_n = sum_{nu=0}^infty <n|D(-d)|nu> <nu|D(d)|0> / (Delta - nu * omega_0 + i*Gamma)
    !>          其中重叠矩阵元采用相联拉盖尔多项式 L_p^(alpha) 严格解析计算
    !> \param[in] omega_0_ev 声子/振动特征频率 (eV)
    !> \param[in] s_factor Huang-Rhys 无量纲电子-声子耦合常数 S
    !> \param[in] gamma_core_ev 核心空穴展宽 HWHM (eV)
    !> \param[in] detuning_ev 入射光偏离共振中心能量 Delta = w1 - w_res (eV)
    !> \param[in] n_max_loss 最大计算振动态阶数 n
    !> \param[out] loss_intensity 各阶振动损失峰 (n = 0, 1, ..., n_max_loss) 相对截面
    subroutine calc_huang_rhys_vibrational_rixs(omega_0_ev, s_factor, gamma_core_ev, detuning_ev, &
                                                n_max_loss, loss_intensity)
        real(dp), intent(in)  :: omega_0_ev, s_factor, gamma_core_ev, detuning_ev
        integer,  intent(in)  :: n_max_loss
        real(dp), intent(out) :: loss_intensity(0:n_max_loss)

        integer  :: n, nu, k
        real(dp) :: d_shift, amp_re, amp_im, fc_in, fc_out, denom, denom_re, denom_im
        real(dp) :: fact(0:30)

        loss_intensity = 0.0_dp
        if (n_max_loss < 0 .or. omega_0_ev <= 0.0_dp) return

        fact(0) = 1.0_dp
        do k = 1, 30
            fact(k) = fact(k - 1) * real(k, dp)
        end do

        d_shift = sqrt(max(0.0_dp, s_factor))

        ! 遍历每个终态振动阶数 n
        do n = 0, min(10, n_max_loss)
            amp_re = 0.0_dp
            amp_im = 0.0_dp

            ! 中间核心能级振动量子数 nu 截断至 25 阶
            do nu = 0, 25
                ! <nu|D(d)|0> = (d^nu / sqrt(nu!)) * exp(-S/2)
                fc_in = (d_shift**nu / sqrt(fact(nu))) * exp(-0.5_dp * s_factor)

                ! <n|D(-d)|nu> = <nu|D(d)|n>
                if (nu >= n) then
                    fc_out = sqrt(fact(n) / fact(nu)) * (d_shift**(nu - n)) * exp(-0.5_dp * s_factor) * &
                             eval_assoc_laguerre(n, nu - n, s_factor, fact)
                else
                    fc_out = ((-1.0_dp)**(n - nu)) * sqrt(fact(nu) / fact(n)) * (d_shift**(n - nu)) * &
                             exp(-0.5_dp * s_factor) * eval_assoc_laguerre(nu, n - nu, s_factor, fact)
                end if

                denom_re = detuning_ev - real(nu, dp) * omega_0_ev
                denom_im = max(1.0e-5_dp, gamma_core_ev)
                denom = denom_re**2 + denom_im**2

                amp_re = amp_re + (fc_out * fc_in) * denom_re / denom
                amp_im = amp_im - (fc_out * fc_in) * denom_im / denom
            end do

            loss_intensity(n) = amp_re**2 + amp_im**2
        end do
    end subroutine calc_huang_rhys_vibrational_rixs

    !> \brief 辅助纯函数: 计算相联拉盖尔多项式 L_p^(alpha)(x)
    pure function eval_assoc_laguerre(p, alpha, x, fact) result(val)
        integer,  intent(in) :: p, alpha
        real(dp), intent(in) :: x
        real(dp), intent(in) :: fact(0:30)
        real(dp)             :: val

        integer :: m

        val = 0.0_dp
        do m = 0, p
            val = val + ((-1.0_dp)**m) * (fact(p + alpha) / &
                        (fact(p - m) * fact(alpha + m) * fact(m))) * (x**m)
        end do
    end function eval_assoc_laguerre

end module mod_resonant_xray_scattering
