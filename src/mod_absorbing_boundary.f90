!> \brief 复吸收势边界（CAP）与量子概率通量计算模块
!> \details 用于消除有限网格边界波包反射伪影，并高精度计算光解离通量、逃逸碎片产额与时间分辨光电离几率。
!> \author LiHao
module mod_absorbing_boundary
    use mod_constants, only: dp, PI, EYE
    implicit none
    private

    public :: CAP_SIN2, CAP_POLYNOMIAL
    public :: absorbing_boundary_t
    public :: cap_init
    public :: cap_evaluate
    public :: cap_apply_mask
    public :: calculate_probability_flux
    public :: calculate_norm_inside

    integer, parameter :: CAP_SIN2       = 1
    integer, parameter :: CAP_POLYNOMIAL = 2

    !> \brief 复吸收边界配置与状态结构体
    type :: absorbing_boundary_t
        integer  :: cap_type = CAP_SIN2    !< 吸收势函数形式
        real(dp) :: r_start  = 30.0_dp     !< 吸收起始位置 R_0 (a.u.)
        real(dp) :: r_end    = 40.0_dp     !< 网格最大边界 R_max (a.u.)
        real(dp) :: strength = 1.0_dp      !< 吸收强度系数 A (a.u.)
        integer  :: power    = 2           !< 多项式吸收幂次 (通常为 2 或 3)
    end type absorbing_boundary_t

contains

    !> \brief 初始化复吸收边界参数
    subroutine cap_init(r_start, r_end, strength, cap_type, cap_obj)
        real(dp), intent(in) :: r_start, r_end, strength
        integer, intent(in), optional :: cap_type
        type(absorbing_boundary_t), intent(out) :: cap_obj

        cap_obj%r_start = r_start
        cap_obj%r_end = r_end
        cap_obj%strength = strength
        if (present(cap_type)) cap_obj%cap_type = cap_type
    end subroutine cap_init

    !> \brief 计算位置 r 处的虚部复吸收势 -i * W(r) (Hartree)
    pure function cap_evaluate(r, cap_obj) result(v_abs)
        real(dp), intent(in) :: r
        type(absorbing_boundary_t), intent(in) :: cap_obj
        complex(dp) :: v_abs
        real(dp) :: frac

        v_abs = (0.0_dp, 0.0_dp)
        if (r <= cap_obj%r_start) return

        frac = (r - cap_obj%r_start) / max(1.0e-12_dp, (cap_obj%r_end - cap_obj%r_start))
        frac = min(1.0_dp, max(0.0_dp, frac))

        select case(cap_obj%cap_type)
        case(CAP_SIN2)
            ! W(r) = A * sin^2(pi/2 * frac)
            v_abs = -EYE * cap_obj%strength * (sin(0.5_dp * PI * frac)**2)
        case(CAP_POLYNOMIAL)
            ! W(r) = A * frac^power
            v_abs = -EYE * cap_obj%strength * (frac**cap_obj%power)
        case default
            v_abs = -EYE * cap_obj%strength * (sin(0.5_dp * PI * frac)**2)
        end select
    end function cap_evaluate

    !> \brief 步进中应用平滑吸收掩膜衰减波包: psi = psi * exp(-W * dt)
    subroutine cap_apply_mask(r_grid, dt, cap_obj, psi)
        real(dp), intent(in) :: r_grid(:)
        real(dp), intent(in) :: dt
        type(absorbing_boundary_t), intent(in) :: cap_obj
        complex(dp), intent(inout) :: psi(:)
        integer :: i, n
        real(dp) :: frac, damp

        n = size(r_grid)
        do i = 1, n
            if (r_grid(i) > cap_obj%r_start) then
                frac = (r_grid(i) - cap_obj%r_start) / max(1.0e-12_dp, (cap_obj%r_end - cap_obj%r_start))
                frac = min(1.0_dp, max(0.0_dp, frac))
                select case(cap_obj%cap_type)
                case(CAP_SIN2)
                    damp = exp(-cap_obj%strength * dt * (sin(0.5_dp * PI * frac)**2))
                case default
                    damp = exp(-cap_obj%strength * dt * (frac**cap_obj%power))
                end select
                psi(i) = psi(i) * damp
            end if
        end do
    end subroutine cap_apply_mask

    !> \brief 计算渐近监测面处的量子概率流密度 J(r_d, t) = (1/mu) * Im[psi* * d_psi/dr]
    pure function calculate_probability_flux(r_grid, mass, psi, idx_detect) result(flux)
        real(dp), intent(in) :: r_grid(:)
        real(dp), intent(in) :: mass
        complex(dp), intent(in) :: psi(:)
        integer, intent(in) :: idx_detect
        real(dp) :: flux
        real(dp) :: dr
        complex(dp) :: dpsi_dr

        flux = 0.0_dp
        if (idx_detect <= 1 .or. idx_detect >= size(psi)) return

        ! 二阶中心差分近似导数
        dr = r_grid(idx_detect + 1) - r_grid(idx_detect - 1)
        dpsi_dr = (psi(idx_detect + 1) - psi(idx_detect - 1)) / dr

        flux = (1.0_dp / mass) * aimag(conjg(psi(idx_detect)) * dpsi_dr)
    end function calculate_probability_flux

    !> \brief 计算未吸收内部反应区（R < R_0）的存活总范数: P_inside = int_0^R0 |psi|^2 dr
    pure function calculate_norm_inside(r_grid, psi, r_boundary) result(norm_in)
        real(dp), intent(in) :: r_grid(:)
        complex(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: r_boundary
        real(dp) :: norm_in
        integer :: i, n
        real(dp) :: dr

        norm_in = 0.0_dp
        n = size(r_grid)
        if (n <= 1) return

        dr = r_grid(2) - r_grid(1)
        do i = 1, n
            if (r_grid(i) <= r_boundary) then
                norm_in = norm_in + (abs(psi(i))**2) * dr
            end if
        end do
    end function calculate_norm_inside

end module mod_absorbing_boundary
