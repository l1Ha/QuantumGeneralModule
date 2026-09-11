!> \brief 物理常数与单位转换核心模块
!> \details 提供 CODATA 推荐的基础常数、数学常数及高精度原子单位（a.u.）与国际单位制（SI）双向换算。
!> \author LiHao
module mod_constants
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32, int64
    implicit none
    private

    ! 公共常数导出
    public :: dp, int32, int64
    public :: PI, TWOPI, HALFPI, SQRTPI, EYE
    public :: C_LIGHT, HBAR, H_PLANCK, M_E, CHARGE_E, EPS0, KB, AMU2AU

    ! 换算因子导出
    public :: AU2S, S2AU, AU2FS, FS2AU, AU2PS, PS2AU
    public :: AU2EV, EV2AU, AU2CM, CM2AU, AU2J, J2AU, AU2K, K2AU
    public :: AU2M, M2AU, AU2ANG, ANG2AU, AU2NM, NM2AU
    public :: AU2VM, VM2AU, AU2MV_CM, MV_CM2AU
    public :: AU2DEBYE, DEBYE2AU, AU2W_CM2, W_CM2AU

    ! 纯函数转换接口导出
    public :: to_au, from_au

    ! 数学常数
    real(dp), parameter :: PI     = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: TWOPI  = 2.0_dp * PI
    real(dp), parameter :: HALFPI = 0.5_dp * PI
    real(dp), parameter :: SQRTPI = 1.77245385090551602729816748334114518_dp
    complex(dp), parameter :: EYE = (0.0_dp, 1.0_dp)

    ! 基础物理常数 (CODATA 2018/2022)
    real(dp), parameter :: C_LIGHT  = 2.99792458e8_dp         !< 真空光速 (m/s)
    real(dp), parameter :: HBAR     = 1.054571817e-34_dp       !< 约化普朗克常数 (J*s)
    real(dp), parameter :: H_PLANCK = 6.62607015e-34_dp        !< 普朗克常数 (J*s)
    real(dp), parameter :: M_E      = 9.1093837015e-31_dp      !< 电子静止质量 (kg)
    real(dp), parameter :: CHARGE_E = 1.602176634e-19_dp       !< 基本电荷 (C)
    real(dp), parameter :: EPS0     = 8.8541878128e-12_dp      !< 真空介电常数 (F/m)
    real(dp), parameter :: KB       = 1.380649e-23_dp          !< 玻尔兹曼常数 (J/K)
    real(dp), parameter :: AMU2AU   = 1822.888486209_dp        !< 1 amu 对应的原子单位质量 (m_e)

    ! 时间单位转换
    real(dp), parameter :: AU2S     = 2.4188843265857e-17_dp   !< a.u. 时间 -> 秒 (s)
    real(dp), parameter :: S2AU     = 1.0_dp / AU2S            !< 秒 (s) -> a.u.
    real(dp), parameter :: AU2FS    = 2.4188843265857e-2_dp    !< a.u. -> 飞秒 (fs)
    real(dp), parameter :: FS2AU    = 1.0_dp / AU2FS           !< 飞秒 (fs) -> a.u.
    real(dp), parameter :: AU2PS    = 2.4188843265857e-5_dp    !< a.u. -> 皮秒 (ps)
    real(dp), parameter :: PS2AU    = 1.0_dp / AU2PS           !< 皮秒 (ps) -> a.u.

    ! 能量单位转换
    real(dp), parameter :: AU2EV    = 27.211386245988_dp       !< 1 Hartree -> eV
    real(dp), parameter :: EV2AU    = 1.0_dp / AU2EV           !< eV -> a.u.
    real(dp), parameter :: AU2CM    = 219474.63136320_dp       !< 1 Hartree -> cm^-1 (波数)
    real(dp), parameter :: CM2AU    = 1.0_dp / AU2CM           !< cm^-1 -> a.u.
    real(dp), parameter :: AU2J     = 4.3597447222071e-18_dp   !< 1 Hartree -> 焦耳 (J)
    real(dp), parameter :: J2AU     = 1.0_dp / AU2J            !< J -> a.u.
    real(dp), parameter :: AU2K     = 3.1577502480407e5_dp     !< 1 Hartree -> 开尔文 (K)
    real(dp), parameter :: K2AU     = 1.0_dp / AU2K            !< K -> a.u.

    ! 长度单位转换
    real(dp), parameter :: AU2M     = 0.529177210903e-10_dp    !< 1 Bohr -> 米 (m)
    real(dp), parameter :: M2AU     = 1.0_dp / AU2M            !< 米 (m) -> a.u.
    real(dp), parameter :: AU2ANG   = 0.529177210903_dp        !< 1 Bohr -> 埃 (Angstrom)
    real(dp), parameter :: ANG2AU   = 1.0_dp / AU2ANG          !< 埃 -> a.u.
    real(dp), parameter :: AU2NM    = 0.0529177210903_dp       !< 1 Bohr -> 纳米 (nm)
    real(dp), parameter :: NM2AU    = 1.0_dp / AU2NM           !< 纳米 -> a.u.

    ! 电场强度单位转换
    real(dp), parameter :: AU2VM    = 5.14220674763e11_dp      !< 1 a.u. 电场 -> V/m
    real(dp), parameter :: VM2AU    = 1.0_dp / AU2VM           !< V/m -> a.u.
    real(dp), parameter :: AU2MV_CM = 5142.20674763_dp         !< 1 a.u. 电场 -> MV/cm
    real(dp), parameter :: MV_CM2AU = 1.0_dp / AU2MV_CM        !< MV/cm -> a.u.

    ! 偶极矩单位转换
    real(dp), parameter :: AU2DEBYE = 2.541746473_dp           !< 1 a.u. 偶极矩 -> Debye
    real(dp), parameter :: DEBYE2AU = 1.0_dp / AU2DEBYE        !< Debye -> a.u.

    ! 激光光强单位转换 (I = 1/2 * eps0 * c * E^2)
    real(dp), parameter :: AU2W_CM2 = 3.5094452e16_dp          !< 1 a.u. 光强 -> W/cm^2
    real(dp), parameter :: W_CM2AU  = 1.0_dp / AU2W_CM2        !< W/cm^2 -> a.u.

contains

    !> \brief 统一将特定单位的标量值换算为原子单位 a.u.
    pure function to_au(val, unit_name) result(val_au)
        real(dp), intent(in) :: val
        character(len=*), intent(in) :: unit_name
        real(dp) :: val_au

        select case(trim(adjustl(unit_name)))
        case('fs')
            val_au = val * FS2AU
        case('ps')
            val_au = val * PS2AU
        case('s')
            val_au = val * S2AU
        case('ev', 'eV')
            val_au = val * EV2AU
        case('cm-1', 'cm^-1')
            val_au = val * CM2AU
        case('j', 'J')
            val_au = val * J2AU
        case('k', 'K')
            val_au = val * K2AU
        case('ang', 'angstrom', 'Angstrom')
            val_au = val * ANG2AU
        case('nm')
            val_au = val * NM2AU
        case('m')
            val_au = val * M2AU
        case('amu')
            val_au = val * AMU2AU
        case('v/m', 'V/m')
            val_au = val * VM2AU
        case('mv/cm', 'MV/cm')
            val_au = val * MV_CM2AU
        case('debye', 'Debye')
            val_au = val * DEBYE2AU
        case('w/cm2', 'W/cm^2')
            val_au = val * W_CM2AU
        case default
            val_au = val
        end select
    end function to_au

    !> \brief 统一将原子单位 a.u. 标量值转换为指定实际物理单位
    pure function from_au(val_au, unit_name) result(val)
        real(dp), intent(in) :: val_au
        character(len=*), intent(in) :: unit_name
        real(dp) :: val

        select case(trim(adjustl(unit_name)))
        case('fs')
            val = val_au * AU2FS
        case('ps')
            val = val_au * AU2PS
        case('s')
            val = val_au * AU2S
        case('ev', 'eV')
            val = val_au * AU2EV
        case('cm-1', 'cm^-1')
            val = val_au * AU2CM
        case('j', 'J')
            val = val_au * AU2J
        case('k', 'K')
            val = val_au * AU2K
        case('ang', 'angstrom', 'Angstrom')
            val = val_au * AU2ANG
        case('nm')
            val = val_au * AU2NM
        case('m')
            val = val_au * AU2M
        case('amu')
            val = val_au / AMU2AU
        case('v/m', 'V/m')
            val = val_au * AU2VM
        case('mv/cm', 'MV/cm')
            val = val_au * AU2MV_CM
        case('debye', 'Debye')
            val = val_au * AU2DEBYE
        case('w/cm2', 'W/cm^2')
            val = val_au * AU2W_CM2
        case default
            val = val_au
        end select
    end function from_au

end module mod_constants
