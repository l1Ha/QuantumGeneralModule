!> \file ex29_cold_ion_atom_scattering.f90
!> \brief Example 29: Cold Ion-Atom Hybrid Scattering and Polarization Dynamics
!> \author LiHao
program ex29_cold_ion_atom_scattering
    use mod_constants, only: dp
    use mod_ion_atom_scattering
    implicit none

    type(ion_atom_system_t) :: sys_yb_li, sys_ba_rb
    real(dp) :: e_coll_ev, e_coll_au
    real(dp) :: b_c_au, sigma_l_au, sigma_l_ang2
    real(dp) :: t_limit, heat_rate
    integer  :: stat

    print '(A)', "================================================================"
    print '(A)', " Example 29: Cold Ion-Atom Hybrid Scattering & Micro-heating    "
    print '(A)', "================================================================"

    ! 1. 初始化两组代表性实验杂化冷碰撞体系
    ! Yb+ + 6Li: 极大质量比体系 (m_ion/m_atom ~ 28.4)
    call init_ion_atom_system(sys_yb_li, m_ion_amu=171.0_dp, m_atom_amu=6.015_dp, &
                              charge_e=1.0_dp, alpha_au=164.1_dp, stat=stat)
    ! Ba+ + 87Rb: 接近等质量体系 (m_ion/m_atom ~ 1.58)
    call init_ion_atom_system(sys_ba_rb, m_ion_amu=138.0_dp, m_atom_amu=87.0_dp, &
                              charge_e=1.0_dp, alpha_au=318.8_dp, stat=stat)

    print '(A)', " 1. Characteristic Interaction Scales:"
    print '(A, F8.1, A, ES12.4, A)', "   Yb+/6Li:  R* = ", sys_yb_li%r_star_bohr, &
        " a0, E* = ", sys_yb_li%e_star_kelvin * 1.0e6_dp, " uK"
    print '(A, F8.1, A, ES12.4, A)', "   Ba+/87Rb: R* = ", sys_ba_rb%r_star_bohr, &
        " a0, E* = ", sys_ba_rb%e_star_kelvin * 1.0e6_dp, " uK"

    print '(A)', " 2. Universal Langevin Reaction Rate Coefficients:"
    print '(A, ES14.4, A)', "   K_Langevin (Yb+/6Li)  = ", sys_yb_li%k_langevin_cm3_s, " cm^3/s"
    print '(A, ES14.4, A)', "   K_Langevin (Ba+/87Rb) = ", sys_ba_rb%k_langevin_cm3_s, " cm^3/s"

    ! 3. 碰撞能谱与经典临界碰撞参数
    e_coll_ev = 1.0e-3_dp  ! 1 meV
    e_coll_au = e_coll_ev / 27.211386_dp
    b_c_au = calc_langevin_critical_impact_parameter(sys_yb_li, e_coll_au)
    sigma_l_au = calc_langevin_cross_section(sys_yb_li, e_coll_au)
    sigma_l_ang2 = sigma_l_au * (0.529177210903_dp**2)

    print '(A)', " 3. Collision Cross Sections at E = 1 meV:"
    print '(A, F10.2, A, F10.2, A)', "   Critical Impact Param b_c = ", b_c_au, " a0"
    print '(A, ES14.4, A)', "   Langevin Capture Cross Section = ", sigma_l_ang2, " Angstrom^2"

    ! 4. 射频场中微运动发热与平衡极限
    call calc_rf_micromotion_heating(trap_q=0.20_dp, m_ion=171.0_dp, m_atom=6.015_dp, &
                                    temp_atom_k=1.0e-6_dp, coll_rate_hz=200.0_dp, &
                                    t_limit_k=t_limit, heating_rate_k_s=heat_rate)

    print '(A)', " 4. Paul Trap Micromotion Heating (Yb+/6Li, q = 0.20):"
    print '(A, ES12.4, A, ES12.4, A)', "   Equilibrium T_limit = ", t_limit * 1.0e6_dp, &
        " uK, Initial dE/dt = ", heat_rate, " K/s"

    print '(A)', "================================================================"
    print '(A)', " Example 29 completed successfully.                             "
    print '(A)', "================================================================"

end program ex29_cold_ion_atom_scattering
