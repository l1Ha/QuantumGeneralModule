! ==============================================================================
! GeneralModule Example 25: 2D Corrugated Surface Diffraction & SAR
!
! Scientific Engineering Application:
!  1. He / LiF(001) quantum surface scattering under Hard Corrugated Surface (HCS)
!  2. 2D diffraction channels G = (m, n) and Eikonal diffraction probabilities
!  3. Debye-Waller thermal phonon attenuation factor DW(T) from 50 K to 800 K
!  4. Selective Adsorption Resonance (SAR) Fano profile scan vs incident angle
!
! Reference:
!   Boato et al., J. Phys. C: Solid State Phys. 6, L394 (1973).
!
! Standard: Fortran 2008 (Pure Fortran, zero external dependencies)
! ==============================================================================
program ex25_surface_corrugated_diffraction
    use mod_constants, only: dp, AU2ANG, AU2EV
    use mod_surface_scattering
    implicit none

    type(surface_lattice_t)   :: lif
    type(surface_potential_t) :: he_pot
    type(diffraction_beam_t)  :: channels(25)
    integer  :: n_ch, n_bound, i, u_diff, u_sar, step
    real(dp) :: e_inc, th_inc, phi_inc, m_he
    real(dp) :: temp_k, dw_val, th_deg, fano_val, de_mev
    logical  :: is_near

    print *, "================================================================"
    print *, " Example 25: 2D Corrugated Surface Diffraction & SAR Resonances "
    print *, "================================================================"

    ! 1. Initialize LiF(001) surface lattice and He-LiF Morse potential
    m_he = 4.0026_dp
    call init_surface_lattice(lif, ax_ang=2.84_dp, ay_ang=2.84_dp, &
                              zeta_x_ang=0.06_dp, zeta_y_ang=0.06_dp, &
                              m_sub_amu=25.94_dp, debye_temp_k=730.0_dp)

    call init_surface_potential_morse(he_pot, well_depth_mev=7.5_dp, &
                                      range_inv_ang=1.1_dp, mass_amu=m_he, &
                                      n_bound=n_bound)

    write(*, '(A, F8.3, A, F8.3, A)') " Surface Lattice Constants   : ax = ", lif%ax_bohr * AU2ANG, &
                                       " A, ay = ", lif%ay_bohr * AU2ANG, " A"
    write(*, '(A, F8.3, A)')          " Surface Corrugation zeta    : ", lif%zeta_x_bohr * AU2ANG, " A"
    write(*, '(A, I4, A, F8.2, A)')   " Morse Bound States          : ", n_bound, &
                                       " (Well Depth D = ", he_pot%well_depth_au * AU2EV * 1000.0_dp, " meV)"
    print *, "----------------------------------------------------------------"

    ! 2. 2D Diffraction Probabilities at E = 20 meV, theta = 40 deg
    e_inc   = 0.020_dp ! 20 meV
    th_inc  = 40.0_dp  ! 40 deg
    phi_inc = 0.0_dp   ! along [100] azimuth

    call calc_hcs_diffraction_probabilities(lif, mass_amu=m_he, energy_ev=e_inc, &
                                            theta_i_deg=th_inc, phi_i_deg=phi_inc, &
                                            max_order=1, n_channels=n_ch, channels=channels)

    open(newunit=u_diff, file="ex25_surface_diffraction.dat", status="replace", action="write")
    write(u_diff, '(A)') "# GeneralModule Example 25: He/LiF(001) Surface Diffraction Beams"
    write(u_diff, '(A)') "# E_inc = 20 meV, theta_i = 40 deg, phi_i = 0 deg"
    write(u_diff, '(A)') "# Col 1: Order m"
    write(u_diff, '(A)') "# Col 2: Order n"
    write(u_diff, '(A)') "# Col 3: Open status (1=open, 0=closed)"
    write(u_diff, '(A)') "# Col 4: Outgoing polar angle theta_G (deg)"
    write(u_diff, '(A)') "# Col 5: Outgoing azimuthal angle phi_G (deg)"
    write(u_diff, '(A)') "# Col 6: Diffraction probability P_G"

    do i = 1, n_ch
        write(u_diff, '(2I4, I3, 3ES16.6)') channels(i)%m, channels(i)%n, &
            merge(1, 0, channels(i)%is_open), channels(i)%theta_f_deg, &
            channels(i)%phi_f_deg, channels(i)%probability
    end do
    close(u_diff)
    print *, " Saved 2D diffraction beam spectrum to ex25_surface_diffraction.dat"

    ! 3. Selective Adsorption Resonance (SAR) & Debye-Waller vs Incident Angle / Temperature
    open(newunit=u_sar, file="ex25_surface_sar_scan.dat", status="replace", action="write")
    write(u_sar, '(A)') "# GeneralModule Example 25: SAR Fano Resonances & Debye-Waller Scan"
    write(u_sar, '(A)') "# Col 1: Incident angle theta_i (deg)"
    write(u_sar, '(A)') "# Col 2: Detuning Delta E (meV) to bound state v=0, G=(-1, 0)"
    write(u_sar, '(A)') "# Col 3: Specular Fano resonance modulation ratio"
    write(u_sar, '(A)') "# Col 4: Surface temperature T (K)"
    write(u_sar, '(A)') "# Col 5: Debye-Waller attenuation factor DW(T)"

    do step = 0, 100
        th_deg = 50.0_dp + real(step, dp) * 0.10_dp ! 50 deg to 60 deg
        temp_k = 50.0_dp + real(step, dp) * 7.50_dp ! 50 K to 800 K

        call calc_selective_adsorption_resonance(lif, he_pot, mass_amu=m_he, &
                                                 energy_ev=e_inc, theta_i_deg=th_deg, &
                                                 phi_i_deg=phi_inc, m_res=-1, n_res=0, &
                                                 v_bound=0, is_near_res=is_near, &
                                                 delta_e_mev=de_mev, &
                                                 fano_specular_ratio=fano_val)

        dw_val = calc_surface_debye_waller(lif, mass_amu=m_he, k_iz_au=1.8_dp, &
                                           kz_g_au=1.8_dp, temp_k=temp_k)

        write(u_sar, '(5ES16.6)') th_deg, de_mev, fano_val, temp_k, dw_val
    end do
    close(u_sar)
    print *, " Saved SAR Fano resonance & Debye-Waller scan to ex25_surface_sar_scan.dat"
    print *, " Example 25 completed successfully."

end program ex25_surface_corrugated_diffraction
