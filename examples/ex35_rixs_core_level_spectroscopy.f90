!> \file ex35_rixs_core_level_spectroscopy.f90
!> \brief Example 35: Resonant Inelastic X-ray Scattering (RIXS) and Core-Level Spectroscopy
!> \author LiHao
program ex35_rixs_core_level_spectroscopy
    use mod_constants, only: dp
    use mod_resonant_xray_scattering
    implicit none

    type(rixs_system_t) :: sys
    real(dp), dimension(1) :: e_core, gamma_core, d_in
    real(dp), dimension(4) :: e_fin, gamma_fin
    real(dp), dimension(4, 1) :: d_out
    real(dp), dimension(5) :: w_scan
    real(dp), dimension(5) :: loss_scan
    real(dp), dimension(5, 5) :: rixs_2d
    real(dp), dimension(0:3) :: vib_intensity
    real(dp) :: w_in, loss_val, sigma_rixs, sigma_xas
    integer  :: iw, il, stat

    print '(A)', "================================================================"
    print '(A)', " Example 35: Resonant Inelastic X-ray Scattering (Cu L3 RIXS)   "
    print '(A)', "================================================================"

    ! 1. Initialize Cu L3-edge system (2p3/2 -> 3d -> dd excitations)
    e_core = [931.5_dp]        ! Cu L3 resonance (eV)
    gamma_core = [0.35_dp]     ! Core hole lifetime HWHM (eV)
    d_in = [1.0_dp]            ! Dipole absorption element

    ! Final states: Elastic, d(xy), d(xz/yz), d(3z^2-r^2)
    e_fin = [0.0_dp, 1.45_dp, 1.80_dp, 2.25_dp]
    gamma_fin = [0.05_dp, 0.08_dp, 0.08_dp, 0.08_dp]
    d_out(1, 1) = 0.85_dp      ! Elastic
    d_out(2, 1) = 0.45_dp      ! d(xy)
    d_out(3, 1) = 0.70_dp      ! d(xz/yz)
    d_out(4, 1) = 0.55_dp      ! d(3z^2-r^2)

    call init_rixs_system(sys, e_init=0.0_dp, e_inter=e_core, gamma_core=gamma_core, &
                          d_in=d_in, e_fin=e_fin, gamma_fin=gamma_fin, d_out=d_out, stat=stat)

    print '(A)', " 1. X-ray Absorption Spectroscopy (XAS) Near Cu L3 Edge:"
    print '(A)', "    Photon Energy w1 (eV) | XAS Cross Section (arb. units)"
    print '(A)', "   -------------------------------------------------------"
    w_scan = [930.5_dp, 931.0_dp, 931.5_dp, 932.0_dp, 932.5_dp]
    do iw = 1, 5
        w_in = w_scan(iw)
        sigma_xas = calc_xas_cross_section(sys, w_in)
        print '(4X, F10.2, 12X, F10.4)', w_in, sigma_xas
    end do

    print '(A)', ""
    print '(A)', " 2. On-Resonance (w1 = 931.5 eV) RIXS dd-Excitation Spectrum:"
    print '(A)', "    Energy Loss Omega (eV) | Assigned Transition | RIXS Intensity"
    print '(A)', "   --------------------------------------------------------------"
    loss_scan = [0.0_dp, 1.45_dp, 1.80_dp, 2.25_dp, 2.80_dp]
    do il = 1, 5
        loss_val = loss_scan(il)
        sigma_rixs = calc_kramers_heisenberg_cross_section(sys, 931.5_dp, loss_val)
        if (il == 1) then
            print '(4X, F10.2, 14X, A15, 4X, F10.3)', loss_val, "Elastic (Ground)", sigma_rixs
        else if (il == 2) then
            print '(4X, F10.2, 14X, A15, 4X, F10.3)', loss_val, "d(xy) orbital", sigma_rixs
        else if (il == 3) then
            print '(4X, F10.2, 14X, A15, 4X, F10.3)', loss_val, "d(xz/yz) orbital", sigma_rixs
        else if (il == 4) then
            print '(4X, F10.2, 14X, A15, 4X, F10.3)', loss_val, "d(3z^2-r^2) orb", sigma_rixs
        else
            print '(4X, F10.2, 14X, A15, 4X, F10.3)', loss_val, "Background", sigma_rixs
        end if
    end do

    print '(A)', ""
    print '(A)', " 3. 2D RIXS Map Cross Section Matrix (Loss vs w1):"
    call calc_rixs_2d_map(sys, 5, w_scan, 5, loss_scan, rixs_2d)
    print '(A)', "    Loss (eV) \\ w1(eV) |  930.5  |  931.0  |  931.5  |  932.0  |  932.5 "
    print '(A)', "   ------------------------------------------------------------------"
    do il = 1, 5
        print '(4X, F6.2, 12X, 5(F8.3, 1X))', loss_scan(il), (rixs_2d(il, iw), iw = 1, 5)
    end do

    print '(A)', ""
    print '(A)', " 4. Electron-Phonon Coupling (Huang-Rhys S = 0.40, w0 = 70 meV):"
    call calc_huang_rhys_vibrational_rixs(omega_0_ev=0.070_dp, s_factor=0.40_dp, &
                                         gamma_core_ev=0.35_dp, detuning_ev=0.0_dp, &
                                         n_max_loss=3, loss_intensity=vib_intensity)
    print '(A)', "    Phonon Loss Peak n | Energy Loss (meV) | Relative RIXS Intensity"
    print '(A)', "   ------------------------------------------------------------------"
    do il = 0, 3
        print '(4X, I3, 16X, F8.1, 14X, ES11.4)', il, real(il, dp) * 70.0_dp, vib_intensity(il)
    end do

    print '(A)', "================================================================"
    print '(A)', " Example 35 completed successfully.                             "
    print '(A)', "================================================================"

end program ex35_rixs_core_level_spectroscopy
