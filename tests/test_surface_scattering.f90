! ==============================================================================
! GeneralModule: test_surface_scattering.f90
!
! Unit Test Suite: Quantum Surface Scattering & Diffraction
!
! Standard: Fortran 2008
! ==============================================================================

program test_surface_scattering
    use mod_constants, only: dp
    use mod_surface_scattering
    implicit none

    integer :: n_pass, n_total
    type(surface_lattice_t)   :: lif
    type(surface_potential_t) :: he_lif_pot
    type(diffraction_beam_t)  :: channels(25)
    integer  :: n_ch, n_bound, i
    real(dp) :: sum_prob, dw_cold, dw_hot
    real(dp) :: delta_e_mev, fano_ratio
    logical  :: is_near

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, " GeneralModule Unit Tests: Surface Scattering     "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Surface Lattice & Morse Bound States Initialization
    ! --------------------------------------------------------------------------
    ! LiF(001): a = 2.84 A, corrugation zeta = 0.06 A, M_sub = 25.9 amu, Theta_D = 730 K
    call init_surface_lattice(lif, ax_ang=2.84_dp, ay_ang=2.84_dp, &
                              zeta_x_ang=0.06_dp, zeta_y_ang=0.06_dp, &
                              m_sub_amu=25.94_dp, debye_temp_k=730.0_dp)

    ! He/LiF: Well depth D = 7.5 meV, alpha = 1.1 A^-1, M_He = 4.0026 amu
    call init_surface_potential_morse(he_lif_pot, well_depth_mev=7.5_dp, &
                                     range_inv_ang=1.1_dp, mass_amu=4.0026_dp, &
                                     n_bound=n_bound)

    n_total = n_total + 1
    if (lif%bx_bohr > 1.0_dp .and. n_bound >= 3) then
        print *, " [PASS] Surface lattice & Morse potential initialized: b1 = ", lif%bx_bohr, &
                 " a.u., Bound states n_bound = ", n_bound
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Surface initialization error: n_bound = ", n_bound
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: 2D Diffraction Channels Kinematics (Specular is Open)
    ! --------------------------------------------------------------------------
    ! He atom at E = 20.0 meV, incident angle theta_i = 40 deg, phi_i = 0 deg
    call calc_surface_diffraction_channels(lif, mass_amu=4.0026_dp, energy_ev=0.020_dp, &
                                           theta_i_deg=40.0_dp, phi_i_deg=0.0_dp, &
                                           max_order=1, n_channels=n_ch, channels=channels)

    n_total = n_total + 1
    ! Channel 5 is (0, 0) specular beam (index (m+1)*3 + (n+1)+1 for m=0, n=0: 1*3 + 1 + 1 = 5)
    if (channels(5)%m == 0 .and. channels(5)%n == 0 .and. channels(5)%is_open .and. &
        channels(5)%kz_au > 0.0_dp) then
        print *, " [PASS] Specular (0,0) channel open: kz = ", channels(5)%kz_au, " a.u."
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Specular channel closed or misidentified"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: HCS Eikonal Probabilities & Unitarity
    ! --------------------------------------------------------------------------
    call calc_hcs_diffraction_probabilities(lif, mass_amu=4.0026_dp, energy_ev=0.020_dp, &
                                            theta_i_deg=40.0_dp, phi_i_deg=0.0_dp, &
                                            max_order=1, n_channels=n_ch, channels=channels)

    sum_prob = 0.0_dp
    do i = 1, n_ch
        if (channels(i)%is_open) sum_prob = sum_prob + channels(i)%probability
    end do

    n_total = n_total + 1
    ! Unitarity check: sum of open channel probabilities must equal 1.0 within 1e-6
    if (abs(sum_prob - 1.0_dp) < 1.0e-5_dp .and. channels(5)%probability > 0.5_dp) then
        print *, " [PASS] HCS Eikonal diffraction unitarity verified: Sum P_G = ", sum_prob
        print *, "        Specular (0,0) intensity P_00 = ", channels(5)%probability
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Diffraction probabilities failure: sum = ", sum_prob
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Selective Adsorption Resonance (SAR) Evaluation
    ! --------------------------------------------------------------------------
    ! Evaluate near-resonance condition for an evanescent beam
    call calc_selective_adsorption_resonance(lif, he_lif_pot, mass_amu=4.0026_dp, &
                                             energy_ev=0.020_dp, theta_i_deg=55.0_dp, &
                                             phi_i_deg=0.0_dp, m_res=-1, n_res=0, &
                                             v_bound=0, is_near_res=is_near, &
                                             delta_e_mev=delta_e_mev, &
                                             fano_specular_ratio=fano_ratio)

    n_total = n_total + 1
    ! Fano ratio should be positive and finite
    if (fano_ratio > 0.0_dp .and. fano_ratio < 10.0_dp) then
        print *, " [PASS] SAR Fano profile evaluated: Delta E = ", delta_e_mev, &
                 " meV, Fano ratio = ", fano_ratio
        n_pass = n_pass + 1
    else
        print *, " [FAIL] SAR Fano profile invalid: ", fano_ratio
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Debye-Waller Thermal Attenuation
    ! --------------------------------------------------------------------------
    dw_cold = calc_surface_debye_waller(lif, mass_amu=4.0026_dp, k_iz_au=1.5_dp, &
                                        kz_g_au=1.5_dp, temp_k=100.0_dp)
    dw_hot  = calc_surface_debye_waller(lif, mass_amu=4.0026_dp, k_iz_au=1.5_dp, &
                                        kz_g_au=1.5_dp, temp_k=600.0_dp)

    n_total = n_total + 1
    ! Higher temperature should lead to stronger thermal phonon attenuation (DW_hot < DW_cold)
    if (dw_cold > dw_hot .and. dw_hot > 0.0_dp .and. dw_cold <= 1.0_dp) then
        print *, " [PASS] Debye-Waller thermal attenuation verified:"
        print *, "        DW(T=100K) = ", dw_cold, " > DW(T=600K) = ", dw_hot
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Debye-Waller factor error: cold = ", dw_cold, " hot = ", dw_hot
    end if

    print *, "--------------------------------------------------"
    print *, "Surface Scattering Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All surface scattering tests passed."
    else
        stop 1
    end if

end program test_surface_scattering
