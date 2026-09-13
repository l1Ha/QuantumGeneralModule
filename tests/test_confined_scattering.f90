! ==============================================================================
! GeneralModule: test_confined_scattering.f90
!
! Unit Test Suite: Ultracold Confined Quantum Scattering & CIR
!
! Standard: Fortran 2008
! ==============================================================================

program test_confined_scattering
    use mod_constants, only: dp, PI, TWOPI, AMU2AU, HBAR
    use mod_confined_scattering
    implicit none

    integer :: n_pass, n_total
    integer :: stat
    type(waveguide_1d_t) :: wg
    type(planar_2d_t)    :: pl
    type(cir_result_t)   :: res_weak, res_cir, res_tg
    real(dp) :: mu_rb87, omega_perp, a_s_weak, a_s_cir
    real(dp) :: e_b, a_cir_x, a_cir_y
    real(dp) :: a_2d, sigma_2d, k_wave
    complex(dp) :: f_2d
    real(dp) :: g_1d_expected

    n_pass = 0
    n_total = 0

    print *, "=================================================="
    print *, "  GeneralModule Unit Tests: Confined & CIR        "
    print *, "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Waveguide Initialization & Oscillator Length
    ! --------------------------------------------------------------------------
    ! 87Rb + 87Rb reduced mass mu = 43.5 amu
    mu_rb87 = 43.5_dp * AMU2AU
    ! Transverse trap frequency: omega_perp = 2 * pi * 50 kHz ~ 7.61e-12 a.u.
    omega_perp = TWOPI * 50.0e3_dp * (2.4188843265857e-17_dp)

    call init_waveguide_1d(omega_perp, mu_rb87, wg, stat)

    n_total = n_total + 1
    if (stat == 0) then
        print *, " [PASS] 1D Waveguide initialized successfully"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Waveguide initialization failed"
    end if

    n_total = n_total + 1
    if (wg%a_perp_au > 100.0_dp) then
        print *, " [PASS] Transverse oscillator length a_perp = ", wg%a_perp_au, " a0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Unexpected a_perp length: ", wg%a_perp_au
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Olshanii CIR Resonance Pole Position a_CIR = a_perp / C
    ! --------------------------------------------------------------------------
    n_total = n_total + 1
    if (abs(wg%a_cir_au / wg%a_perp_au - 0.968428_dp) < 1.0e-4_dp) then
        print *, " [PASS] CIR resonance pole ratio a_CIR / a_perp = ", &
                 wg%a_cir_au / wg%a_perp_au, " (Exact: 0.96843)"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] CIR pole ratio error: ", wg%a_cir_au / wg%a_perp_au
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Weak Interaction Limit (a_s << a_perp)
    ! --------------------------------------------------------------------------
    a_s_weak = 10.0_dp  ! a_s = 10 a_0 << a_perp (~ 1000 a0)
    call calc_olshanii_cir_parameters(wg, a_s_weak, 1.0e-4_dp, res_weak, stat)

    g_1d_expected = (2.0_dp * a_s_weak) / (mu_rb87 * (wg%a_perp_au**2))

    n_total = n_total + 1
    if (stat == 0 .and. abs(res_weak%g_1d_au / g_1d_expected - 1.0_dp) < 0.02_dp) then
        print *, " [PASS] Weak interaction matches 2*a_s / (mu*a_perp^2) within 2%"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Weak interaction limit mismatch: ", &
                 res_weak%g_1d_au, " vs ", g_1d_expected
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Near CIR Resonance Divergence & Flaggings
    ! --------------------------------------------------------------------------
    a_s_cir = wg%a_cir_au  ! Exactly on CIR resonance
    call calc_olshanii_cir_parameters(wg, a_s_cir, 1.0e-4_dp, res_cir, stat)

    n_total = n_total + 1
    if (res_cir%is_near_cir .and. res_cir%is_tonks_regime) then
        print *, " [PASS] CIR resonance correctly identified and enters Tonks-Girardeau"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] CIR resonance detection failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Confinement-Induced Dimer Binding Energy E_b
    ! --------------------------------------------------------------------------
    e_b = calc_confined_dimer_binding_energy(wg, 200.0_dp)

    n_total = n_total + 1
    if (e_b < 0.0_dp) then
        print *, " [PASS] Confinement-induced dimer bound state E_b = ", e_b, " a.u. < 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Bound state energy is non-negative: ", e_b
    end if

    ! --------------------------------------------------------------------------
    ! Test 6: Tonks-Girardeau Parameter Scaling with 1D Density
    ! --------------------------------------------------------------------------
    ! At very dilute density, gamma_LL should be large (strongly interacting)
    call calc_olshanii_cir_parameters(wg, 100.0_dp, 1.0e-7_dp, res_tg, stat)

    n_total = n_total + 1
    if (res_tg%gamma_ll > 10.0_dp .and. res_tg%is_tonks_regime) then
        print *, " [PASS] Dilute 1D gas enters Tonks-Girardeau regime: gamma_LL = ", &
                 res_tg%gamma_ll, " > 10"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Tonks-Girardeau regime criteria failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 7: Anisotropic Waveguide CIR Doublet Splitting
    ! --------------------------------------------------------------------------
    ! Anisotropic frequencies: omega_x = 60 kHz, omega_y = 40 kHz (eta = 1.5)
    call calc_anisotropic_cir_poles(omega_perp * 1.2_dp, omega_perp * 0.8_dp, &
                                    mu_rb87, a_cir_x, a_cir_y, stat)

    n_total = n_total + 1
    if (stat == 0 .and. a_cir_x /= a_cir_y) then
        print *, " [PASS] Anisotropic CIR poles split: a_CIR,x = ", a_cir_x, &
                 ", a_CIR,y = ", a_cir_y
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Anisotropic CIR pole splitting failed"
    end if

    n_total = n_total + 1
    if (a_cir_x < a_cir_y) then
        print *, " [PASS] Stiffer axis x has smaller CIR pole: ", a_cir_x, " < ", a_cir_y
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Anisotropic CIR ordering failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 8: Quasi-2D Planar Trap Initialization
    ! --------------------------------------------------------------------------
    call init_planar_2d(omega_perp * 2.0_dp, mu_rb87, pl, stat)

    n_total = n_total + 1
    if (stat == 0 .and. pl%a_z_au > 0.0_dp) then
        print *, " [PASS] Quasi-2D planar trap initialized: a_z = ", pl%a_z_au, " a0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Quasi-2D initialization failed"
    end if

    ! --------------------------------------------------------------------------
    ! Test 9: Quasi-2D Scattering Length and Cross Section
    ! --------------------------------------------------------------------------
    k_wave = 1.0e-3_dp  ! Low energy k = 1e-3 a.u.
    call calc_quasi_2d_scattering(pl, a_s_3d_au=100.0_dp, k_wave_au=k_wave, &
                                  a_2d_au=a_2d, f_2d_amp=f_2d, &
                                  sigma_2d_au=sigma_2d, stat=stat)

    n_total = n_total + 1
    if (stat == 0 .and. a_2d > 0.0_dp) then
        print *, " [PASS] Quasi-2D scattering length a_2D = ", a_2d, " a0 > 0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Quasi-2D scattering length failed"
    end if

    n_total = n_total + 1
    if (sigma_2d > 0.0_dp .and. abs(f_2d) > 0.0_dp) then
        print *, " [PASS] Quasi-2D total cross section sigma_2D = ", sigma_2d, " a0"
        n_pass = n_pass + 1
    else
        print *, " [FAIL] Quasi-2D cross section non-positive"
    end if

    print *, "--------------------------------------------------"
    print *, "Confined & CIR Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print *, "SUCCESS: All confined quantum scattering tests passed."
    else
        stop 1
    end if

end program test_confined_scattering
