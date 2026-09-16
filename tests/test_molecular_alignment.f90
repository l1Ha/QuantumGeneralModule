!> \file test_molecular_alignment.f90
!> \brief Unit tests for laser-induced molecular alignment, orientation, and superrotors
!> \author LiHao
program test_molecular_alignment
    use mod_constants, only: dp
    use mod_molecular_alignment
    implicit none

    type(rotor_molecule_t) :: n2_mol
    real(dp), dimension(0:8)   :: cos2_diag
    real(dp), dimension(0:6)   :: cos2_off2
    real(dp), dimension(0:7)   :: cos_off1
    real(dp), dimension(60)    :: t_grid, cos2_trace
    integer  :: j_super, n_pass, n_total, stat
    logical  :: is_dissoc
    real(dp) :: e_rot_ev

    n_pass = 0
    n_total = 0

    print '(A)', "=================================================="
    print '(A)', "  GeneralModule Unit Tests: Molecular Alignment   "
    print '(A)', "=================================================="

    ! --------------------------------------------------------------------------
    ! Test 1: Linear Rotor Initialization & Revival Period (N2)
    ! --------------------------------------------------------------------------
    call init_rotor_molecule(n2_mol, "N2", b_rot_cm1=1.9982_dp, delta_alpha_au=6.70_dp, &
                             dipole_debye=0.0_dp, stat=stat)
    n_total = n_total + 1
    ! N2 revival period T_rev ~ 8.347 ps
    if (stat == 0 .and. abs(n2_mol%t_rev_ps - 8.3466_dp) < 0.01_dp) then
        print '(A, F7.3, A)', " [PASS] N2 molecule initialized: T_rev = ", n2_mol%t_rev_ps, " ps"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Molecular rotor initialization error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 2: Analytical Spherical Harmonic Matrix Elements <J' M| cos^2 theta |J M>
    ! --------------------------------------------------------------------------
    call calc_cos2_matrix_elements(j_max=8, m_proj=0, cos2_diag=cos2_diag, cos2_off2=cos2_off2)
    call calc_cos_matrix_elements(j_max=8, m_proj=0, cos_off1=cos_off1)

    n_total = n_total + 1
    ! <0,0|cos^2 theta|0,0> = 1/3, <2,0|cos^2 theta|0,0> = sqrt(4/45) ~ 0.298142
    ! <1,0|cos theta|0,0> = 1/sqrt(3) ~ 0.577350
    if (abs(cos2_diag(0) - 1.0_dp/3.0_dp) < 1.0e-12_dp .and. &
        abs(cos2_off2(0) - sqrt(4.0_dp / 45.0_dp)) < 1.0e-12_dp .and. &
        abs(cos_off1(0) - 1.0_dp / sqrt(3.0_dp)) < 1.0e-12_dp) then
        print '(A, F8.5, A, F8.5, A)', " [PASS] Matrix elements verified: <0|cos^2|0> = ", &
            cos2_diag(0), ", <2|cos^2|0> = ", cos2_off2(0)
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Analytical angular matrix elements incorrect"
    end if

    ! --------------------------------------------------------------------------
    ! Test 3: Laser-Induced Non-Adiabatic Alignment Dynamics
    ! --------------------------------------------------------------------------
    ! Short 100 fs pulse on N2 at T = 0 K
    call simulate_laser_induced_alignment(n2_mol, laser_i0_wcm2=2.0e13_dp, pulse_fwhm_fs=100.0_dp, &
                                         temp_k=0.0_dp, j_max=12, n_time=60, t_span_ps=3.0_dp, &
                                         t_grid_ps=t_grid, cos2_trace=cos2_trace, stat=stat)
    n_total = n_total + 1
    ! Post-pulse alignment should peak above isotropic 0.333
    if (stat == 0 .and. maxval(cos2_trace) > 0.40_dp .and. cos2_trace(1) > 0.32_dp) then
        print '(A, F8.4, A, F8.4)', " [PASS] Laser alignment transient: Peak <cos^2> = ", &
            maxval(cos2_trace), ", Initial = ", cos2_trace(1)
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Laser alignment dynamic simulation error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 4: Optical Centrifuge Accelerated Superrotor J State
    ! --------------------------------------------------------------------------
    ! Chirp acceleration 0.4 THz/ps over 60 ps
    call calc_optical_centrifuge_kick(n2_mol, chirp_beta_thz2=0.4_dp, pulse_duration_ps=60.0_dp, &
                                     j_superrotor=j_super, stat=stat)
    n_total = n_total + 1
    ! For N2, omega = 24 THz -> J ~ 32
    if (stat == 0 .and. j_super > 20 .and. j_super < 50) then
        print '(A, I3, A)', " [PASS] Optical centrifuge accelerated to extreme superrotor state J = ", &
            j_super, " hbar"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Optical centrifuge state kick error"
    end if

    ! --------------------------------------------------------------------------
    ! Test 5: Superrotor Centrifugal Dissociation Criterion
    ! --------------------------------------------------------------------------
    ! High J state J = 180 on a typical molecule with D_e = 4.0 eV
    call calc_superrotor_dissociation(n2_mol, j_rot=200, d_e_ev=4.0_dp, &
                                      is_dissociated=is_dissoc, rot_energy_ev=e_rot_ev)
    n_total = n_total + 1
    if (is_dissoc .and. e_rot_ev > 4.0_dp) then
        print '(A, F8.3, A)', " [PASS] Centrifugal bond-breaking verified: E_rot(J=200) = ", &
            e_rot_ev, " eV > D_e"
        n_pass = n_pass + 1
    else
        print '(A)', " [FAIL] Centrifugal dissociation threshold evaluation error"
    end if

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', " Molecular Alignment Tests: ", n_pass, " / ", n_total, " PASSED."
    if (n_pass == n_total) then
        print '(A)', " SUCCESS: All molecular alignment & superrotor tests passed."
    else
        error stop " Test failures detected in test_molecular_alignment."
    end if

end program test_molecular_alignment
