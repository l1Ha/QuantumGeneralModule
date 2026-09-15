! ==============================================================================
! GeneralModule Example 22: Bicircular Laser Fields and Chiral PECD
!
! Features:
!  1. Bicircular counter-rotating omega - 2*omega field synthesis (C3 trefoil)
!  2. Chiral tetrahedral model molecules: R- and S-enantiomers
!  3. Chirality measure pseudoscalar invariant chi inversion
!  4. Energy-resolved Photoelectron Circular Dichroism (PECD) beta1(E) & G_PECD(E)
!  5. 2D Photoelectron Angular Distribution (PAD) generation
!
! Standard: Fortran 2008
! ==============================================================================
program ex22_bicircular_pecd_chiral
    use mod_constants, only: dp, PI, FS2AU, AU2FS
    use mod_bicircular_pecd
    implicit none

    type(bicircular_field_t) :: field
    type(chiral_tetrahedral_molecule_t) :: mol_r, mol_s
    integer :: u_out, u_traj, i, ie
    integer, parameter :: N_PTS = 200, N_E = 30
    real(dp) :: t_arr(N_PTS), ex_arr(N_PTS), ey_arr(N_PTS), ax_arr(N_PTS), ay_arr(N_PTS)
    real(dp) :: e_grid(N_E), beta1_r(N_E), g_r(N_E), beta1_s(N_E), g_s(N_E)
    real(dp) :: chi_r, chi_s
    integer  :: n_fold

    print *, "================================================================"
    print *, " Example 22: Bicircular Fields and Chiral Molecule PECD         "
    print *, "================================================================"

    ! 1. Synthesize counter-rotating omega + 2*omega bicircular laser field
    call init_bicircular_field(field, omega1_au=0.057_dp, r_freq=2.0_dp, &
                               i1_wcm2=1.0e14_dp, i2_wcm2=5.0e13_dp, &
                               h1=1, h2=-1, phi1=0.0_dp, phi2=0.0_dp, &
                               fwhm_fs=25.0_dp, envelope_type=1)

    n_fold = calc_dynamical_symmetry_fold(field%h1, field%h2, freq_ratio=2)
    write(*, '(A, I2, A)') " Field Dynamical Symmetry Fold: C", n_fold, " (Trefoil / 3-fold)"

    ! Generate Lissajous parametric trajectory
    call calc_bicircular_trajectory(field, N_PTS, 20.0_dp * FS2AU, &
                                    t_arr, ex_arr, ey_arr, ax_arr, ay_arr)

    open(newunit=u_traj, file="ex22_bicircular_field.dat", status="replace", action="write")
    write(u_traj, '(A)') "# GeneralModule Example 22: Bicircular Field Trajectory"
    write(u_traj, '(A)') "# Col 1: Time t (fs)"
    write(u_traj, '(A)') "# Col 2: Ex (a.u.)"
    write(u_traj, '(A)') "# Col 3: Ey (a.u.)"
    write(u_traj, '(A)') "# Col 4: Ax (a.u.)"
    write(u_traj, '(A)') "# Col 5: Ay (a.u.)"

    do i = 1, N_PTS
        write(u_traj, '(5ES16.6)') t_arr(i) * AU2FS, ex_arr(i), ey_arr(i), ax_arr(i), ay_arr(i)
    end do
    close(u_traj)
    print *, " Bicircular field trajectory saved to ex22_bicircular_field.dat"
    print *, "----------------------------------------------------------------"

    ! 2. Initialize R and S enantiomer tetrahedral molecules
    call init_chiral_tetrahedral_molecule(mol_r, "R")
    call init_chiral_tetrahedral_molecule(mol_s, "S")

    chi_r = calc_chirality_measure(mol_r)
    chi_s = calc_chirality_measure(mol_s)
    write(*, '(A, F10.3)') " Chirality Measure chi(R) : ", chi_r
    write(*, '(A, F10.3)') " Chirality Measure chi(S) : ", chi_s
    print *, "----------------------------------------------------------------"

    ! 3. Energy-Resolved PECD and Forward-Backward Asymmetry
    call calc_pecd_energy_resolved(mol_r, 1.0_dp, 20.0_dp, N_E, e_grid, beta1_r, g_r)
    call calc_pecd_energy_resolved(mol_s, 1.0_dp, 20.0_dp, N_E, e_grid, beta1_s, g_s)

    open(newunit=u_out, file="ex22_chiral_pecd.dat", status="replace", action="write")
    write(u_out, '(A)') "# GeneralModule Example 22: Energy-Resolved PECD Spectrum"
    write(u_out, '(A)') "# Col 1: Photoelectron Energy (eV)"
    write(u_out, '(A)') "# Col 2: beta1 (R-enantiomer)"
    write(u_out, '(A)') "# Col 3: G_PECD (R-enantiomer, %)"
    write(u_out, '(A)') "# Col 4: beta1 (S-enantiomer)"
    write(u_out, '(A)') "# Col 5: G_PECD (S-enantiomer, %)"

    write(*, '(A)') " Energy (eV) |   beta1(R)   | G_PECD(R) (%) |   beta1(S)   | G_PECD(S) (%)"
    write(*, '(A)') "-------------+--------------+---------------+--------------+--------------"

    do ie = 1, N_E
        write(u_out, '(5ES16.6)') e_grid(ie), beta1_r(ie), g_r(ie) * 100.0_dp, &
                                  beta1_s(ie), g_s(ie) * 100.0_dp
        if (mod(ie, 5) == 1) then
            write(*, '(F12.2, " | ", ES12.4, " | ", F13.3, " | ", ES12.4, " | ", F13.3)') &
                e_grid(ie), beta1_r(ie), g_r(ie) * 100.0_dp, beta1_s(ie), g_s(ie) * 100.0_dp
        end if
    end do
    close(u_out)

    print *, "----------------------------------------------------------------"
    print *, " SUCCESS: PECD spectra saved to ex22_chiral_pecd.dat"
    print *, "================================================================"

end program ex22_bicircular_pecd_chiral
