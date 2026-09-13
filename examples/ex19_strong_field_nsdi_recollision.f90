! ==============================================================================
! GeneralModule Example 19: Strong-Field NSDI & Recollision Momentum Dynamics
!
! Features:
!  1. Corkum three-step electron recollision trajectories and 3.17 Up cutoff
!  2. Lotz electron-impact ionization cross section sigma(e, 2e)
!  3. 2D correlated parallel electron momentum distribution P(pz1, pz2)
!  4. Strong positive correlation in 1st/3rd quadrants (COLTRIMS signature)
!  5. Intensity-dependent nonsequential double ionization "knee structure"
!
! Standard: Fortran 2008
! ==============================================================================
program ex19_strong_field_nsdi_recollision
    use mod_constants, only: dp, PI, HALFPI, AU2EV
    use mod_strong_field_nsdi
    implicit none

    type(nsdi_laser_t)  :: laser
    type(nsdi_target_t) :: target
    real(dp) :: phi_0, phi_r, e_rec, ratio_up, corr_c
    real(dp) :: grid_p(31), dist_2d(31, 31)
    real(dp) :: intensities(15), y_nsdi(15), y_sdi(15)
    integer  :: u_out, i, j, stat
    logical  :: ok

    print *, "================================================================"
    print *, " Example 19: Strong-Field NSDI & Correlated Momentum Spectrum   "
    print *, "================================================================"

    ! 1. Initialize 800 nm, 2.5e14 W/cm^2 laser pulse and Argon target
    call init_nsdi_laser(800.0_dp, 2.5e14_dp, laser, stat)
    call init_nsdi_target("Ar", target, stat)

    write(*, '(A, F8.1, A)') " Laser wavelength              : ", laser%wavelength_nm, " nm"
    write(*, '(A, ES10.2, A)')" Peak intensity                : ", laser%intensity_w_cm2, " W/cm^2"
    write(*, '(A, F8.4, A)') " Ponderomotive energy Up       : ", laser%up_au * AU2EV, " eV"
    write(*, '(A, F8.4, A)') " Classical 3.17 Up cutoff      : ", 3.173_dp * laser%up_au * AU2EV, " eV"
    write(*, '(A, F8.4, A)') " Target Ar+ second IP (Ip2)    : ", target%ip2_au * AU2EV, " eV"
    print *, "----------------------------------------------------------------"

    ! 2. Scan recollision return energy vs birth phase phi_0
    open(newunit=u_out, file="ex19_nsdi_recollision.dat", status="replace", action="write")
    write(u_out, '(A)') "# GeneralModule Example 19: Recollision Trajectory Dynamics"
    write(u_out, '(A)') "# Col 1: Birth Phase phi_0 (deg)"
    write(u_out, '(A)') "# Col 2: Return Phase phi_r (deg)"
    write(u_out, '(A)') "# Col 3: Return Kinetic Energy E_rec / Up"
    write(u_out, '(A)') "# Col 4: Return Kinetic Energy E_rec (eV)"

    do i = 1, 40
        phi_0 = 0.05_dp + real(i - 1, dp) * (HALFPI - 0.10_dp) / 39.0_dp
        call calc_recollision_trajectory(phi_0, laser%up_au, phi_r, e_rec, ok)
        if (ok) then
            ratio_up = e_rec / laser%up_au
            write(u_out, '(4ES16.6)') phi_0 * 180.0_dp / PI, phi_r * 180.0_dp / PI, &
                                      ratio_up, e_rec * AU2EV
        end if
    end do
    close(u_out)

    ! 3. Generate 2D correlated parallel momentum spectrum P(pz1, pz2)
    call calc_nsdi_2d_momentum_dist(laser, target, 31, 2.5_dp, grid_p, dist_2d, corr_c)
    write(*, '(A, F8.4)') " Two-electron momentum correlation coefficient C_corr: ", corr_c

    open(newunit=u_out, file="ex19_nsdi_momentum_2d.dat", status="replace", action="write")
    write(u_out, '(A)') "# 2D Correlated Parallel Momentum Distribution P(pz1, pz2)"
    write(u_out, '(A)') "# Col 1: pz1 (a.u.)"
    write(u_out, '(A)') "# Col 2: pz2 (a.u.)"
    write(u_out, '(A)') "# Col 3: Probability Density P(pz1, pz2)"

    do i = 1, 31
        do j = 1, 31
            write(u_out, '(3ES16.6)') grid_p(i), grid_p(j), dist_2d(i, j)
        end do
        write(u_out, '(A)') ""  ! Blank line for gnuplot pm3d
    end do
    close(u_out)

    ! 4. Double ionization yield curve and knee structure
    call calc_double_ion_yield_curve(800.0_dp, target, 15, 1.5e14_dp, 8.0e14_dp, &
                                     intensities, y_nsdi, y_sdi)

    open(newunit=u_out, file="ex19_nsdi_knee.dat", status="replace", action="write")
    write(u_out, '(A)') "# Double Ionization Yield Curve (Nonsequential Knee Structure)"
    write(u_out, '(A)') "# Col 1: Intensity I (W/cm^2)"
    write(u_out, '(A)') "# Col 2: NSDI Yield"
    write(u_out, '(A)') "# Col 3: SDI (Sequential) Yield"
    write(u_out, '(A)') "# Col 4: Total Yield (NSDI + SDI)"

    write(*, '(A)') " Intensity (W/cm^2) |    NSDI Yield     |     SDI Yield     |  NSDI / SDI Ratio"
    write(*, '(A)') "--------------------+-------------------+-------------------+-------------------"

    do i = 1, 15
        write(u_out, '(4ES16.6)') intensities(i), y_nsdi(i), y_sdi(i), y_nsdi(i) + y_sdi(i)
        if (mod(i, 3) == 1) then
            write(*, '(ES19.2, " | ", ES17.4, " | ", ES17.4, " | ", ES17.4)') &
                intensities(i), y_nsdi(i), y_sdi(i), y_nsdi(i) / max(y_sdi(i), 1.0e-30_dp)
        end if
    end do
    close(u_out)

    print *, "----------------------------------------------------------------"
    print *, " Output data saved to ex19_nsdi_recollision.dat, "
    print *, " ex19_nsdi_momentum_2d.dat, and ex19_nsdi_knee.dat"
    print *, " Example 19 completed successfully."

end program ex19_strong_field_nsdi_recollision
