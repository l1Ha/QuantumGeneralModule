! ==============================================================================
! GeneralModule Example 11: Ultracold Three-Body Efimov Recombination
!
! Features:
!  1. Identical bosons (133Cs) Efimov transcendental root s_0 and scaling lambda
!  2. Three-body recombination loss rate K_3(a) for a > 0 (Interference minimum)
!  3. Efimov resonance peak for a < 0 (Universal Efimov trimer resonance)
!  4. Unitary limit finite-temperature loss saturation: K_3 ~ T^(-2)
!
! Standard: Fortran 2008
! ==============================================================================
program ex11_three_body_efimov_recombination
    use mod_constants, only: dp, PI, AMU2AU
    use mod_three_body_recombination
    implicit none

    real(dp) :: mass_cs, s0_val, lambda_val, k3_t_au, k3_t_si
    real(dp) :: a_val, temp_uk
    type(efimov_param_t)   :: param
    type(three_body_loss_t):: res_pos, res_neg
    integer :: u_out, i, n_pts, stat

    print *, "================================================================"
    print *, " Example 11: Ultracold Three-Body Efimov Recombination Rates   "
    print *, "================================================================"

    ! 1. Initialize 133Cs identical boson Efimov system
    mass_cs = 132.90545_dp * AMU2AU
    s0_val = solve_efimov_s0_identical_bosons()
    lambda_val = calc_efimov_scale_factor(s0_val)

    param%s0 = s0_val
    param%scale_factor = lambda_val
    param%a_star = 200.0_dp
    param%a_minus = -100.0_dp
    param%eta_star = 0.06_dp

    write(*, '(A, F10.6)') " Efimov universal root s_0          : ", s0_val
    write(*, '(A, F10.4)') " Discrete scaling factor lambda     : ", lambda_val
    write(*, '(A, F10.2, A)') " Three-body parameter a_* (a > 0)  : ", param%a_star, " a_0"
    write(*, '(A, F10.2, A)') " Resonant trimer pole a_- (a < 0)  : ", param%a_minus, " a_0"
    print *, "----------------------------------------------------------------"

    ! 2. Open output file for data plotting
    open(newunit=u_out, file="three_body_efimov_recombination.dat", status="replace", action="write")
    write(u_out, '(A)') "# Three-Body Efimov Recombination Loss Rate K_3(a)"
    write(u_out, '(A)') "# Col 1: |a| (Bohr a_0)"
    write(u_out, '(A)') "# Col 2: K_3(a > 0) [cm^6/s]"
    write(u_out, '(A)') "# Col 3: K_3(a < 0) [cm^6/s]"

    ! Compute K_3 across a range of scattering lengths
    n_pts = 60
    do i = 1, n_pts
        ! Logarithmic grid from 20 to 10000 a_0
        a_val = 20.0_dp * 10.0_dp**(real(i - 1, dp) * (log10(10000.0_dp / 20.0_dp) / real(n_pts - 1, dp)))

        call calc_three_body_recombination_a_positive(a_val, mass_cs, param, res_pos, stat)
        call calc_three_body_recombination_a_negative(-a_val, mass_cs, param, res_neg, stat)

        write(u_out, '(3ES18.8)') a_val, res_pos%k3_si, res_neg%k3_si
    end do
    close(u_out)
    print *, " Recombination rate profiles written to: three_body_efimov_recombination.dat"

    ! 3. Finite-Temperature Unitary Limit Power Law K_3 ~ T^(-2)
    print *, "----------------------------------------------------------------"
    print *, " Finite-Temperature Loss at Unitary Limit (|a| -> infinity):"
    write(*, '(A12, A24)') " T (uK)", " K_3(T) [cm^6/s]"
    do i = 1, 5
        temp_uk = 0.5_dp * real(i, dp)
        call calc_unitary_three_body_loss_temperature(temp_uk * 1.0e-6_dp, mass_cs, param%eta_star, &
                                                      k3_t_au, k3_t_si, stat)
        write(*, '(F12.2, ES24.6)') temp_uk, k3_t_si
    end do

    print *, "================================================================"
    print *, " Example 11 completed successfully.                             "
    print *, "================================================================"
end program ex11_three_body_efimov_recombination
