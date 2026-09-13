! ==============================================================================
! GeneralModule Example 09: 各向异性磁偶极与电偶极量子散射及自旋弛豫模拟
! 理论参考文献:
!   - H. T. C. Stoof et al., Phys. Rev. B 38, 4688 (1988)
!   - A. J. Moerdijk & B. J. Verhaar, Phys. Rev. A 53, 4341 (1996)
!   - J. L. Bohn, M. Cavagnero, C. Ticknor, New J. Phys. 11, 055039 (2009)
!   - K.-K. Ni et al., Science 322, 231 (2008)
! 演示内容：
! 1. 87Rb 磁阱中弱场寻优态原子的非弹性磁偶极自旋弛豫 (Dipolar Relaxation)
! 2. s-d 分波各向异性耦合强度与外加静磁场 B (1 ~ 50 Gauss) 扫描
! 3. 极性分子 KRb 外加直流电场 (0 ~ 20 kV/cm) 下 Stark 诱导偶极与特征偶极力程 a_d
! ==============================================================================
program ex09_dipolar_relaxation_scattering
    use mod_constants, only: dp, PI, TO_AU, FROM_AU
    use mod_dipolar_scattering
    use mod_io_utils, only: print_banner
    implicit none

    integer, parameter :: N_B = 21, N_E_FIELD = 21
    real(dp) :: b_scan(N_B), sigma_rel(N_B), k_rel(N_B)
    real(dp) :: e_field_scan(N_E_FIELD), d_ind_scan(N_E_FIELD), ad_scan(N_E_FIELD)
    type(dipolar_relaxation_result_t) :: res
    type(polar_molecule_t) :: krb
    real(dp) :: e_coll_au, b_gauss, ef_val, d_val, ad_val
    integer  :: i, u, st

    call print_banner("GeneralModule Example 09: Dipolar Quantum Scattering", 70)

    ! --------------------------------------------------------------------------
    ! 1. 87Rb 原子磁偶极自旋弛豫 (s-波入，d-波出) 随外磁场扫描
    ! --------------------------------------------------------------------------
    print *, ">> [1/2] Simulating 87Rb Inelastic Dipolar Relaxation in Magnetic Field..."
    print *, "   Collision Energy: 1.0 uK (~ 3.17e-11 a.u.)"
    print *, "   Initial Channel : |S=1, Ms=1, L=0, Ml=0> (Weak-field seeking)"
    print *, "   Exit Channels   : |S=1, Ms=0, L=2, Ml=1> & |S=1, Ms=-1, L=2, Ml=2>"
    print *, "----------------------------------------------------------------------"
    print '(A12, A18, A20, A20)', "B (Gauss)", "Delta E (uK)", "Cross Sec (a0^2)", "K_rel (cm^3/s)"
    print *, "----------------------------------------------------------------------"

    e_coll_au = 1.0e-6_dp * 3.166811563e-6_dp ! 1 uK in a.u.

    open(newunit=u, file="dipolar_relaxation_scan.dat", status="replace", action="write")
    write(u, '(A)') "# 87Rb Dipolar Relaxation vs Magnetic Field B"
    write(u, '(A)') "# B(Gauss)  Delta_E(uK)  sigma_rel(a0^2)  K_rel(cm^3/s)"

    do i = 1, N_B
        b_gauss = 1.0_dp + real(i - 1, dp) * (49.0_dp / real(N_B - 1, dp))
        call calc_dipolar_relaxation_cross_section(86.909_dp, b_gauss, e_coll_au, 8.0_dp, res, st)

        b_scan(i)    = b_gauss
        sigma_rel(i) = res%cross_section_au
        k_rel(i)     = res%rate_coeff_cm3_s

        ! 将释放的能量转换为 uK
        write(u, '(4ES16.7)') b_gauss, res%e_released / 3.166811563e-12_dp, &
                              res%cross_section_au, res%rate_coeff_cm3_s

        if (mod(i, 4) == 1 .or. i == N_B) then
            print '(F10.2, 3ES20.6)', b_gauss, res%e_released / 3.166811563e-12_dp, &
                                      res%cross_section_au, res%rate_coeff_cm3_s
        end if
    end do
    close(u)
    print *, ">> 87Rb dipolar relaxation scan saved to: dipolar_relaxation_scan.dat"

    ! --------------------------------------------------------------------------
    ! 2. 极性分子 KRb 在直流电场下的诱导偶极矩与特征偶极长度 a_d 扫描
    ! --------------------------------------------------------------------------
    print *
    print *, ">> [2/2] Simulating Polar Molecule (KRb) Stark Induced Dipole & Dipole Length a_d..."
    krb%name         = "KRb"
    krb%mass_amu     = 126.37_dp
    krb%b_rot_cm1    = 0.037_dp
    krb%dipole_debye = 0.566_dp
    print '(A, F8.4, A, F8.4, A)', "   KRb Rotational Constant B_rot: ", krb%b_rot_cm1, &
                                   " cm^-1, Permanent Dipole: ", krb%dipole_debye, " Debye"
    print *, "----------------------------------------------------------------------"
    print '(A15, A20, A20, A15)', "E_field (kV/cm)", "Induced Dipole (D)", "Ratio d_ind/d0", "a_d (a0)"
    print *, "----------------------------------------------------------------------"

    open(newunit=u, file="polar_molecule_stark_scan.dat", status="replace", action="write")
    write(u, '(A)') "# KRb Polar Molecule Stark Induced Dipole vs Electric Field"
    write(u, '(A)') "# E_field(kV/cm)  d_ind(Debye)  d_ind/d0  a_d(a0)"

    do i = 1, N_E_FIELD
        ef_val = real(i - 1, dp) * (20.0_dp / real(N_E_FIELD - 1, dp))
        call calc_stark_induced_dipole(krb, ef_val, d_val, st)
        ad_val = calc_dipole_length_scale(krb%mass_amu, d_val)

        e_field_scan(i) = ef_val
        d_ind_scan(i)   = d_val
        ad_scan(i)      = ad_val

        write(u, '(4ES16.7)') ef_val, d_val, d_val / krb%dipole_debye, ad_val

        if (mod(i, 4) == 1 .or. i == N_E_FIELD) then
            print '(F12.2, 3ES20.6)', ef_val, d_val, d_val / krb%dipole_debye, ad_val
        end if
    end do
    close(u)
    print *, ">> KRb Stark scan data saved to: polar_molecule_stark_scan.dat"

    print *
    print *, "======================================================================"
    print *, "   SUCCESS: Example 09 (Dipolar Scattering) Completed Successfully!   "
    print *, "======================================================================"

end program ex09_dipolar_relaxation_scattering
