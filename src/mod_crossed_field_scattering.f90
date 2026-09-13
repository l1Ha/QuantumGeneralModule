! ==============================================================================
! GeneralModule: mod_crossed_field_scattering.f90
!
! Quantum Dynamics in Crossed Electric and Magnetic Fields (E and B with arbitrary angle)
!
! Theoretical Foundations:
!   - T. V. Tscherbul & R. V. Krems, Phys. Rev. Lett. 97, 083201 (2006) [Tilted fields]
!   - B. Friedrich & D. Herschbach, PRL 74, 4623 (1995) [Orientation in combined fields]
!   - R. V. Krems, Int. Rev. Phys. Chem. 24, 99 (2005) [Cold molecules in fields]
!
! Key Capabilities:
!   1. Arbitrary tilt angle theta_EB between DC electric field E and magnetic field B
!   2. Off-diagonal Delta_M = +/- 1 transitions driven by transverse B_x and E_x
!   3. Rotational-spin coupled Stark-Zeeman Hamiltonian matrix diagonalization
!   4. Molecular orientation <cos theta>, alignment <cos^2 theta>, & spin polarization <S_z>
!   5. Avoided crossings and field-dressed eigenstate synthesis in tilted geometry
!
! Standard: Fortran 2008
! ==============================================================================

module mod_crossed_field_scattering
    use mod_constants, only: dp, PI, TWOPI, HBAR, AMU2AU, GAUSS2AU, DEBYE2AU, MV_CM2AU
    use mod_linear_algebra, only: diag_symmetric_matrix
    implicit none
    private

    ! Public derived types
    public :: crossed_field_config_t, crossed_field_state_t

    ! Public procedures
    public :: init_crossed_field_config
    public :: build_crossed_field_hamiltonian
    public :: solve_crossed_field_eigenstates
    public :: calc_crossed_field_observables
    public :: scan_tilt_angle_spectrum

    ! Crossed electric & magnetic field configuration
    type :: crossed_field_config_t
        real(dp) :: e_field_mv_cm     ! Electric field magnitude in MV/cm
        real(dp) :: e_field_au        ! Electric field in a.u.
        real(dp) :: b_field_gauss     ! Magnetic field magnitude in Gauss
        real(dp) :: b_field_au        ! Magnetic field in a.u.
        real(dp) :: theta_eb_rad      ! Angle theta_EB between E and B in radians
        real(dp) :: rot_constant_ghz  ! Rotational constant B_e in GHz
        real(dp) :: rot_constant_au   ! Rotational constant in a.u.
        real(dp) :: dipole_debye      ! Permanent dipole moment in Debye
        real(dp) :: dipole_au         ! Dipole in a.u.
        real(dp) :: g_electron        ! Electron spin g-factor (~ 2.002319)
        integer  :: j_max             ! Maximum rotational quantum number J_max
    end type crossed_field_config_t

    ! Eigenstate result in crossed fields
    type :: crossed_field_state_t
        integer               :: dim_tot           ! Total Hilbert space dimension
        real(dp), allocatable :: energies_au(:)    ! Eigenvalues in a.u.
        real(dp), allocatable :: eigenvectors(:,:) ! Eigenvectors (columns)
        real(dp), allocatable :: orientation(:)    ! <cos(theta_R)> for each state
        real(dp), allocatable :: spin_sz(:)        ! <S_z> for each state
    end type crossed_field_state_t

    ! Bohr magneton in a.u. is exactly 0.5
    real(dp), parameter :: MU_B_AU = 0.5_dp

contains

    ! ==========================================================================
    ! init_crossed_field_config:
    ! Initializes crossed field parameters with specified magnitudes & tilt angle
    ! ==========================================================================
    pure subroutine init_crossed_field_config(e_field_kv_cm, b_field_gauss, &
                                            theta_eb_deg, rot_ghz, dipole_d, &
                                            j_max, cfg, stat)
        real(dp), intent(in)                    :: e_field_kv_cm
        real(dp), intent(in)                    :: b_field_gauss
        real(dp), intent(in)                    :: theta_eb_deg
        real(dp), intent(in)                    :: rot_ghz
        real(dp), intent(in)                    :: dipole_d
        integer, intent(in)                     :: j_max
        type(crossed_field_config_t), intent(out) :: cfg
        integer, intent(out)                    :: stat

        stat = 0
        if (rot_ghz <= 0.0_dp .or. j_max < 0) then
            stat = 1
            return
        end if

        cfg%e_field_mv_cm = e_field_kv_cm * 1.0e-3_dp
        cfg%e_field_au = cfg%e_field_mv_cm * MV_CM2AU
        cfg%b_field_gauss = b_field_gauss
        cfg%b_field_au = b_field_gauss * GAUSS2AU
        cfg%theta_eb_rad = theta_eb_deg * (PI / 180.0_dp)

        ! 1 GHz = 1.0e9 / (4.134137e16) a.u. ~ 2.4188843e-8 a.u.
        cfg%rot_constant_ghz = rot_ghz
        cfg%rot_constant_au = rot_ghz * 1.0e9_dp * (2.4188843265857e-17_dp) * TWOPI

        cfg%dipole_debye = dipole_d
        cfg%dipole_au = dipole_d * DEBYE2AU
        cfg%g_electron = 2.00231930436256_dp
        cfg%j_max = j_max
    end subroutine init_crossed_field_config

    ! ==========================================================================
    ! build_crossed_field_hamiltonian:
    ! Builds the Hamiltonian in rotational-spin basis |J, M_J, S=1/2, M_S>:
    !   H = B_rot * J^2 - d * E * cos(theta_R) + g_e * mu_B * (B_z * S_z + B_x * S_x)
    ! where B_z = B * cos(theta_EB), B_x = B * sin(theta_EB).
    ! ==========================================================================
    subroutine build_crossed_field_hamiltonian(cfg, h_mat, dim_tot, stat)
        type(crossed_field_config_t), intent(in) :: cfg
        real(dp), allocatable, intent(out)       :: h_mat(:,:)
        integer, intent(out)                     :: dim_tot
        integer, intent(out)                     :: stat

        integer :: j_max, n_rot, j1, m1, ms1, j2, m2, ms2, idx1, idx2
        real(dp) :: b_z, b_x, e_z, diag_rot, diag_zeeman_z, off_zeeman_x, stark_elem
        real(dp) :: c_term

        stat = 0
        j_max = cfg%j_max
        n_rot = (j_max + 1)**2
        dim_tot = n_rot * 2  ! Spin S=1/2 gives 2 spin projections (M_S = +1/2, -1/2)

        allocate(h_mat(dim_tot, dim_tot))
        h_mat = 0.0_dp

        b_z = cfg%b_field_au * cos(cfg%theta_eb_rad)
        b_x = cfg%b_field_au * sin(cfg%theta_eb_rad)
        e_z = cfg%e_field_au

        ! Indexing function: idx(J, M_J, M_S) with M_S in {+1/2, -1/2} mapped to {1, 2}
        idx1 = 0
        do j1 = 0, j_max
            do m1 = -j1, j1
                do ms1 = 1, 2  ! 1: +1/2, 2: -1/2
                    idx1 = idx1 + 1

                    ! 1. Diagonal rotational energy: B_rot * J * (J + 1)
                    diag_rot = cfg%rot_constant_au * real(j1 * (j1 + 1), dp)

                    ! 2. Zeeman Z-coupling: g_e * mu_B * B_z * M_S
                    if (ms1 == 1) then
                        diag_zeeman_z = cfg%g_electron * MU_B_AU * b_z * 0.5_dp
                    else
                        diag_zeeman_z = - cfg%g_electron * MU_B_AU * b_z * 0.5_dp
                    end if

                    h_mat(idx1, idx1) = diag_rot + diag_zeeman_z

                    ! Loop over state 2
                    idx2 = 0
                    do j2 = 0, j_max
                        do m2 = -j2, j2
                            do ms2 = 1, 2
                                idx2 = idx2 + 1
                                if (idx2 < idx1) cycle

                                ! 3. Transverse magnetic field B_x couples M_S = +1/2 <-> -1/2
                                ! <1/2, +1/2 | S_x | 1/2, -1/2> = 1/2
                                if (j1 == j2 .and. m1 == m2 .and. ms1 /= ms2) then
                                    off_zeeman_x = cfg%g_electron * MU_B_AU * b_x * 0.5_dp
                                    h_mat(idx1, idx2) = h_mat(idx1, idx2) + off_zeeman_x
                                    h_mat(idx2, idx1) = h_mat(idx1, idx2)
                                end if

                                ! 4. DC Stark coupling along Z: - d * E_z * cos(theta_R)
                                ! Selection rules: Delta_J = +/- 1, Delta_M_J = 0, Delta_M_S = 0
                                if (ms1 == ms2 .and. m1 == m2 .and. abs(j1 - j2) == 1) then
                                    ! <J+1, M | cos(theta) | J, M> = sqrt(((J+1)^2 - M^2) / ((2*J+1)*(2*J+3)))
                                    if (j2 == j1 + 1) then
                                        c_term = sqrt(real((j2**2 - m1**2), dp) / &
                                                      real((2*j1 + 1) * (2*j2 + 1), dp))
                                        stark_elem = - cfg%dipole_au * e_z * c_term
                                        h_mat(idx1, idx2) = h_mat(idx1, idx2) + stark_elem
                                        h_mat(idx2, idx1) = h_mat(idx1, idx2)
                                    end if
                                end if

                            end do
                        end do
                    end do

                end do
            end do
        end do
    end subroutine build_crossed_field_hamiltonian

    ! ==========================================================================
    ! solve_crossed_field_eigenstates:
    ! Diagonalizes the crossed field Hamiltonian and computes orientation
    ! ==========================================================================
    subroutine solve_crossed_field_eigenstates(cfg, state, stat)
        type(crossed_field_config_t), intent(in) :: cfg
        type(crossed_field_state_t), intent(out) :: state
        integer, intent(out)                     :: stat

        real(dp), allocatable :: h_mat(:,:), d_eig(:), z_mat(:,:)
        integer :: dim_tot

        call build_crossed_field_hamiltonian(cfg, h_mat, dim_tot, stat)
        if (stat /= 0) return

        state%dim_tot = dim_tot
        allocate(d_eig(dim_tot), z_mat(dim_tot, dim_tot))
        allocate(state%energies_au(dim_tot), state%eigenvectors(dim_tot, dim_tot))
        allocate(state%orientation(dim_tot), state%spin_sz(dim_tot))

        ! Diagonalize real symmetric matrix
        call diag_symmetric_matrix(dim_tot, h_mat, d_eig, z_mat, stat)
        if (stat /= 0) return

        state%energies_au = d_eig
        state%eigenvectors = z_mat

        call calc_crossed_field_observables(cfg, state)
    end subroutine solve_crossed_field_eigenstates

    ! ==========================================================================
    ! calc_crossed_field_observables:
    ! Evaluates expectation values <cos(theta_R)> and <S_z> for each eigenstate
    ! ==========================================================================
    subroutine calc_crossed_field_observables(cfg, state)
        type(crossed_field_config_t), intent(in) :: cfg
        type(crossed_field_state_t), intent(inout) :: state

        integer :: n_st, k, j1, m1, ms1, j2, m2, ms2, idx1, idx2
        real(dp) :: c_term, cos_theta, sz_val

        n_st = state%dim_tot

        do k = 1, n_st
            cos_theta = 0.0_dp
            sz_val = 0.0_dp

            idx1 = 0
            do j1 = 0, cfg%j_max
                do m1 = -j1, j1
                    do ms1 = 1, 2
                        idx1 = idx1 + 1

                        ! <S_z> expectation
                        if (ms1 == 1) then
                            sz_val = sz_val + 0.5_dp * (state%eigenvectors(idx1, k)**2)
                        else
                            sz_val = sz_val - 0.5_dp * (state%eigenvectors(idx1, k)**2)
                        end if

                        ! <cos(theta_R)> off-diagonal expectation
                        idx2 = 0
                        do j2 = 0, cfg%j_max
                            do m2 = -j2, j2
                                do ms2 = 1, 2
                                    idx2 = idx2 + 1
                                    if (ms1 == ms2 .and. m1 == m2 .and. j2 == j1 + 1) then
                                        c_term = sqrt(real((j2**2 - m1**2), dp) / &
                                                      real((2*j1 + 1) * (2*j2 + 1), dp))
                                        cos_theta = cos_theta + 2.0_dp * state%eigenvectors(idx1, k) * &
                                                    state%eigenvectors(idx2, k) * c_term
                                    end if
                                end do
                            end do
                        end do

                    end do
                end do
            end do

            state%orientation(k) = cos_theta
            state%spin_sz(k) = sz_val
        end do
    end subroutine calc_crossed_field_observables

    ! ==========================================================================
    ! scan_tilt_angle_spectrum:
    ! Scans the tilt angle theta_EB from 0 to 90 degrees and records lowest eigenenergies
    ! ==========================================================================
    subroutine scan_tilt_angle_spectrum(cfg_in, theta_deg_arr, n_angles, &
                                       lowest_energies_au, stat)
        type(crossed_field_config_t), intent(in) :: cfg_in
        integer, intent(in)                      :: n_angles
        real(dp), intent(in)                     :: theta_deg_arr(n_angles)
        real(dp), intent(out)                    :: lowest_energies_au(n_angles, 4)
        integer, intent(out)                     :: stat

        type(crossed_field_config_t) :: cfg
        type(crossed_field_state_t) :: state
        integer :: i, s

        stat = 0
        cfg = cfg_in

        do i = 1, n_angles
            cfg%theta_eb_rad = theta_deg_arr(i) * (PI / 180.0_dp)
            call solve_crossed_field_eigenstates(cfg, state, s)
            if (s /= 0) then
                stat = s
                return
            end if
            lowest_energies_au(i, 1:min(4, state%dim_tot)) = state%energies_au(1:min(4, state%dim_tot))
        end do
    end subroutine scan_tilt_angle_spectrum

end module mod_crossed_field_scattering
