!> \brief GeneralModule 统一入口总模块
!> \details 汇聚并重新导出算法库内的全部子模块，允许开发者或 AI 直接以 `use general_module` 访问全部接口。
!> \author LiHao
module general_module
    use mod_constants
    use mod_special_functions
    use mod_linear_algebra
    use mod_dvr_grid
    use mod_laser_pulse
    use mod_absorbing_boundary
    use mod_thermal_ensemble
    use mod_wavepacket_propagator
    use mod_coulomb_atomic
    use mod_hhg_spectra
    use mod_chebyshev_propagator
    use mod_multistate_coupling
    use mod_rovibrational
    use mod_io_utils
    implicit none

    public :: dp, int32, int64
    public :: PI, TWOPI, HALFPI, SQRTPI, EYE
    public :: to_au, from_au

    public :: legendre_poly, assoc_legendre_poly, wigner_3j, clebsch_gordan
    public :: rot_matrix_cos_theta, rot_matrix_cos2_theta

    public :: diag_symmetric_matrix, fft_1d, fft_2d

    public :: dvr_1d_t, dvr_legendre_t, dvr_sinc_init, dvr_legendre_init, fgh_solve_bound_states
    public :: dvr_expectation_value, dvr_matrix_element

    public :: pulse_config_t, pulse_envelope, pulse_electric_field, pulse_electric_field_2d
    public :: pulse_vector_potential, pulse_stark_shift, pulse_generate_timeseries
    public :: create_gaussian_pulse, create_sin2_pulse, create_chirped_pulse
    public :: PULSE_GAUSSIAN, PULSE_SIN2, PULSE_FLATTOP, PULSE_CHIRP, PULSE_TWOCOLOR, PULSE_THZ_TRAIN

    public :: absorbing_boundary_t, cap_init, cap_evaluate, cap_apply_mask, calculate_probability_flux, calculate_norm_inside
    public :: CAP_SIN2, CAP_POLYNOMIAL

    public :: boltzmann_rotational_weights, boltzmann_vibrational_weights, thermal_average_1d, thermal_average_2d, bose_einstein_factor

    public :: propagate_split_operator_1d, propagate_split_operator_2d, rk4_step, solve_bloch_two_level, abm4_step

    ! 强场原子物理模型
    public :: atom_config_t, get_atom_config
    public :: soft_core_coulomb_potential, soft_core_coulomb_derivative
    public :: keldysh_parameter, ponderomotive_energy, hhg_cutoff_energy, quiver_radius, adk_ionization_rate

    ! 高次谐波与光子能谱
    public :: calculate_dipole_length, calculate_dipole_acceleration, hhg_power_spectrum
    public :: gabor_transform_point, lewenstein_sfa_dipole

    ! 切比雪夫与能谱滤波
    public :: chebyshev_propagate_step, window_operator_pes

    ! 多态非绝热动力学
    public :: propagate_split_operator_2channel, landau_zener_probability, calculate_channel_populations

    ! 分子转振动力学与态跃迁
    public :: calc_franck_condon_factors, calc_vibrational_dipole_matrix, calc_rotational_constants_bv
    public :: rovibrational_state_index, rovibrational_state_unindex
    public :: build_rovibrational_hamiltonian, build_rovibrational_dipole_matrix
    public :: build_rovibrational_polarizability_matrix, create_stirap_pulses

    ! 标准化科学 I/O 工具
    public :: save_data_table_1d, save_data_table_2d, save_matrix_dat
    public :: print_banner, print_progress_bar

end module general_module
