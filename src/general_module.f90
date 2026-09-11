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
    use mod_interpolation
    use mod_photofragment_flux
    use mod_open_quantum
    use mod_optimal_control
    use mod_ti_scattering
    use mod_td_scattering
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

    ! 势能面高精度样条插值与外推
    public :: spline_1d_t, spline_init, spline_eval, spline_eval_deriv, spline_eval_deriv2
    public :: spline_clean, interpolate_pes_to_grid, BC_NATURAL, BC_CLAMPED

    ! 光解离碎片能谱与自相关吸收谱
    public :: calculate_autocorrelation, calculate_absorption_spectrum
    public :: energy_resolved_flux_amplitude, calculate_ker_spectrum, calculate_branching_ratios

    ! 开放量子系统与耗散动力学
    public :: lindblad_dissipator, rk4_lindblad_step, calculate_quantum_purity
    public :: calculate_von_neumann_entropy, calculate_quantum_coherence
    public :: create_relaxation_jump_op, create_dephasing_jump_op

    ! 量子最优控制理论 (Krotov 算法)
    public :: oct_fidelity, oct_shape_function, oct_krotov_step, oct_optimize_pulse

    ! 非含时散射理论与超冷碰撞 (Time-Independent Scattering)
    public :: scattering_state_t, ere_result_t, resonance_info_t
    public :: PARTICLE_DISTINGUISHABLE, PARTICLE_IDENTICAL_BOSON
    public :: PARTICLE_IDENTICAL_FERMION_POLARIZED, PARTICLE_IDENTICAL_FERMION_UNPOLAR
    public :: riccati_bessel_neumann
    public :: calc_scattering_length_numerov, calc_scattering_length_logder
    public :: calc_phase_shift_single_l, calc_partial_wave_cross_sections
    public :: optical_theorem_cross_section, calc_differential_cross_section
    public :: calc_differential_cross_section_identical
    public :: calc_transport_cross_sections, calc_cross_section_spectrum
    public :: calc_generalized_cross_sections, calc_differential_legendre_expansion
    public :: fit_effective_range_expansion
    public :: gribakin_flambaum_length, van_der_waals_mean_length
    public :: analyze_shape_resonance, calc_coupled_channel_smatrix_2x2

    ! 含时波包散射理论与 S-矩阵 (Time-Dependent Scattering)
    public :: td_scattering_channel_t
    public :: gaussian_wavepacket_1d, gaussian_momentum_amplitude
    public :: accumulate_flux_amplitude, calculate_td_transmission
    public :: calculate_td_smatrix_element, project_wavepacket_to_smatrix
    public :: multichannel_td_smatrix_elements
    public :: wavepacket_centroid_position, wavepacket_wigner_delay
    public :: calculate_td_differential_cross_section_2d

end module general_module
