!> \brief GeneralModule 统一入口总模块
!> \details 汇聚并重新导出算法库内的全部子模块，
!>          允许开发者或 AI 直接以 `use general_module` 访问全部接口。
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
    use mod_field_scattering
    use mod_dipolar_scattering
    use mod_photoassociation
    use mod_three_body_recombination
    use mod_confined_scattering
    use mod_autoionization_fano
    use mod_crossed_field_scattering
    use mod_triatomic_geometry
    use mod_spinor_bec
    use mod_hyperspherical_reactive
    use mod_dipolar_droplets_lhy
    use mod_strong_field_nsdi
    use mod_feshbach_bound_states
    implicit none

    public :: dp, int32, int64
    public :: PI, TWOPI, HALFPI, SQRTPI, EYE
    public :: to_au, from_au
    public :: GAUSS2AU, AU2GAUSS, AU2TESLA, TESLA2AU

    public :: legendre_poly, assoc_legendre_poly, wigner_3j, clebsch_gordan
    public :: wigner_3j_half, clebsch_gordan_half, wigner_6j_half, wigner_9j_half
    public :: rot_matrix_cos_theta, rot_matrix_cos2_theta

    public :: diag_symmetric_matrix, fft_1d, fft_2d, inv_real_matrix, inv_complex_matrix

    public :: dvr_1d_t, dvr_legendre_t, dvr_sinc_init, dvr_legendre_init, fgh_solve_bound_states
    public :: dvr_expectation_value, dvr_matrix_element

    public :: pulse_config_t, pulse_envelope, pulse_electric_field, pulse_electric_field_2d
    public :: pulse_vector_potential, pulse_stark_shift, pulse_generate_timeseries
    public :: create_gaussian_pulse, create_sin2_pulse, create_chirped_pulse
    public :: PULSE_GAUSSIAN, PULSE_SIN2, PULSE_FLATTOP, PULSE_CHIRP, PULSE_TWOCOLOR, PULSE_THZ_TRAIN

    public :: absorbing_boundary_t, cap_init, cap_evaluate, cap_apply_mask
    public :: calculate_probability_flux, calculate_norm_inside
    public :: CAP_SIN2, CAP_POLYNOMIAL

    public :: boltzmann_rotational_weights, boltzmann_vibrational_weights
    public :: thermal_average_1d, thermal_average_2d, bose_einstein_factor

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
    public :: scattering_state_t, ere_result_t, resonance_info_t, multichannel_result_t
    public :: PARTICLE_DISTINGUISHABLE, PARTICLE_IDENTICAL_BOSON
    public :: PARTICLE_IDENTICAL_FERMION_POLARIZED, PARTICLE_IDENTICAL_FERMION_UNPOLAR
    public :: NORM_ENERGY, NORM_MOMENTUM, NORM_UNIT_AMPLITUDE
    public :: riccati_bessel_neumann
    public :: calc_scattering_length_numerov, calc_scattering_length_logder
    public :: calc_phase_shift_single_l, calc_scattering_wavefunction_ti
    public :: calc_scattering_wavefunction_1d_cartesian
    public :: calc_partial_wave_cross_sections
    public :: optical_theorem_cross_section, calc_differential_cross_section
    public :: calc_differential_cross_section_identical
    public :: calc_transport_cross_sections, calc_cross_section_spectrum
    public :: calc_generalized_cross_sections, calc_differential_legendre_expansion
    public :: fit_effective_range_expansion
    public :: gribakin_flambaum_length, van_der_waals_mean_length
    public :: analyze_shape_resonance, calc_coupled_channel_smatrix_2x2
    public :: calc_multichannel_close_coupling_logder, calc_feshbach_resonance_scan
    public :: segmented_grid_t, create_segmented_grid
    public :: calc_scattering_length_segmented_numerov, calc_scattering_wavefunction_segmented_ti
    public :: calc_phase_shift_segmented, calc_multichannel_close_coupling_segmented_logder

    ! 含时波包散射理论与 S-矩阵 (Time-Dependent Scattering)
    public :: td_scattering_channel_t
    public :: gaussian_wavepacket_1d, gaussian_momentum_amplitude
    public :: accumulate_flux_amplitude, calculate_td_transmission
    public :: calculate_td_smatrix_element, project_wavepacket_to_smatrix
    public :: multichannel_td_smatrix_elements
    public :: wavepacket_centroid_position, wavepacket_wigner_delay
    public :: calculate_td_differential_cross_section_2d
    public :: accumulate_wavefunction_spectral_projection, extract_td_scattering_wavefunction

    ! 外加电磁场超冷量子碰撞散射 (Field-Dressed & Multi-Basis Scattering)
    public :: BASIS_UNCOUPLED, BASIS_F_COUPLED, BASIS_TOTAL_SPIN, BASIS_FIELD_DRESSED
    public :: MU_B_AU, MU_N_AU, GHZ2AU, AU2GHZ
    public :: cold_atom_t, field_channel_t, field_feshbach_result_t
    public :: get_cold_atom_preset, calc_breit_rabi_energies
    public :: build_field_collision_channels, calc_basis_transform_matrix
    public :: build_asymptotic_hamiltonian, build_spin_exchange_matrix
    public :: build_field_potential_matrix, calc_magnetic_feshbach_resonance_scan
    public :: fit_feshbach_resonance_parameters

    ! 各向异性磁偶极与电偶极超冷散射 (Dipolar Quantum Scattering)
    public :: FINE_STRUCT_ALPHA, ELECTRON_GS, BOHR_MAGNETON_AU
    public :: polar_molecule_t, dipolar_channel_t, dipolar_relaxation_result_t
    public :: calc_mddi_spatial_matrix_element, calc_mddi_spin_matrix_element
    public :: calc_mddi_total_matrix_element, calc_dipolar_relaxation_cross_section
    public :: calc_dipolar_relaxation_thermal_rate, calc_stark_induced_dipole
    public :: calc_eddi_matrix_element, calc_dipole_length_scale
    public :: build_dipolar_coupled_potential_matrix

    ! 超冷光缔合谱学与自由-束缚态量子跃迁 (Ultracold Photoassociation)
    public :: pa_config_t, pa_transition_result_t
    public :: calc_free_bound_fc_overlap, calc_pa_stimulated_width
    public :: calc_pa_cross_section, calc_pa_thermal_rate_coefficient
    public :: calc_pa_spectrum_scan, calc_two_photon_raman_association_coupling

    ! 超冷三体复合与 Efimov 少体物理 (Three-Body Recombination & Efimov Physics)
    public :: efimov_param_t, three_body_loss_t
    public :: solve_efimov_s0_identical_bosons, solve_efimov_s0_heteronuclear
    public :: calc_efimov_scale_factor, calc_three_body_recombination_a_positive
    public :: calc_three_body_recombination_a_negative, calc_three_body_recombination_universal
    public :: calc_unitary_three_body_loss_temperature, scan_efimov_recombination_spectrum

    ! 低维微观光阱受限量子散射与约束诱导共振 (Confinement-Induced Resonances)
    public :: waveguide_1d_t, planar_2d_t, cir_result_t
    public :: init_waveguide_1d, init_planar_2d, calc_olshanii_cir_parameters
    public :: calc_confined_dimer_binding_energy, calc_anisotropic_cir_poles
    public :: calc_quasi_2d_scattering, calc_lieb_liniger_parameter

    ! 自电离体系与 Fano 共振理论 (Autoionization & Complex Coordinate Rotation)
    public :: fano_profile_t, ccr_resonance_t
    public :: calc_fano_profile, calc_fano_cross_section_spectrum
    public :: calc_autoionization_lifetime, calc_fano_q_parameter
    public :: solve_ccr_resonance_model, calc_time_domain_autoionization_decay

    ! 交叉外电磁场强耦合量子散射 (Crossed Electric & Magnetic Fields)
    public :: crossed_field_config_t, crossed_field_state_t
    public :: init_crossed_field_config, build_crossed_field_hamiltonian
    public :: solve_crossed_field_eigenstates, calc_crossed_field_observables
    public :: scan_tilt_angle_spectrum

    ! 三原子反应散射、Jacobi 坐标与圆锥交叉几何相位 (Triatomic & Geometric Phase)
    public :: jacobi_coord_t, internuclear_dist_t, leps_param_t, conical_intersection_t
    public :: jacobi_to_internuclear, internuclear_to_jacobi
    public :: calc_leps_potential, init_default_h3_leps
    public :: calc_conical_intersection_adiabats, calc_berry_phase_around_ci

    ! 超冷旋量玻色爱因斯坦凝聚与多体自旋动力学 (Spinor BEC Dynamics)
    public :: spinor_param_t, spinor_state_t
    public :: init_spinor_preset, calc_spinor_interaction_couplings
    public :: calc_quadratic_zeeman_shift, propagate_spinor_sma_rk4
    public :: calc_spinor_energy, simulate_spin_mixing_dynamics

    ! 三原子超球面反应散射动力学 (Hyperspherical Reactive Scattering)
    public :: reaction_mass_t, transition_state_t
    public :: init_reaction_mass, jacobi_to_hyperspherical, hyperspherical_to_jacobi
    public :: calc_eckart_transmission, calc_cumulative_reaction_probability
    public :: calc_canonical_rate_constant, calc_tst_wigner_rate

    ! 超冷偶极量子液滴与李-黄-杨量子涨落修正 (Dipolar Droplets & LHY)
    public :: dipolar_droplet_param_t
    public :: init_dipolar_droplet_param, calc_pelster_lima_q5
    public :: calc_equilibrium_droplet_density, calc_droplet_chemical_potential
    public :: calc_critical_atom_number, calc_egpe_energy_density

    ! 强场非顺序双电离与电子重碰撞动力学 (Strong-Field NSDI & Recollision)
    public :: nsdi_laser_t, nsdi_target_t, nsdi_result_t
    public :: init_nsdi_laser, init_nsdi_target
    public :: calc_ponderomotive_energy, calc_keldysh_gamma, calc_adk_rate
    public :: calc_recollision_trajectory, calc_lotz_cross_section
    public :: calc_nsdi_drift_momenta, calc_nsdi_2d_momentum_dist
    public :: calc_double_ion_yield_curve

    ! 磁与光 Feshbach 共振与弱束缚分子态 (Feshbach Bound States & Loss)
    public :: mfr_param_t, ofr_param_t
    public :: init_mfr_preset, init_ofr_param
    public :: calc_mfr_scattering_length, calc_mfr_r_star, calc_mfr_s_res
    public :: calc_mfr_bound_energy_universal, calc_mfr_bound_energy_coupled
    public :: calc_mfr_closed_channel_fraction
    public :: calc_ofr_complex_scattering_length, calc_ofr_inelastic_loss_rate

end module general_module
