"""
pygenmod: Python Interface and High-Level Scientific Tools for GeneralModule
"""

from .constants import (
    PI, TWOPI, HALFPI, SQRTPI,
    C_LIGHT, HBAR, H_PLANCK, M_E, CHARGE_E, EPS0, KB, AMU2AU,
    AU2S, S2AU, AU2FS, FS2AU, AU2PS, PS2AU,
    AU2EV, EV2AU, AU2CM, CM2AU, AU2J, J2AU, AU2K, K2AU,
    AU2M, M2AU, AU2ANG, ANG2AU, AU2NM, NM2AU,
    AU2VM, VM2AU, AU2MV_CM, MV_CM2AU,
    AU2DEBYE, DEBYE2AU, AU2W_CM2, W_CM2AU,
    to_au, from_au
)

from .pulse import (
    PULSE_GAUSSIAN, PULSE_SIN2, PULSE_FLATTOP, PULSE_CHIRP, PULSE_TWOCOLOR, PULSE_THZ_TRAIN,
    PulseConfig, pulse_envelope, pulse_electric_field, pulse_stark_shift
)

from .dvr import (
    SincDVR, dvr_sinc_init, fgh_solve_bound_states
)

from .coulomb import (
    AtomConfig, get_atom_config, soft_core_coulomb_potential, soft_core_coulomb_derivative,
    keldysh_parameter, ponderomotive_energy, hhg_cutoff_energy, quiver_radius, adk_ionization_rate
)

from .hhg import (
    calculate_dipole_acceleration, hhg_power_spectrum
)

from .multistate import (
    landau_zener_probability, calculate_channel_populations
)

from .rovibrational import (
    rovibrational_state_index, rovibrational_state_unindex,
    rot_matrix_cos_theta, calc_franck_condon_factors,
    calc_rotational_constants_bv, build_rovibrational_hamiltonian,
    build_rovibrational_dipole_matrix
)

from .scattering import (
    van_der_waals_mean_length, gribakin_flambaum_length,
    square_well_scattering_length_exact, calc_scattering_length_numerov,
    plot_scattering_length_wavefunction
)

from .visualizer import (
    set_publication_style, plot_wavefunctions, plot_pulses, plot_alignment_dynamics
)

__version__ = "1.3.0"
__all__ = [
    "PI", "TWOPI", "HALFPI", "SQRTPI",
    "to_au", "from_au",
    "PulseConfig", "pulse_envelope", "pulse_electric_field", "pulse_stark_shift",
    "SincDVR", "dvr_sinc_init", "fgh_solve_bound_states",
    "AtomConfig", "get_atom_config", "soft_core_coulomb_potential", "soft_core_coulomb_derivative",
    "keldysh_parameter", "ponderomotive_energy", "hhg_cutoff_energy", "quiver_radius", "adk_ionization_rate",
    "calculate_dipole_acceleration", "hhg_power_spectrum",
    "landau_zener_probability", "calculate_channel_populations",
    "rovibrational_state_index", "rovibrational_state_unindex",
    "rot_matrix_cos_theta", "calc_franck_condon_factors",
    "calc_rotational_constants_bv", "build_rovibrational_hamiltonian",
    "build_rovibrational_dipole_matrix",
    "van_der_waals_mean_length", "gribakin_flambaum_length",
    "square_well_scattering_length_exact", "calc_scattering_length_numerov",
    "plot_scattering_length_wavefunction",
    "plot_wavefunctions", "plot_pulses", "plot_alignment_dynamics"
]
