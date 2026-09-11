"""
Physical constants and unit conversion utilities in atomic units (a.u.) and SI units.
Directly mirrors `mod_constants.f90` for seamless Fortran <-> Python interoperability.
"""
import math

# Mathematical constants
PI = math.pi
TWOPI = 2.0 * math.pi
HALFPI = 0.5 * math.pi
SQRTPI = math.sqrt(math.pi)

# Fundamental constants (CODATA 2018/2022)
C_LIGHT = 2.99792458e8         # Speed of light in vacuum (m/s)
HBAR = 1.054571817e-34         # Reduced Planck constant (J*s)
H_PLANCK = 6.62607015e-34      # Planck constant (J*s)
M_E = 9.1093837015e-31         # Electron rest mass (kg)
CHARGE_E = 1.602176634e-19     # Elementary charge (C)
EPS0 = 8.8541878128e-12        # Vacuum permittivity (F/m)
KB = 1.380649e-23              # Boltzmann constant (J/K)
AMU2AU = 1822.888486209        # 1 amu in electron masses (m_e)

# Time conversions
AU2S = 2.4188843265857e-17     # a.u. -> s
S2AU = 1.0 / AU2S
AU2FS = 2.4188843265857e-2     # a.u. -> fs
FS2AU = 1.0 / AU2FS
AU2PS = 2.4188843265857e-5     # a.u. -> ps
PS2AU = 1.0 / AU2PS

# Energy conversions
AU2EV = 27.211386245988        # Hartree -> eV
EV2AU = 1.0 / AU2EV
AU2CM = 219474.63136320        # Hartree -> cm^-1
CM2AU = 1.0 / AU2CM
AU2J = 4.3597447222071e-18     # Hartree -> J
J2AU = 1.0 / AU2J
AU2K = 3.1577502480407e5       # Hartree -> K
K2AU = 1.0 / AU2K

# Length conversions
AU2M = 0.529177210903e-10      # Bohr -> m
M2AU = 1.0 / AU2M
AU2ANG = 0.529177210903        # Bohr -> Angstrom
ANG2AU = 1.0 / AU2ANG
AU2NM = 0.0529177210903        # Bohr -> nm
NM2AU = 1.0 / AU2NM

# Electric field conversions
AU2VM = 5.14220674763e11       # a.u. -> V/m
VM2AU = 1.0 / AU2VM
AU2MV_CM = 5142.20674763       # a.u. -> MV/cm
MV_CM2AU = 1.0 / AU2MV_CM

# Dipole moment conversions
AU2DEBYE = 2.541746473         # a.u. -> Debye
DEBYE2AU = 1.0 / AU2DEBYE

# Intensity conversions
AU2W_CM2 = 3.5094452e16        # a.u. -> W/cm^2
W_CM2AU = 1.0 / AU2W_CM2

_CONVERSIONS_TO_AU = {
    'fs': FS2AU,
    'ps': PS2AU,
    's': S2AU,
    'ev': EV2AU,
    'cm-1': CM2AU,
    'cm^-1': CM2AU,
    'j': J2AU,
    'k': K2AU,
    'ang': ANG2AU,
    'angstrom': ANG2AU,
    'nm': NM2AU,
    'm': M2AU,
    'amu': AMU2AU,
    'v/m': VM2AU,
    'mv/cm': MV_CM2AU,
    'debye': DEBYE2AU,
    'w/cm2': W_CM2AU,
    'w/cm^2': W_CM2AU,
}

_CONVERSIONS_FROM_AU = {
    'fs': AU2FS,
    'ps': AU2PS,
    's': AU2S,
    'ev': AU2EV,
    'cm-1': AU2CM,
    'cm^-1': AU2CM,
    'j': AU2J,
    'k': AU2K,
    'ang': AU2ANG,
    'angstrom': AU2ANG,
    'nm': AU2NM,
    'm': AU2M,
    'amu': 1.0 / AMU2AU,
    'v/m': AU2VM,
    'mv/cm': AU2MV_CM,
    'debye': AU2DEBYE,
    'w/cm2': AU2W_CM2,
    'w/cm^2': AU2W_CM2,
}


def to_au(val: float, unit_name: str) -> float:
    """Convert scalar or array value from given unit to atomic units (a.u.)."""
    unit_lower = unit_name.strip().lower()
    factor = _CONVERSIONS_TO_AU.get(unit_lower, 1.0)
    return val * factor


def from_au(val_au: float, unit_name: str) -> float:
    """Convert value from atomic units (a.u.) to specified target unit."""
    unit_lower = unit_name.strip().lower()
    factor = _CONVERSIONS_FROM_AU.get(unit_lower, 1.0)
    return val_au * factor
