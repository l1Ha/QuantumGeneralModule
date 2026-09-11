"""
Strong field atomic models, soft-core potentials, and strong-field ionization parameters.
Mirrors `mod_coulomb_atomic.f90`.
"""
import numpy as np
import math
from dataclasses import dataclass
from .constants import AU2EV, EV2AU


@dataclass
class AtomConfig:
    name: str = "H"
    ip_au: float = 0.5
    soft_core_a: float = 1.0
    z_eff: float = 1.0
    l_quantum: int = 0
    m_quantum: int = 0


_ATOMS = {
    'H': AtomConfig("H", ip_au=0.5, soft_core_a=1.0, z_eff=1.0, l_quantum=0, m_quantum=0),
    'HE': AtomConfig("He", ip_au=0.90357, soft_core_a=0.751, z_eff=1.0, l_quantum=0, m_quantum=0),
    'NE': AtomConfig("Ne", ip_au=0.79248, soft_core_a=0.885, z_eff=1.0, l_quantum=1, m_quantum=0),
    'AR': AtomConfig("Ar", ip_au=0.57915, soft_core_a=1.340, z_eff=1.0, l_quantum=1, m_quantum=0),
    'KR': AtomConfig("Kr", ip_au=0.51446, soft_core_a=1.480, z_eff=1.0, l_quantum=1, m_quantum=0),
    'XE': AtomConfig("Xe", ip_au=0.44577, soft_core_a=1.620, z_eff=1.0, l_quantum=1, m_quantum=0),
}


def get_atom_config(name: str) -> AtomConfig:
    """Retrieve pre-configured SAE parameters for common noble gas atoms and hydrogen."""
    key = name.strip().upper()
    return _ATOMS.get(key, AtomConfig("H", ip_au=0.5, soft_core_a=1.0))


def soft_core_coulomb_potential(x: np.ndarray, soft_a: float, z_eff: float = 1.0) -> np.ndarray:
    """1D soft-core Coulomb potential: V(x) = -Z / sqrt(x^2 + a^2)."""
    x = np.asarray(x, dtype=float)
    return -z_eff / np.sqrt(x**2 + soft_a**2)


def soft_core_coulomb_derivative(x: np.ndarray, soft_a: float, z_eff: float = 1.0) -> np.ndarray:
    """1D soft-core Coulomb derivative: dV/dx = Z * x / (x^2 + a^2)^(3/2)."""
    x = np.asarray(x, dtype=float)
    return z_eff * x / ((x**2 + soft_a**2)**1.5)


def keldysh_parameter(omega: float, e_peak: float, ip_au: float) -> float:
    """Keldysh parameter: gamma = omega * sqrt(2 * Ip) / E_0."""
    if e_peak <= 1e-15:
        return 1e10
    return omega * np.sqrt(2.0 * ip_au) / e_peak


def ponderomotive_energy(e_peak: float, omega: float) -> float:
    """Ponderomotive energy: U_p = E_0^2 / (4 * omega^2)."""
    if omega <= 1e-15:
        return 0.0
    return (e_peak**2) / (4.0 * omega**2)


def hhg_cutoff_energy(ip_au: float, e_peak: float, omega: float) -> float:
    """High harmonic generation cutoff law: E_cutoff = Ip + 3.17 * Up."""
    return ip_au + 3.1725955 * ponderomotive_energy(e_peak, omega)


def quiver_radius(e_peak: float, omega: float) -> float:
    """Excursion amplitude in laser field: alpha_0 = E_0 / omega^2."""
    if omega <= 1e-15:
        return 0.0
    return e_peak / (omega**2)


def adk_ionization_rate(e_field: float, ip_au: float, z_eff: float = 1.0,
                        l: int = 0, m: int = 0) -> float:
    """Ammosov-Delone-Krainov (ADK) tunneling ionization rate for hydrogenic / complex atoms."""
    ef = abs(e_field)
    if ef <= 1e-6 or ip_au <= 1e-6:
        return 0.0

    n_star = z_eff / np.sqrt(2.0 * ip_au)
    c_nl_sq = (2.0**(2.0 * n_star)) / (n_star * math.gamma(n_star + 1.0) * math.gamma(n_star))
    factor_f = (2.0 * l + 1) * math.factorial(l + abs(m)) / (
        (2**abs(m)) * math.factorial(abs(m)) * math.factorial(l - abs(m))
    )

    e_crit = 2.0 * (2.0 * ip_au)**1.5
    w_adk = c_nl_sq * factor_f * ip_au * (
        (e_crit / ef)**(2.0 * n_star - abs(m) - 1.0)
    ) * np.exp(-e_crit / (3.0 * ef))

    return float(w_adk)
