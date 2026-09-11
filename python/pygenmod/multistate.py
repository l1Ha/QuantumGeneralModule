"""
Multi-state non-adiabatic dynamics and Landau-Zener crossing utilities.
Mirrors `mod_multistate_coupling.f90`.
"""
import numpy as np
from .constants import TWOPI


def landau_zener_probability(v12_crossing: float, velocity: float, delta_slope: float) -> float:
    """Calculate Landau-Zener classic non-adiabatic transition probability."""
    denom = abs(velocity * delta_slope)
    if denom <= 1e-15:
        return 0.0
    return float(np.exp(-TWOPI * (v12_crossing**2) / denom))


def calculate_channel_populations(psi1: np.ndarray, psi2: np.ndarray, dx: float):
    """Calculate integrated norm on channel 1, channel 2, and branch ratio."""
    pop1 = float(np.sum(np.abs(psi1)**2) * dx)
    pop2 = float(np.sum(np.abs(psi2)**2) * dx)
    tot = pop1 + pop2
    ratio = pop2 / tot if tot > 1e-15 else 0.0
    return pop1, pop2, ratio
