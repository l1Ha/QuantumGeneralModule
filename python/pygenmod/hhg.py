"""
High harmonic generation (HHG) analysis and dipole calculation utilities.
Mirrors `mod_hhg_spectra.f90`.
"""
import numpy as np
from .constants import TWOPI


def calculate_dipole_acceleration(dv_dx: np.ndarray, dx: float,
                                  psi: np.ndarray, e_field: float) -> float:
    """Calculate instantaneous Ehrenfest dipole acceleration a(t)."""
    norm_sq = np.sum(np.abs(psi)**2) * dx
    force_pot = np.sum(np.abs(psi)**2 * dv_dx) * dx
    return float(-force_pot - e_field * norm_sq)


def hhg_power_spectrum(a_t: np.ndarray, dt: float):
    """
    Compute harmonic power spectrum |a(omega)|^2 using windowed FFT.
    Returns (omega_array, power_spectrum).
    """
    a_t = np.asarray(a_t, dtype=float)
    n = len(a_t)
    window = np.hanning(n)
    signal = a_t * window

    fft_val = np.fft.fft(signal) * dt
    power = np.abs(fft_val[:n // 2])**2

    freqs = np.fft.fftfreq(n, d=dt)[:n // 2] * TWOPI
    return freqs, power
