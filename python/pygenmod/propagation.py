"""
Split-operator propagation with higher-order compositions.

Mirrors `mod_wavepacket_propagator` (src/mod_wavepacket_propagator.f90):
the second-order Strang step  S2(dt) = e^{-iV dt/2} F^{-1} e^{-iT dt} F e^{-iV dt/2}
and the fourth-order symmetric compositions of Suzuki (5 sub-steps) and
Yoshida (3 sub-steps). FFT kinetic operator by numpy.fft, so the discrete
Hamiltonian matches the Fortran split-operator exactly.
"""
import numpy as np

__all__ = ["split_strang_step", "split_suzuki4_step", "split_yoshida4_step"]

#: Suzuki 1990 fractal composition: S4 = S2(p)^2 S2((1-4p)) S2(p)^2
SUZUKI_P = 1.0 / (4.0 - 4.0 ** (1.0 / 3.0))          # 0.4144907717943757
#: Yoshida 1990 three-factor composition: S4 = S2(x1) S2(x0) S2(x1)
YOSHIDA_X1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))        # 1.3512071919596578
YOSHIDA_X0 = -(2.0 ** (1.0 / 3.0)) * YOSHIDA_X1      # -1.7024143839193153


def split_strang_step(psi, v_pot, dx, mass, dt):
    """One Strang step: psi <- exp(-i H dt) psi (H = FFT kinetic + V)."""
    n = psi.size
    k = np.fft.fftfreq(n, d=dx / (2.0 * np.pi))      # matches Fortran p indexing
    psi = psi * np.exp(-0.5j * v_pot * dt)
    psi = np.fft.ifft(np.exp(-1j * (k ** 2) * dt / (2.0 * mass)) * np.fft.fft(psi))
    psi = psi * np.exp(-0.5j * v_pot * dt)
    return psi


def split_suzuki4_step(psi, v_pot, dx, mass, dt):
    """One Suzuki 4th-order step (5 Strang sub-steps)."""
    w1 = SUZUKI_P * dt
    w0 = (1.0 - 4.0 * SUZUKI_P) * dt
    for w in (w1, w1, w0, w1, w1):
        psi = split_strang_step(psi, v_pot, dx, mass, w)
    return psi


def split_yoshida4_step(psi, v_pot, dx, mass, dt):
    """One Yoshida 4th-order step (3 Strang sub-steps)."""
    for w in (YOSHIDA_X1 * dt, YOSHIDA_X0 * dt, YOSHIDA_X1 * dt):
        psi = split_strang_step(psi, v_pot, dx, mass, w)
    return psi
