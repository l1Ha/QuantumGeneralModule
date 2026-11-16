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


def adaptive_split_propagate(psi, v_pot, dx, mass, t_total, tol,
                             order4=False, h0=1e-2, h_min=1e-8, h_max=1.0):
    """Adaptive-step split propagation with local error control.

    Mirrors `adaptive_split_propagate_1d` (src/mod_wavepacket_propagator.f90):
    step doubling estimates the local error (one step h vs two steps h/2);
    a PI-type controller  h_new = h * clamp((tol/err)**(1/(p+1)), 0.2, 5)
    accepts/rejects steps (p = 2 for Strang, 4 for Yoshida4).

    Returns
    -------
    psi : propagated state
    info : dict with n_steps, n_reject and the per-accepted-step h history.
    """
    step = split_yoshida4_step if order4 else split_strang_step
    p_ord = 4 if order4 else 2
    safety, fac_min, fac_max = 0.9, 0.2, 5.0
    h = min(max(h0, h_min), h_max)
    t_now = 0.0
    n_acc = n_rej = 0
    hist = []
    while t_now < t_total - 1e-15:
        if t_now + h > t_total:
            h = t_total - t_now
        a = step(psi, v_pot, dx, mass, h)
        b = step(step(psi, v_pot, dx, mass, 0.5 * h), v_pot, dx, mass, 0.5 * h)
        err = float(np.linalg.norm(a - b) * np.sqrt(dx))
        if err <= tol or h <= h_min * 1.000001:
            psi = b
            t_now += h
            n_acc += 1
            hist.append(h)
            ratio = safety * (tol / max(err, 1e-300)) ** (1.0 / (p_ord + 1)) if err > 0 else fac_max
            h = min(max(h * min(max(ratio, fac_min), fac_max), h_min), h_max)
        else:
            n_rej += 1
            ratio = safety * (tol / err) ** (1.0 / (p_ord + 1))
            h = min(max(h * min(max(ratio, fac_min), 0.9), h_min), h_max)
    return psi, {"n_steps": n_acc, "n_reject": n_rej, "h_history": hist}
