"""
Laser pulse synthesizers and field generators.
Mirrors `mod_laser_pulse.f90`.
"""
import numpy as np
from dataclasses import dataclass
from .constants import PI, TWOPI, HALFPI, FS2AU


PULSE_GAUSSIAN = 1
PULSE_SIN2 = 2
PULSE_FLATTOP = 3
PULSE_CHIRP = 4
PULSE_TWOCOLOR = 5
PULSE_THZ_TRAIN = 6


@dataclass
class PulseConfig:
    shape_type: int = PULSE_GAUSSIAN
    field_peak: float = 0.05
    freq_central: float = 0.057
    duration: float = 50.0 * FS2AU
    t_center: float = 0.0
    chirp_rate: float = 0.0
    cep_phase: float = 0.0
    t_flattop: float = 0.0
    two_color_ratio: float = 0.0
    two_color_phase: float = 0.0
    train_count: int = 1
    train_delay: float = 0.0


def pulse_envelope(t: np.ndarray, cfg: PulseConfig) -> np.ndarray:
    """Compute instantaneous envelope f(t) in [0, 1]. Vectorized over numpy array."""
    t = np.asarray(t, dtype=float)
    dt = t - cfg.t_center
    tau = max(1e-12, cfg.duration)

    if cfg.shape_type in (PULSE_GAUSSIAN, PULSE_CHIRP, PULSE_TWOCOLOR):
        return np.exp(-2.772588722239781 * (dt / tau)**2)

    elif cfg.shape_type == PULSE_SIN2:
        mask = np.abs(dt) <= tau
        env = np.zeros_like(dt)
        env[mask] = np.sin(PI * (dt[mask] + tau) / (2.0 * tau))**2
        return env

    elif cfg.shape_type == PULSE_FLATTOP:
        env = np.zeros_like(dt)
        abs_dt = np.abs(dt)
        half_flat = 0.5 * cfg.t_flattop
        t_ramp = tau

        flat_mask = abs_dt <= half_flat
        ramp_mask = (abs_dt > half_flat) & (abs_dt <= half_flat + t_ramp)

        env[flat_mask] = 1.0
        env[ramp_mask] = 0.5 * (1.0 + np.cos(PI * (abs_dt[ramp_mask] - half_flat) / t_ramp))
        return env

    elif cfg.shape_type == PULSE_THZ_TRAIN:
        env = np.zeros_like(dt)
        for k in range(cfg.train_count):
            k_dt = t - (cfg.t_center + k * cfg.train_delay)
            env += np.exp(-2.772588722239781 * (k_dt / tau)**2)
        return np.clip(env, 0.0, 1.0)

    else:
        return np.exp(-2.772588722239781 * (dt / tau)**2)


def pulse_electric_field(t: np.ndarray, cfg: PulseConfig) -> np.ndarray:
    """Compute instantaneous electric field E(t) (a.u.). Vectorized."""
    t = np.asarray(t, dtype=float)
    dt = t - cfg.t_center
    env = pulse_envelope(t, cfg)

    if cfg.shape_type == PULSE_CHIRP:
        phase = cfg.freq_central * dt + 0.5 * cfg.chirp_rate * dt**2 + cfg.cep_phase
        return cfg.field_peak * env * np.cos(phase)

    elif cfg.shape_type == PULSE_TWOCOLOR:
        e1 = cfg.field_peak * env * np.cos(cfg.freq_central * dt + cfg.cep_phase)
        e2 = cfg.field_peak * cfg.two_color_ratio * env * np.cos(
            2.0 * cfg.freq_central * dt + cfg.two_color_phase
        )
        return e1 + e2

    elif cfg.shape_type == PULSE_THZ_TRAIN:
        efield = np.zeros_like(dt)
        for k in range(cfg.train_count):
            k_dt = t - (cfg.t_center + k * cfg.train_delay)
            env_k = np.exp(-2.772588722239781 * (k_dt / cfg.duration)**2)
            efield += cfg.field_peak * env_k * np.sin(cfg.freq_central * k_dt + cfg.cep_phase)
        return efield

    else:
        return cfg.field_peak * env * np.cos(cfg.freq_central * dt + cfg.cep_phase)


def pulse_stark_shift(t: np.ndarray, cfg: PulseConfig,
                      alpha_parallel: float, alpha_perp: float,
                      theta: float) -> np.ndarray:
    """Compute instantaneous AC Stark shift Delta_S(t, theta) (a.u.)."""
    ef = pulse_electric_field(t, cfg)
    costh = np.cos(theta)
    sinth = np.sin(theta)
    alpha_eff = alpha_parallel * costh**2 + alpha_perp * sinth**2
    return -0.25 * alpha_eff * (ef**2)
