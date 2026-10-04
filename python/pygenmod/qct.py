"""
Quasi-classical trajectory (QCT) reactive scattering for atom-diatom systems.

Python reference implementation of textbook Chapter 17 (准经典轨迹动力学).
It mirrors the recommended Fortran interface `mod_qct_dynamics`
(`qct_config_t`, `qct_trajectory_t`, `qct_result_t`, `qct_init_trajectory`,
`qct_propagate_step`, `qct_analyze_final_state`, `run_qct_ensemble`,
`qct_cross_section`, `qct_thermal_rate`); the Fortran module itself is still on
the roadmap (book section 18.8 / Chapter 12).

All quantities are in atomic units (hbar = 1). The default potential energy
surface is the classic LEPS surface for H + H2 (Karplus-Porter-Sharma), which
is the benchmark of the original QCT literature.

References:
    M. Karplus, R. N. Porter, R. D. Sharma, J. Chem. Phys. 43, 3259 (1965).
    D. G. Truhlar, J. T. Muckerman, "Quasiclassical trajectory calculation of
    cross sections", Methods in Computational Physics 10 (1971).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Tuple

import numpy as np
from scipy.optimize import brentq

from .constants import EV2AU

__all__ = [
    "LEPSParameters", "leps_energy_gradient", "morse_vibrational_energy",
    "ebk_internal_energy",
    "QCTConfig", "QCTTrajectory", "QCTResult", "InitialSampler",
    "qct_init_trajectory", "qct_propagate_step", "propagate_trajectory",
    "qct_analyze_final_state", "run_qct_ensemble",
    "qct_cross_section", "stratified_cross_section", "opacity_function",
    "state_resolved_cross_sections", "differential_cross_section",
    "qct_thermal_rate", "thermal_population", "wilson_interval",
]


# ---------------------------------------------------------------------------
# Potential energy surface: LEPS for H + H2
# ---------------------------------------------------------------------------

@dataclass
class LEPSParameters:
    """LEPS parameters for a symmetric (H3-type) surface.

    d_e   : Morse well depth of one diatomic pair (Hartree)
    beta  : Morse range parameter (bohr^-1)
    r_e   : Morse equilibrium distance (bohr)
    """
    d_e: float = 4.746 * EV2AU
    beta: float = 1.942 * 0.529177  # 1.942 A^-1 in bohr^-1
    r_e: float = 0.7416 / 0.529177  # 0.7416 A in bohr

    def q_and_j(self, r):
        """Coulomb integral Q(r) and exchange integral J(r) of the LEPS form."""
        x = r - self.r_e
        em = np.exp(-self.beta * x)
        q = 0.5 * self.d_e * (1.5 * em * em - em)
        j = 0.25 * self.d_e * (em * em - 6.0 * em)
        return q, j

    def dq_and_dj(self, r):
        """Derivatives dQ/dr and dJ/dr."""
        x = r - self.r_e
        em = np.exp(-self.beta * x)
        dq = 0.5 * self.d_e * self.beta * (em - 3.0 * em * em)
        dj = 0.25 * self.d_e * self.beta * (6.0 * em - 2.0 * em * em)
        return dq, dj


_PAIRS = ((0, 1), (1, 2), (2, 0))


def _pair_distances(q: np.ndarray):
    """Pair distances and unit vectors for atom triplets shaped (..., 3, 3)."""
    diff = q[..., :, None, :] - q[..., None, :, :]          # (..., 3, 3, 3)
    r = np.sqrt(np.einsum("...i,...i->...", diff, diff))     # (..., 3, 3)
    r = np.maximum(r, 1.0e-10)
    unit = diff / r[..., :, None]
    return r, unit


def leps_energy_gradient(q: np.ndarray, par: Optional[LEPSParameters] = None
                         ) -> Tuple[np.ndarray, np.ndarray]:
    """LEPS energy and cartesian forces for one or many atom triplets.

    Parameters
    ----------
    q   : array (..., 3, 3); positions of three atoms on the last two axes.
    par : LEPS parameters (default: H + H2).

    Returns
    -------
    (V, F) with V shape (...) in Hartree and F = -dV/dq of shape (..., 3, 3).

    The London form is
        V = Q1 + Q2 + Q3 - sqrt(Delta),
        Delta = J1^2 + J2^2 + J3^2 - J1 J2 - J2 J3 - J3 J1
              = (1/2)[(J1-J2)^2 + (J2-J3)^2 + (J3-J1)^2] >= 0.
    """
    if par is None:
        par = LEPSParameters()
    q = np.asarray(q, dtype=float)
    r, unit = _pair_distances(q)
    rp = np.stack([r[..., a, b] for a, b in _PAIRS], axis=-1)      # (..., 3)
    qq, jj = par.q_and_j(rp)
    delta = 0.5 * ((jj[..., 0] - jj[..., 1]) ** 2 +
                   (jj[..., 1] - jj[..., 2]) ** 2 +
                   (jj[..., 2] - jj[..., 0]) ** 2)
    sq = np.sqrt(np.maximum(delta, 1.0e-300))
    v = qq.sum(axis=-1) - sq
    dq, dj = par.dq_and_dj(rp)
    # dV/dr_i = Q_i' - (2 J_i - J_j - J_k)/(2 sqrt(Delta)) J_i'
    #   with 2 J_i - J_j - J_k = 3 J_i - (J1+J2+J3)
    coef = (3.0 * jj - jj.sum(axis=-1, keepdims=True)) / (2.0 * sq[..., None])
    dv = dq - coef * dj                                            # (..., 3)
    force = np.zeros_like(q)
    for k, (a, b) in enumerate(_PAIRS):
        g = dv[..., k][..., None]
        force[..., a, :] -= g * unit[..., a, b, :]
        force[..., b, :] += g * unit[..., a, b, :]
    return v, force


def morse_vibrational_energy(v: float, par: Optional[LEPSParameters] = None,
                             mu_r: float = 0.5 * 1837.15) -> float:
    """Exact Morse vibrational energy G(v) (Hartree); v may be non-integer.

    G(v) = w0 (v + 1/2) - beta^2/(2 mu) (v + 1/2)^2 with w0 = beta sqrt(2 D / mu).
    """
    if par is None:
        par = LEPSParameters()
    w0 = par.beta * math.sqrt(2.0 * par.d_e / mu_r)
    w = v + 0.5
    return min(w0 * w - par.beta ** 2 / (2.0 * mu_r) * w * w, par.d_e)


# ---------------------------------------------------------------------------
# Book interface: config / trajectory / result containers (17.8)
# ---------------------------------------------------------------------------

@dataclass
class QCTConfig:
    """Mirrors `qct_config_t` of book section 17.8."""
    e_coll: float = 1.0 * EV2AU     # collision energy (Hartree)
    b_max: float = 3.0              # impact parameter cutoff (bohr)
    r_start: float = 9.0            # initial atom--diatom separation (bohr)
    r_end: float = 11.0             # product separation criterion (bohr)
    dt: float = 4.0                 # Velocity-Verlet time step (a.u.)
    n_traj: int = 200               # trajectories per ensemble
    max_steps: int = 20000          # integration step budget
    v_initial: int = 0              # initial vibrational quantum number
    j_initial: int = 0              # initial rotational quantum number
    seed: int = 20260101            # RNG seed (reproducibility)
    enforce_zpe: bool = False       # passive ZPE constraint at binning
    init_mode: str = "action"       # "action" (EBK) or "energy" (G + F assignment)
    n_b_bins: int = 20              # stratification bins for the opacity function
    mass: float = 1837.15           # atom mass; symmetric H3-type system (m_e)


@dataclass
class QCTTrajectory:
    """Mirrors `qct_trajectory_t` of book section 17.8."""
    q: np.ndarray                   # final positions (3, 3)
    p: np.ndarray                   # final momenta (3, 3)
    time: float = 0.0
    b_impact: float = 0.0
    e_total: float = 0.0
    j_total: np.ndarray = None      # conserved total angular momentum vector
    r_dir_initial: np.ndarray = None  # initial direction of the relative motion
    v_final: int = -1
    j_final: int = -1
    e_internal: float = 0.0
    reactive: bool = False
    converged: bool = False


@dataclass
class QCTResult:
    """Mirrors `qct_result_t` of book section 17.8."""
    cross_section: float = 0.0
    stat_error: float = 0.0
    n_traj: int = 0
    n_reactive: int = 0
    wilson_lo: float = 0.0
    wilson_hi: float = 0.0
    opacity_b: np.ndarray = None
    opacity_p: np.ndarray = None
    state_cross_section: Dict[Tuple[int, int], float] = field(default_factory=dict)
    trajectories: list = None


# ---------------------------------------------------------------------------
# Initial conditions: EBK action assignment and invariant-torus sampling (17.2)
# ---------------------------------------------------------------------------

def _effective_potential(r, j2, par, mu_r):
    em = np.exp(-par.beta * (r - par.r_e))
    return par.d_e * (em * em - 2.0 * em) + j2 / (2.0 * mu_r * r * r)


def _effective_potential_bottom(j2, par, mu_r):
    """Location and value of the minimum of V_eff (shifted outward by the
    centrifugal term, so it is NOT at r_e for j > 0)."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda r: _effective_potential(r, j2, par, mu_r),
                          bounds=(0.3 * par.r_e, 4.0 * par.r_e), method="bounded",
                          options={"xatol": 1.0e-12})
    return float(res.x), float(res.fun)


def _radial_turning_points(e_int, j2, par, mu_r):
    """Inner and outer classical turning points of V_eff(r) = V_Morse + j^2/(2 mu r^2)."""
    f = lambda r: _effective_potential(r, j2, par, mu_r) - e_int
    r_bot, _ = _effective_potential_bottom(j2, par, mu_r)
    r_in = brentq(f, 1.0e-3, r_bot, xtol=1.0e-14, rtol=8.9e-16)
    r_out = r_bot
    while f(r_out) < 0.0 and r_out < 50.0 * par.r_e:
        r_out += 0.05 * par.r_e
    r_out = brentq(f, r_bot, r_out, xtol=1.0e-14, rtol=8.9e-16)
    return r_in, r_out


def _radial_action(e_int, j2, par, mu_r) -> float:
    from scipy.integrate import quad
    r_in, r_out = _radial_turning_points(e_int, j2, par, mu_r)
    fun = lambda r: math.sqrt(max(2.0 * mu_r * (e_int - _effective_potential(r, j2, par, mu_r)), 0.0))
    mid = 0.5 * (r_in + r_out)
    lhs, _ = quad(fun, r_in, mid, limit=200)
    rhs, _ = quad(fun, mid, r_out, limit=200)
    return 2.0 * (lhs + rhs)


def ebk_internal_energy(v: int, j: int, par: Optional[LEPSParameters] = None,
                        mu_r: float = 0.5 * 1837.15) -> float:
    """EBK/WKB action assignment (book section 17.2): solve

        J(E) = 2 int_{r_-}^{r_+} sqrt(2 mu (E - V_eff(r))) dr = 2 pi (v + 1/2).

    For a pure Morse oscillator the EBK condition is exact, so for j = 0 this
    reproduces G(v) to round-off. For j > 0 the centrifugal term couples
    rotation and vibration (this is the approximation the book warns about).
    """
    if par is None:
        par = LEPSParameters()
    j2 = float(j * (j + 1))     # |j|^2 = hbar^2 j(j+1), hbar = 1
    if j2 / (2.0 * mu_r * par.r_e ** 2) >= par.d_e:
        raise ValueError(f"j = {j} has no bound orbit on this Morse potential")
    _, vmin = _effective_potential_bottom(j2, par, mu_r)
    target = 2.0 * math.pi * (v + 0.5)
    fun = lambda e: _radial_action(e, j2, par, mu_r) - target
    e_lo = vmin + 1.0e-9 * par.d_e
    e_hi = -1.0e-7 * abs(par.d_e)
    if fun(e_lo) > 0.0 or fun(e_hi) < 0.0:
        raise ValueError("EBK action bracket failed; state may not exist")
    return brentq(fun, e_lo, e_hi, xtol=1.0e-12 * par.d_e, rtol=1.0e-13)


def _verlet_1d(r, p, dt, n_steps, force_fn):
    """Scalar Velocity-Verlet for the 1D radial orbit (fast pure-Python path)."""
    a = force_fn(r)
    for _ in range(n_steps):
        p += 0.5 * dt * a
        r += dt * p
        a = force_fn(r)
        p += 0.5 * dt * a
    return r, p


def _radial_orbit_table(e_int, j2, par, mu_r, n_table=4096):
    """Integrate the 1D radial motion for one period.

    Uniform sampling in time along a periodic orbit is uniform sampling of the
    angle variable on the invariant torus (book section 17.2); this avoids the
    biased "sample r uniformly" shortcut the book explicitly warns against.

    Returns (r_grid, p_grid) of length n_table covering one period [0, T).
    """
    r_in, r_out = _radial_turning_points(e_int, j2, par, mu_r)
    w0 = par.beta * math.sqrt(2.0 * par.d_e / mu_r)
    t_guess = 2.0 * math.pi / w0

    def force(r):
        return -(_effective_potential(r + 1.0e-9, j2, par, mu_r) -
                 _effective_potential(r - 1.0e-9, j2, par, mu_r)) / 2.0e-9

    # find the period: start at the outer turning point, integrate until the
    # orbit first returns to r_out with positive radial momentum
    r, p = r_out * (1.0 - 1.0e-9), 0.0
    t = 0.0
    dt_o = t_guess / 20000.0
    left = False
    while not (left and r >= r_out and p > 0.0):
        if r < r_out - 1.0e-6:
            left = True
        r, p = _verlet_1d(r, p, dt_o, 1, force)
        t += dt_o
        if t > 50.0 * t_guess:
            raise RuntimeError("radial orbit did not close")
    n_tab = max(n_table, 64)
    dt_tab = t / n_tab
    r, p = r_out * (1.0 - 1.0e-9), 0.0
    rs = np.empty(n_tab)
    ps = np.empty(n_tab)
    for k in range(n_tab):
        rs[k] = r
        ps[k] = p
        r, p = _verlet_1d(r, p, dt_tab, 1, force)
    return rs, ps


@dataclass
class InitialSampler:
    """Pre-computed invariant-torus data shared by a whole ensemble.

    All trajectories of an ensemble share (v, j), so the radial orbit table
    and the EBK energy are computed once and reused (uniform-in-time sampling
    of the stored orbit then gives each trajectory its own phase).
    """
    cfg: QCTConfig
    par: LEPSParameters
    e_int: float = 0.0
    r_tab: np.ndarray = None
    p_tab: np.ndarray = None

    @classmethod
    def build(cls, cfg: QCTConfig, par: Optional[LEPSParameters] = None) -> "InitialSampler":
        if par is None:
            par = LEPSParameters()
        mu_r = 0.5 * cfg.mass
        if cfg.init_mode == "action":
            e_int = ebk_internal_energy(cfg.v_initial, cfg.j_initial, par, mu_r)
        elif cfg.init_mode == "energy":
            # H_int = G(v) + F_v(j), referred to the dissociation limit; the
            # dynamics work with absolute energies whose bottom sits at -D_e.
            e_int = (-par.d_e + morse_vibrational_energy(cfg.v_initial, par, mu_r) +
                     cfg.j_initial * (cfg.j_initial + 1) /
                     (2.0 * mu_r * par.r_e ** 2))
        else:
            raise ValueError(f"unknown init_mode: {cfg.init_mode}")
        j2 = float(cfg.j_initial * (cfg.j_initial + 1))
        r_tab, p_tab = _radial_orbit_table(e_int, j2, par, mu_r)
        return cls(cfg=cfg, par=par, e_int=e_int, r_tab=r_tab, p_tab=p_tab)

    def sample(self, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray, float]:
        """Return (q, p, b) for one trajectory: positions, momenta, impact parameter."""
        cfg = self.cfg
        m = cfg.mass
        mu_r = 0.5 * m
        mu_r_tot = 2.0 * m / 3.0
        j_mag = math.sqrt(cfg.j_initial * (cfg.j_initial + 1))

        # uniform time along the orbit = uniform vibrational angle
        idx = int(rng.integers(0, len(self.r_tab)))
        r0, p_r0 = float(self.r_tab[idx]), float(self.p_tab[idx])

        # rotation: random angular-momentum direction and orbital phase
        j_hat = _random_unit_vector(rng)
        j_vec = j_mag * j_hat
        u0 = _vector_orthogonal_to(j_hat, rng)
        phi_rot = 2.0 * math.pi * rng.random()
        c, s = math.cos(phi_rot), math.sin(phi_rot)
        axis = u0 * c + np.cross(j_hat, u0) * s
        axis = axis / np.linalg.norm(axis)          # stays orthogonal to j_hat

        # collision geometry: b uniform on the disk (b = b_max sqrt(xi))
        b = cfg.b_max * math.sqrt(rng.random())
        rx = math.sqrt(max(cfg.r_start ** 2 - b ** 2, 0.0))
        r_vec = np.array([rx, b, 0.0])
        speed = math.sqrt(2.0 * cfg.e_coll / mu_r_tot)
        r_dot = np.array([-speed, 0.0, 0.0])        # inward along x

        p_r_vec = p_r0 * axis + np.cross(j_vec, axis) / r0

        # Jacobi -> cartesian with total momentum zero (equal masses):
        #   v_A = (2/3) Rdot,  v_{B,C} = -(1/3) Rdot -/+ rdot/2
        v_a = (2.0 / 3.0) * r_dot
        v_rel = p_r_vec / mu_r
        v_b = -(1.0 / 3.0) * r_dot - 0.5 * v_rel
        v_c = -(1.0 / 3.0) * r_dot + 0.5 * v_rel
        q = np.stack([r_vec, -0.5 * axis * r0, 0.5 * axis * r0])
        p = m * np.stack([v_a, v_b, v_c])
        return q, p, b


def qct_init_trajectory(cfg: QCTConfig, rng: np.random.Generator,
                        par: Optional[LEPSParameters] = None,
                        sampler: Optional[InitialSampler] = None) -> QCTTrajectory:
    """Sample one initial condition on the invariant torus (book section 17.2)."""
    if par is None:
        par = LEPSParameters()
    if sampler is None:
        sampler = InitialSampler.build(cfg, par)
    q, p, b = sampler.sample(rng)
    traj = QCTTrajectory(q=q, p=p, b_impact=b)
    traj.e_total = _total_energy(q, p, cfg.mass, par)
    traj.j_total = _total_angular_momentum(q, p, cfg.mass)
    traj.r_dir_initial = np.array([-1.0, 0.0, 0.0])
    return traj


def _random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=3)
    n = np.linalg.norm(v)
    while n < 1.0e-12:
        v = rng.normal(size=3)
        n = np.linalg.norm(v)
    return v / n


def _vector_orthogonal_to(j_hat, rng):
    for cand in np.eye(3):
        c = cand - j_hat * np.dot(cand, j_hat)
        if np.linalg.norm(c) > 0.1:
            return c / np.linalg.norm(c)
    return _random_unit_vector(rng)


# ---------------------------------------------------------------------------
# Propagation: Velocity-Verlet (book section 17.3)
# ---------------------------------------------------------------------------

def _total_energy(q, p, m, par) -> float:
    kinetic = float(np.sum(p * p) / (2.0 * m))
    v, _ = leps_energy_gradient(q, par)
    return kinetic + float(v)


def _total_angular_momentum(q, p, m) -> np.ndarray:
    return np.cross(q, p).sum(axis=0)


def qct_propagate_step(q: np.ndarray, p: np.ndarray, dt: float, m: float,
                       par: Optional[LEPSParameters] = None,
                       force: Optional[np.ndarray] = None):
    """One Velocity-Verlet step in place; returns (p, f).

    p_{n+1/2} = p_n + dt/2 F(q_n);  q_{n+1} = q_n + dt p_{n+1/2}/m;
    p_{n+1}   = p_{n+1/2} + dt/2 F(q_{n+1}).

    ``force`` may carry F(q_n) from the previous step to reuse; the returned
    value is F(q_{n+1}), so a trajectory loop needs exactly ONE force
    evaluation per step (half the cost of a naive implementation).
    """
    if force is None:
        _, f = leps_energy_gradient(q, par)
    else:
        f = force
    p += 0.5 * dt * f
    q += dt * p / m
    _, f = leps_energy_gradient(q, par)
    p += 0.5 * dt * f
    return p, f


def _product_jacobi(q: np.ndarray):
    """Identify the product diatom (smallest-separation pair) and single atom."""
    r, _ = _pair_distances(q[None, ...])
    r = r[0]
    pairs = ((0, 1), (1, 2), (0, 2))
    seps = [r[a, b] for a, b in pairs]
    pair = pairs[int(np.argmin(seps))]
    atom = ({0, 1, 2} - set(pair)).pop()
    return pair, atom


def propagate_trajectory(traj: QCTTrajectory, cfg: QCTConfig,
                         par: Optional[LEPSParameters] = None) -> QCTTrajectory:
    """Propagate a single trajectory until the product separation exceeds r_end.

    Reuses the force between Velocity-Verlet half-kicks: one LEPS gradient
    evaluation per step.
    """
    if par is None:
        par = LEPSParameters()
    if cfg.r_end <= cfg.r_start:
        raise ValueError("cfg.r_end must exceed cfg.r_start")
    m = cfg.mass
    q = traj.q.copy()
    p = traj.p.copy()
    _, f = leps_energy_gradient(q, par)
    for _ in range(cfg.max_steps):
        p, f = qct_propagate_step(q, p, cfg.dt, m, par, force=f)
        traj.time += cfg.dt
        pair, atom = _product_jacobi(q)
        com = 0.5 * (q[pair[0]] + q[pair[1]])
        if float(np.linalg.norm(q[atom] - com)) > cfg.r_end:
            traj.converged = True
            break
    traj.q = q
    traj.p = p
    return traj


def _propagate_ensemble(q0: np.ndarray, p0: np.ndarray, cfg: QCTConfig,
                        par: LEPSParameters):
    """Vectorized Velocity-Verlet over all trajectories at once.

    Equivalent to calling :func:`propagate_trajectory` on each trajectory
    (same force field, same integrator); trajectories are frozen once their
    product separation exceeds r_end. Returns (q_f, p_f, times, converged).

    Hardware-efficiency notes: the three pair separations are computed
    directly (no (N, 3, 3, 3) distance tensor), the product-pair lookup is a
    table gather instead of a per-step Python loop, and the loop exits as
    soon as every trajectory has terminated.
    """
    if cfg.r_end <= cfg.r_start:
        raise ValueError("cfg.r_end must exceed cfg.r_start")
    m = cfg.mass
    n = q0.shape[0]
    q = q0.copy()
    p = p0.copy()
    q_f = q0.copy()
    p_f = p0.copy()
    active = np.ones(n, dtype=bool)
    times = np.zeros(n)
    converged = np.zeros(n, dtype=bool)
    pair_list = ((0, 1), (1, 2), (0, 2))
    n1_tab = np.array([pr[0] for pr in pair_list])
    n2_tab = np.array([pr[1] for pr in pair_list])
    nat_tab = np.array([({0, 1, 2} - set(pr)).pop() for pr in pair_list])
    rows = np.arange(n)
    _, f = leps_energy_gradient(q, par)
    for _ in range(cfg.max_steps):
        p += 0.5 * cfg.dt * f
        q += cfg.dt * p / m
        _, f = leps_energy_gradient(q, par)
        p += 0.5 * cfg.dt * f
        times[active] += cfg.dt
        # three pair separations, computed directly (memory-lean)
        d01 = q[:, 0, :] - q[:, 1, :]
        d12 = q[:, 1, :] - q[:, 2, :]
        d02 = q[:, 0, :] - q[:, 2, :]
        seps = np.stack([np.einsum("ij,ij->i", d01, d01),
                         np.einsum("ij,ij->i", d12, d12),
                         np.einsum("ij,ij->i", d02, d02)], axis=1)
        np.sqrt(seps, out=seps)
        pair_idx = np.argmin(seps, axis=1)
        n1 = n1_tab[pair_idx]
        n2 = n2_tab[pair_idx]
        natom = nat_tab[pair_idx]
        com = 0.5 * (q[rows, n1, :] + q[rows, n2, :])
        r_sep = np.linalg.norm(q[rows, natom, :] - com, axis=1)
        done = active & (r_sep > cfg.r_end)
        if done.any():
            # snapshot the terminal state; frozen trajectories must not be
            # integrated further, otherwise q keeps drifting with stale p
            q_f[done] = q[done]
            p_f[done] = p[done]
            converged |= done
            active &= ~done
            f[~active] = 0.0
            if not active.any():
                break
    still = ~converged
    q_f[still] = q[still]
    p_f[still] = p[still]
    return q_f, p_f, times, converged


# ---------------------------------------------------------------------------
# Final-state analysis (book section 17.4)
# ---------------------------------------------------------------------------

def _product_internal_state(q, p, pair, atom, m, par, r_dir_init=None):
    """Product Jacobi analysis: (e_int, |j'|, e_trans, cos(theta)).

    The internal energy is evaluated directly from the product diatom
    coordinates, e_int = |p_r|^2/(2 mu_r) + V_Morse(r): at R > r_end the LEPS
    surface has already reduced to the isolated-pair Morse, and this avoids
    the catastrophic cancellation of E_total - E_trans at high collision
    energies (the classic precision pitfall of QCT product analysis).
    """
    r_vec = q[pair[1]] - q[pair[0]]
    p_r_vec = 0.5 * (p[pair[1]] - p[pair[0]])            # mu_r = m/2
    com_p = 0.5 * (p[pair[0]] + p[pair[1]])
    mu_r_tot = 2.0 * m / 3.0
    p_big = mu_r_tot * (p[atom] / m - com_p / m)         # P = mu_R (v_atom - v_com)
    j_vec = np.cross(r_vec, p_r_vec)
    r_len = float(np.linalg.norm(r_vec))
    em = math.exp(-par.beta * (r_len - par.r_e))
    v_morse = par.d_e * (em * em - 2.0 * em)
    mu_r = 0.5 * m
    e_int = float(p_r_vec @ p_r_vec) / (2.0 * mu_r) + v_morse
    # translational energy of the relative motion: |P|^2/(2 mu_R); the
    # l^2/(2 mu R^2) term of the Jacobi Hamiltonian is the tangential part of
    # |P|^2 and must not be added again.
    e_trans = float(p_big @ p_big / (2.0 * mu_r_tot))
    if r_dir_init is None:
        r_dir_init = np.array([-1.0, 0.0, 0.0])
    pf = np.linalg.norm(p_big)
    cos_theta = float(np.dot(r_dir_init, p_big) / pf) if pf > 0.0 else 1.0
    return e_int, float(np.linalg.norm(j_vec)), e_trans, cos_theta


def qct_analyze_final_state(traj: QCTTrajectory, cfg: QCTConfig,
                            par: Optional[LEPSParameters] = None) -> QCTTrajectory:
    """Histogram binning of the final state (book section 17.4).

    j' = round(-1/2 + sqrt(1/4 + |j'|^2));  v' follows by removing the
    rotational energy F_v'(j') and inverting the Morse levels G(v). Trajectories
    whose internal energy falls below the vibrational zero-point energy are
    marked v' = -1 (ZPE leakage); under the passive ZPE constraint
    (cfg.enforce_zpe) they are discarded from state-resolved statistics.
    """
    if par is None:
        par = LEPSParameters()
    m = cfg.mass
    pair, atom = _product_jacobi(traj.q)
    traj.reactive = (tuple(sorted(pair)) != (1, 2))
    e_int, j_mag, _, _ = _product_internal_state(traj.q, traj.p, pair, atom, m, par,
                                                 traj.r_dir_initial)
    traj.e_internal = e_int
    traj.j_final = int(round(-0.5 + math.sqrt(0.25 + j_mag * j_mag)))
    mu_r = 0.5 * m
    rot = traj.j_final * (traj.j_final + 1) / (2.0 * mu_r * par.r_e ** 2)
    # e_int is absolute (bottom at -D_e); refer it to the dissociation limit
    g = e_int + par.d_e - rot
    w0 = par.beta * math.sqrt(2.0 * par.d_e / mu_r)
    c = par.beta ** 2 / (2.0 * mu_r)
    disc = w0 * w0 - 4.0 * c * g
    if g < morse_vibrational_energy(0.0, par, mu_r) or disc < 0.0:
        traj.v_final = -1        # below ZPE: not assignable
    else:
        w = (w0 - math.sqrt(disc)) / (2.0 * c)
        traj.v_final = int(round(w - 0.5))
    return traj


# ---------------------------------------------------------------------------
# Ensemble driver and cross sections (book sections 17.4-17.5)
# ---------------------------------------------------------------------------

def run_qct_ensemble(cfg: QCTConfig, par: Optional[LEPSParameters] = None,
                     progress: Optional[Callable[[int, int], None]] = None) -> QCTResult:
    """Run an ensemble of QCT trajectories and estimate cross sections.

    Initial conditions are sampled on the invariant torus of the selected
    (v, j) state; the ensemble is propagated with the vectorized Velocity-Verlet
    integrator (identical dynamics to :func:`propagate_trajectory`).
    """
    if par is None:
        par = LEPSParameters()
    rng = np.random.default_rng(cfg.seed)
    sampler = InitialSampler.build(cfg, par)
    trajs = []
    for _ in range(cfg.n_traj):
        trajs.append(qct_init_trajectory(cfg, rng, par, sampler))
    q0 = np.stack([t.q for t in trajs])
    p0 = np.stack([t.p for t in trajs])
    qf, pf, times, converged = _propagate_ensemble(q0, p0, cfg, par)
    for i, t in enumerate(trajs):
        t.q, t.p, t.time, t.converged = qf[i], pf[i], float(times[i]), bool(converged[i])
        qct_analyze_final_state(t, cfg, par)
        if progress is not None:
            progress(i + 1, cfg.n_traj)
    return _assemble_result(trajs, cfg)


def _assemble_result(trajs, cfg) -> QCTResult:
    n = len(trajs)
    n_rxn = sum(1 for t in trajs if t.reactive)
    res = QCTResult(n_traj=n, n_reactive=n_rxn)
    prob = n_rxn / n if n else 0.0
    res.cross_section = math.pi * cfg.b_max ** 2 * prob
    res.stat_error = (math.pi * cfg.b_max ** 2 *
                      math.sqrt(max(prob * (1.0 - prob), 0.0) / n) if n else 0.0)
    lo, hi = wilson_interval(n_rxn, n)
    res.wilson_lo = math.pi * cfg.b_max ** 2 * lo
    res.wilson_hi = math.pi * cfg.b_max ** 2 * hi
    res.trajectories = trajs
    res.opacity_b, res.opacity_p = opacity_function(trajs, cfg)
    res.state_cross_section = state_resolved_cross_sections(trajs, cfg)
    return res


def opacity_function(trajs, cfg, n_bins: Optional[int] = None):
    """Opacity function P_rxn(b) stratified into collision-parameter bins."""
    n_bins = n_bins or cfg.n_b_bins
    edges = np.linspace(0.0, cfg.b_max, n_bins + 1)
    center = 0.5 * (edges[:-1] + edges[1:])
    p_rxn = np.zeros(n_bins)
    for s in range(n_bins):
        sel = [t for t in trajs if edges[s] <= t.b_impact < edges[s + 1]]
        if sel:
            p_rxn[s] = sum(1 for t in sel if t.reactive) / len(sel)
    return center, p_rxn


def stratified_cross_section(trajs, cfg) -> float:
    """sigma = 2 pi sum_s Delta b_s b_s P_s (book section 17.4)."""
    edges = np.linspace(0.0, cfg.b_max, cfg.n_b_bins + 1)
    total = 0.0
    for s in range(cfg.n_b_bins):
        db = edges[s + 1] - edges[s]
        b_s = 0.5 * (edges[s] + edges[s + 1])
        sel = [t for t in trajs if edges[s] <= t.b_impact < edges[s + 1]]
        if sel:
            p_s = sum(1 for t in sel if t.reactive) / len(sel)
            total += db * b_s * p_s
    return 2.0 * math.pi * total


def qct_cross_section(n_rxn: int, n_tot: int, b_max: float):
    """Monte Carlo estimator sigma = pi b_max^2 N_rxn / N_tot with binomial error."""
    if n_tot <= 0:
        return 0.0, 0.0
    p = n_rxn / n_tot
    sigma = math.pi * b_max ** 2 * p
    err = math.pi * b_max ** 2 * math.sqrt(max(p * (1.0 - p), 0.0) / n_tot)
    return sigma, err


def state_resolved_cross_sections(trajs, cfg, enforce_zpe: Optional[bool] = None):
    """sigma_{v'j'} = pi b_max^2 N_{v'j'} / N_tot (book section 17.4)."""
    enforce_zpe = cfg.enforce_zpe if enforce_zpe is None else enforce_zpe
    pool = [t for t in trajs if not (t.reactive and t.v_final < 0)] if enforce_zpe else trajs
    n_tot = len(pool)
    out: Dict[Tuple[int, int], float] = {}
    if n_tot == 0:
        return out
    for t in pool:
        if t.reactive and t.v_final >= 0:
            key = (t.v_final, t.j_final)
            out[key] = out.get(key, 0.0) + 1.0
    for key in out:
        out[key] *= math.pi * cfg.b_max ** 2 / n_tot
    return out


def differential_cross_section(trajs, cfg, n_theta: int = 18,
                               b_max: Optional[float] = None,
                               par: Optional[LEPSParameters] = None):
    """Monte Carlo d(sigma)/d(Omega) histogram over the scattering angle.

    Uses d sigma/d Omega = (pi b_max^2 / N) * N_bin / (2 pi sin(theta) dtheta):
    the sin(theta) weight converts the histogram into a solid-angle density,
    and the multi-branch structure of theta(b) (rainbow scattering) is summed
    automatically by the binning (book section 17.4). Bins touching
    theta = 0 or pi are returned as NaN: the forward/backward directions need
    separate treatment, as the book notes. ``par`` should match the surface
    the ensemble was run with.
    """
    b_max = cfg.b_max if b_max is None else b_max
    edges = np.linspace(0.0, math.pi, n_theta + 1)
    dtheta = edges[1] - edges[0]
    theta_c = 0.5 * (edges[:-1] + edges[1:])
    counts = np.zeros(n_theta)
    n = len(trajs)
    for t in trajs:
        if not t.reactive:
            continue
        pair, atom = _product_jacobi(t.q)
        _, _, _, cos_t = _product_internal_state(t.q, t.p, pair, atom, cfg.mass,
                                                 par, t.r_dir_initial)
        ang = math.acos(max(-1.0, min(1.0, cos_t)))
        counts[min(int(ang / dtheta), n_theta - 1)] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        dsdo = (math.pi * b_max ** 2 / max(n, 1) * counts /
                (2.0 * math.pi * np.sin(theta_c) * dtheta))
    dsdo[(theta_c < dtheta) | (theta_c > math.pi - dtheta)] = np.nan
    return theta_c, dsdo


def qct_thermal_rate(energies, sigmas, temperature: float, mu_r: float):
    """Thermal rate coefficient k(T) from sigma(E) (book section 17.4):

        k = [8 / (pi mu (kT)^3)]^(1/2) * int E sigma(E) exp(-E/kT) dE.

    ``energies`` (Hartree) and ``sigmas`` (bohr^2) are sampled on a grid; the
    result is in a.u. of rate (bohr^3 per a.u. time), k_B = 1 in a.u.
    """
    e = np.asarray(energies, dtype=float)
    s = np.asarray(sigmas, dtype=float)
    order = np.argsort(e)
    e, s = e[order], s[order]
    integrand = e * s * np.exp(-e / temperature)
    integral = np.trapz(integrand, e)
    pref = math.sqrt(8.0 / (math.pi * mu_r * temperature ** 3))
    return pref * integral


def thermal_population(states, temperature: float):
    """Boltzmann population p_vj = (2j+1) exp(-E_vj/T) / Z over given states.

    ``states`` is a sequence of (v, j, e_vj) tuples with energies in Hartree.
    """
    e = np.array([s[2] for s in states], dtype=float)
    deg = np.array([2 * s[1] + 1 for s in states], dtype=float)
    w = deg * np.exp(-(e - e.min()) / temperature)
    return w / w.sum()


def wilson_interval(n_success: int, n_total: int, z: float = 1.959963984540054):
    """Wilson score interval for a binomial proportion (book section 17.5).

    Recommended over the plain normal approximation for low-probability events.
    """
    if n_total == 0:
        return 0.0, 1.0
    p = n_success / n_total
    z2 = z * z
    denom = 1.0 + z2 / n_total
    center = (p + z2 / (2.0 * n_total)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n_total + z2 / (4.0 * n_total ** 2)) / denom
    return max(center - half, 0.0), min(center + half, 1.0)
