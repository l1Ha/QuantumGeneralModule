"""
Discrete Variable Representation (DVR) and Fourier Grid Hamiltonian (FGH) solver.
Mirrors `mod_dvr_grid.f90`.
"""
import numpy as np
from dataclasses import dataclass
from .constants import PI


@dataclass
class SincDVR:
    x_min: float
    x_max: float
    n_points: int
    mass: float
    dx: float
    x: np.ndarray
    t_mat: np.ndarray


def dvr_sinc_init(x_min: float, x_max: float, n_points: int, mass: float = 1.0) -> SincDVR:
    """Initialize Colbert-Miller Sinc-DVR grid and analytical kinetic energy matrix."""
    dx = (x_max - x_min) / float(n_points + 1)
    x = x_min + np.arange(1, n_points + 1) * dx

    i_idx, j_idx = np.meshgrid(np.arange(1, n_points + 1), np.arange(1, n_points + 1), indexing='ij')
    diff = i_idx - j_idx

    factor = 1.0 / (2.0 * mass * dx**2)
    t_mat = np.zeros((n_points, n_points), dtype=float)

    # Diagonal
    diag_mask = (diff == 0)
    t_mat[diag_mask] = factor * (PI**2 / 3.0)

    # Off-diagonal
    off_diag = ~diag_mask
    sign = np.where(diff[off_diag] % 2 == 0, 1.0, -1.0)
    t_mat[off_diag] = factor * sign * (2.0 / (diff[off_diag]**2))

    return SincDVR(
        x_min=x_min,
        x_max=x_max,
        n_points=n_points,
        mass=mass,
        dx=dx,
        x=x,
        t_mat=t_mat
    )


def fgh_solve_bound_states(dvr: SincDVR, v_pot: np.ndarray):
    """
    Solve bound states of arbitrary 1D potential using Fourier Grid Hamiltonian (FGH).
    Returns (eigenvalues, normalized wavefunctions).
    Wavefunctions are normalized such that sum(|psi|^2) * dx = 1.
    """
    v_pot = np.asarray(v_pot, dtype=float)
    if len(v_pot) != dvr.n_points:
        raise ValueError(f"v_pot length ({len(v_pot)}) must match dvr n_points ({dvr.n_points})")

    h_mat = dvr.t_mat + np.diag(v_pot)
    eig_vals, eig_vecs = np.linalg.eigh(h_mat)

    # Normalize: psi(x_i) = eig_vecs(i, v) / sqrt(dx)
    wavefunctions = eig_vecs / np.sqrt(dvr.dx)

    return eig_vals, wavefunctions
