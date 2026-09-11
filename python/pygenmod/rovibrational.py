"""
Rovibrational Dynamics and Matrix Generation Module
"""
import numpy as np

def rovibrational_state_index(v: int, j: int, j_max: int) -> int:
    """Map (v, J) quantum numbers to 0-based linear index."""
    return v * (j_max + 1) + j

def rovibrational_state_unindex(idx: int, j_max: int) -> tuple[int, int]:
    """Unindex linear index to (v, J)."""
    v = idx // (j_max + 1)
    j = idx % (j_max + 1)
    return v, j

def rot_matrix_cos_theta(j: int, j_prime: int, m: int = 0) -> float:
    """Compute rigid rotor dipole matrix element <j, m | cos(theta) | j_prime, m>."""
    if abs(j - j_prime) != 1 or abs(m) > min(j, j_prime):
        return 0.0
    j_max = max(j, j_prime)
    return np.sqrt((j_max**2 - m**2) / ((2.0 * j_max - 1.0) * (2.0 * j_max + 1.0)))

def calc_franck_condon_factors(chi_a: np.ndarray, chi_b: np.ndarray, dx: float) -> np.ndarray:
    """Compute Franck-Condon factor matrix FC(v, v') = |<chi_a,v | chi_b,v'>|^2."""
    n_va = chi_a.shape[1]
    n_vb = chi_b.shape[1]
    fc = np.zeros((n_va, n_vb))
    for va in range(n_va):
        for vb in range(n_vb):
            overlap = np.sum(chi_a[:, va] * chi_b[:, vb]) * dx
            fc[va, vb] = overlap ** 2
    return fc

def calc_rotational_constants_bv(chi: np.ndarray, r_grid: np.ndarray, dx: float, mass: float) -> np.ndarray:
    """Compute rotational constants B_v for each vibrational state."""
    nv = chi.shape[1]
    inv_r2 = 1.0 / (2.0 * mass * np.maximum(0.1, r_grid ** 2))
    b_v = np.zeros(nv)
    for v in range(nv):
        b_v[v] = np.sum((chi[:, v] ** 2) * inv_r2) * dx
    return b_v

def build_rovibrational_hamiltonian(v_max: int, j_max: int, e_vib: np.ndarray, b_v: np.ndarray) -> np.ndarray:
    """Build diagonal field-free rovibrational Hamiltonian E(v, J) = E_vib(v) + B_v * J * (J + 1)."""
    n_states = (v_max + 1) * (j_max + 1)
    h_diag = np.zeros(n_states)
    for v in range(v_max + 1):
        for j in range(j_max + 1):
            k = rovibrational_state_index(v, j, j_max)
            h_diag[k] = e_vib[v] + b_v[v] * j * (j + 1)
    return h_diag

def build_rovibrational_dipole_matrix(v_max: int, j_max: int, dip_vib: np.ndarray) -> np.ndarray:
    """Build full rovibrational transition dipole matrix <v, J | mu | v', J'>."""
    n_states = (v_max + 1) * (j_max + 1)
    dip_mat = np.zeros((n_states, n_states))
    for v in range(v_max + 1):
        for j in range(j_max + 1):
            k1 = rovibrational_state_index(v, j, j_max)
            for vp in range(v_max + 1):
                for jp in range(j_max + 1):
                    k2 = rovibrational_state_index(vp, jp, j_max)
                    if abs(j - jp) == 1:
                        dip_mat[k1, k2] = dip_vib[v, vp] * rot_matrix_cos_theta(j, jp, 0)
    return dip_mat
