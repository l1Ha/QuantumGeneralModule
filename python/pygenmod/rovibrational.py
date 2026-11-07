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

def rot_matrix_cos2_theta(j: int, j_prime: int, m: int = 0) -> float:
    """Compute rigid rotor Raman matrix element <j, m | cos^2(theta) | j_prime, m>.

    Selection rules: Delta j = 0, +/-2 and Delta m = 0 (the angular part of the
    polarizability interaction in the book, section 2.4). The element is
    computed exactly through the completeness product

        <j m|cos^2|j' m> = sum_{j''} <j m|cos|j'' m> <j'' m|cos|j' m>,

    where cos(theta) connects only j'' = j +/- 1 and j'' = j' +/- 1, so the
    sum has at most two terms and reproduces the closed forms exactly
    (e.g. <2,0|cos^2|2,0> = 4/15 + 9/35 = 0.5238095238, <0,0|cos^2|0,0> = 1/3).
    """
    if abs(j - j_prime) not in (0, 2) or abs(m) > min(j, j_prime):
        return 0.0
    total = 0.0
    for j_mid in (j - 1, j + 1):
        if j_mid < abs(m):
            continue
        total += rot_matrix_cos_theta(j, j_mid, m) * rot_matrix_cos_theta(j_prime, j_mid, m)
    return total

ROT_FACTOR_AMU_A2_CM1 = 16.857631  # hbar^2/(2*amu*angstrom^2*hc) in cm^-1


def watson_asymmetric_levels(a_const: float, b_const: float, c_const: float,
                             j_qnum: int) -> np.ndarray:
    """Asymmetric-top rigid rotational levels (Watson limit, book Sec 3.4).

    H = A Ja^2 + B Jb^2 + C Jc^2 in the symmetric-top basis |J K M> (Wang
    matrix): diagonal (B+C)/2*J(J+1) + (A-(B+C)/2)*K^2, off-diagonal
    (B-C)/4*sqrt([J(J+1)-K(K+1)][J(J+1)-(K+1)(K+2)]) for Delta K = +2.
    Returns the 2J+1 ascending levels (same units as the input constants).
    Limits: B=C -> prolate B J(J+1) + (A-B)K^2; A=B -> oblate
    A J(J+1) - (A-C)K_c^2; A=B=C -> (2J+1)-fold degenerate A J(J+1).
    """
    n_dim = 2 * j_qnum + 1
    jj = float(j_qnum * (j_qnum + 1))
    h_mat = np.zeros((n_dim, n_dim))
    for k1 in range(-j_qnum, j_qnum + 1):
        i1 = k1 + j_qnum
        h_mat[i1, i1] = 0.5 * (b_const + c_const) * jj \
            + (a_const - 0.5 * (b_const + c_const)) * float(k1 * k1)
        if k1 + 2 <= j_qnum:
            k2 = k1 + 2
            i2 = k2 + j_qnum
            t1 = jj - float(k1 * (k1 + 1))
            t2 = jj - float((k1 + 1) * (k1 + 2))
            h_val = 0.25 * (b_const - c_const) * np.sqrt(max(t1 * t2, 0.0))
            h_mat[i2, i1] = h_val
            h_mat[i1, i2] = h_val
    return np.linalg.eigvalsh(h_mat)


def watson_asymmetric_levels_cd(a_const: float, b_const: float, c_const: float,
                                d_j: float, d_jk: float, d_k: float,
                                j_qnum: int, d1: float = 0.0,
                                d2: float = 0.0) -> np.ndarray:
    """Asymmetric-top levels including the full Watson A-reduced quartics.

        H_cd = -D_J J^4 - D_JK J^2 Ja^2 - D_K Ja^4
               + d1 J^2 (J_+^2 + J_-^2) + d2 (J_+^4 + J_-^4)

    D > 0 lowers the levels per the standard convention
    E = B J(J+1) - D_J J^2 (J+1)^2. The optional d1/d2 are the
    Delta K = +/-2/4 off-diagonal ladder terms (multi-atomic step 8):
        <K+2|d1 J^2 J_+^2|K> = d1 [J(J+1)]^2
            * sqrt([J(J+1)-K(K+1)][J(J+1)-(K+1)(K+2)])
        <K+4|d2 J_+^4|K> = d2 sqrt(prod_{m=0..3} [J(J+1)-(K+m)(K+m+1)])
    Trace theorem: the d1/d2 terms are purely off-diagonal (Tr = 0), so
    sum(levels) is d1/d2-independent to all orders (Watson A/S reduction
    decouples under the sum rule). The SPFIT (2K+/-3) symmetrized form is
    equivalent in the symmetrized-function basis; its sign convention is an
    external benchmark.
    Fortran mirror: `asymmetric_top_levels_cd` (src/mod_watson_rotor.f90).
    """
    n_dim = 2 * j_qnum + 1
    jj = float(j_qnum * (j_qnum + 1))
    h_mat = np.zeros((n_dim, n_dim))
    for k1 in range(-j_qnum, j_qnum + 1):
        i1 = k1 + j_qnum
        h_mat[i1, i1] = 0.5 * (b_const + c_const) * jj \
            + (a_const - 0.5 * (b_const + c_const)) * float(k1 * k1) \
            - d_j * jj * jj - d_jk * jj * float(k1 * k1) - d_k * float(k1 * k1) ** 2
        if k1 + 2 <= j_qnum:
            k2 = k1 + 2
            i2 = k2 + j_qnum
            t1 = jj - float(k1 * (k1 + 1))
            t2 = jj - float((k1 + 1) * (k1 + 2))
            h_val = 0.25 * (b_const - c_const) * np.sqrt(max(t1 * t2, 0.0)) \
                + d1 * jj * jj * np.sqrt(max(t1 * t2, 0.0))
            h_mat[i2, i1] = h_val
            h_mat[i1, i2] = h_val
        if d2 != 0.0 and k1 + 4 <= j_qnum:
            k4 = k1 + 4
            i4 = k4 + j_qnum
            prod4 = 1.0
            for m in range(4):
                km = k1 + m
                prod4 *= jj - float(km * (km + 1))
            h_val = d2 * np.sqrt(max(prod4, 0.0))
            h_mat[i4, i1] = h_val
            h_mat[i1, i4] = h_val
    return np.linalg.eigvalsh(h_mat)


def inertia_constants_from_geometry(positions: np.ndarray, masses: np.ndarray
                                    ) -> tuple[float, float, float, np.ndarray, np.ndarray]:
    """Principal moments and rotational constants A, B, C (cm^-1) from nuclear
    geometry (multi-atomic step 1, book Sec 3.4). Planar molecules satisfy the
    perpendicular-axis theorem I_c = I_a + I_b (usable as a self-check)."""
    positions = np.asarray(positions, dtype=float)
    masses = np.asarray(masses, dtype=float)
    if positions.shape != (3, len(masses)):
        raise ValueError("positions must have shape (3, N) with N = len(masses)")
    com = (masses[None, :] * positions).sum(axis=1) / masses.sum()
    r_com = positions - com[:, None]
    inertia = np.zeros((3, 3))
    for i in range(len(masses)):
        r2 = float(r_com[:, i] @ r_com[:, i])
        inertia += masses[i] * (r2 * np.eye(3) - np.outer(r_com[:, i], r_com[:, i]))
    moments = np.linalg.eigvalsh(inertia)
    if np.any(moments <= 0.0):
        raise ValueError("non-positive principal moment: degenerate or collinear geometry")
    constants = ROT_FACTOR_AMU_A2_CM1 / moments
    return float(constants[0]), float(constants[1]), float(constants[2]), \
        moments, com


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
