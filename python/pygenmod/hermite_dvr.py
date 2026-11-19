"""
Hermite DVR (harmonic-oscillator basis) — Python mirror of
`mod_hermite_dvr` (src/mod_hermite_dvr.f90).

Physicists' Hermite functions
``psi_n(x) = (2^n n! sqrt(pi))^{-1/2} e^{-x^2/2} H_n(x)``
are the exact eigenbasis of ``H = p^2/(2m) + m omega^2 x^2 / 2`` (hbar = 1).
In mass-frequency scaled coordinates ``q = sqrt(m omega) x`` the kinetic
matrix follows in closed form from the ladder operators (no approximation),

    T_ij = (omega/2) <i|p_q^2|j>,
    <j|p_q^2|j> = j + 1/2,  <j-2|p_q^2|j> = -(1/2) sqrt(j (j-1)),

and potential matrix elements are evaluated by Gauss-Hermite quadrature
(Golub-Welsch nodes/weights), exact for polynomial potentials up to degree
2N-1. The quadrature must use the *de-exponentiated* function values
``psit_n(q) = psi_n(q) e^{+q^2/2}`` because the weights already carry
``e^{-q^2}``; using the raw ``psi`` would silently evaluate with an
``e^{-2 q^2}`` weight and corrupt every matrix element.
"""
import numpy as np
from math import lgamma

__all__ = [
    "gauss_hermite_nodes_weights",
    "hermite_function_values",
    "hermite_ho_kinetic",
    "hermite_dvr_hamiltonian",
]


def gauss_hermite_nodes_weights(n_ord: int):
    """Gauss-Hermite nodes/weights for the weight ``e^{-x^2}``.

    Nodes: Golub-Welsch — eigenvalues of the Jacobi matrix of the monic
    Hermite recurrence ``p_{k+1} = x p_k - (k/2) p_{k-1}`` (zero diagonal,
    off-diagonal ``sqrt(k/2)``). Replaces the earlier Newton root-finding,
    which could skip roots for poor initial guesses.

    Weights: Christoffel identity ``w_k = 1/(N psit_{N-1}(x_k)^2)`` with
    ``psit_n = psi_n e^{+x^2/2}`` (the de-exponentiated Hermite function).
    The textbook Golub-Welsch expression ``sqrt(pi) v_{1k}^2`` from
    ``numpy.linalg.eigh`` must NOT be used here: the first eigenvector
    component of an extreme node is ~1e-48, far below LAPACK's absolute
    noise floor (~1e-16), so the returned weight is polluted (observed
    1e-59 instead of 1e-96 for N=120) and the potential matrix picks up
    catastrophic cancellation. The Fortran mirror is immune because tql2
    accumulates rotations, preserving tiny components' relative accuracy.
    """
    jac = np.zeros((n_ord, n_ord))
    k = np.arange(1, n_ord)
    jac[k - 1, k] = np.sqrt(0.5 * k)
    jac[k, k - 1] = jac[k - 1, k]
    x_nodes = np.linalg.eigvalsh(jac)
    psi_last = np.array([hermite_function_values(x, n_ord - 1)[n_ord - 1]
                         for x in x_nodes])
    psit_last = psi_last * np.exp(0.5 * x_nodes ** 2)
    w_gauss = 1.0 / (n_ord * psit_last ** 2)
    return x_nodes, w_gauss


def hermite_function_values(x, n_max: int) -> np.ndarray:
    """Normalized Hermite functions psi_0..psi_{n_max} at point ``x``.

    Recurrence: ``H_{n+1} = 2 x H_n - 2 n H_{n-1}``, ``H_1' = 2`` (kept for
    interface parity with the Fortran poly+derivative kernel; the derivative
    is not used by the DVR assembly).
    """
    x = float(x)
    hn = np.zeros(n_max + 1)
    dhn = np.zeros(n_max + 1)
    hn[0] = 1.0
    if n_max >= 1:
        hn[1] = 2.0 * x
        dhn[1] = 2.0
    for k in range(1, n_max):
        hn[k + 1] = 2.0 * x * hn[k] - 2.0 * k * hn[k - 1]
        dhn[k + 1] = 2.0 * (k + 1) * hn[k]
    # normalization in log domain: ln norm = 0.5 (n ln2 + lgamma(n+1) + 0.5 ln pi).
    # The direct product 2^n n! overflows for n >~ 165 (2^159 159! ~ 1e329);
    # folding everything into a single exp with the e^{-x^2/2} factor also
    # saves one rounding step.
    n_idx = np.arange(n_max + 1)
    log_norm = 0.5 * (n_idx * np.log(2.0)
                      + np.array([lgamma(m + 1.0) for m in n_idx])
                      + 0.5 * np.log(np.pi))
    return hn * np.exp(-log_norm - 0.5 * x * x)


def hermite_ho_kinetic(n_ord: int, omega: float) -> np.ndarray:
    """HO kinetic matrix in the Hermite basis (ladder-operator closed form).

    ``p_q = i(a^dag - a)/sqrt(2)`` gives
    ``<j|p_q^2|j> = j + 1/2`` (0-based) and
    ``<j±2|p_q^2|j> = -(1/2) sqrt(...)``; ``T_ij = (omega/2) <i|p_q^2|j>``.
    """
    t_mat = np.zeros((n_ord, n_ord))
    j = np.arange(n_ord, dtype=float)          # 0-based quantum number
    diag = 0.5 * omega * (j + 0.5)
    t_mat[np.arange(n_ord), np.arange(n_ord)] = diag
    # <j-2|p_q^2|j> = -(1/2) sqrt(j (j-1)) for j >= 2 (0-based)
    if n_ord > 2:
        jm = np.arange(2, n_ord)
        val = -0.25 * omega * np.sqrt(jm * (jm - 1.0))
        t_mat[jm - 2, jm] = val
        t_mat[jm, jm - 2] = val
    return t_mat


def hermite_dvr_hamiltonian(n_ord: int, omega: float, v_at_nodes: np.ndarray):
    """Assemble the Hermite DVR Hamiltonian.

    Parameters
    ----------
    n_ord : int
        Basis dimension.
    omega : float
        Oscillator frequency; the caller folds the mass scaling into the
        potential (``V(x = q/sqrt(m omega))`` evaluated at the q-space nodes).
    v_at_nodes : ndarray (n_ord,)
        Potential values at the q-space Gauss-Hermite nodes.

    Returns
    -------
    h_mat : ndarray (n_ord, n_ord)
    x_nodes : ndarray (n_ord,)
    """
    v_at_nodes = np.asarray(v_at_nodes, dtype=float)
    x_nodes, w_gauss = gauss_hermite_nodes_weights(n_ord)
    # De-exponentiated values: psit(k, i) = psi_{i-1}(x_k) e^{+x_k^2/2}
    psit = np.empty((n_ord, n_ord))
    for k in range(n_ord):
        psi_k = hermite_function_values(x_nodes[k], n_ord - 1)
        psit[k, :] = psi_k * np.exp(0.5 * x_nodes[k] ** 2)
    # V_ij = sum_k psit[k,i] (w_k V_k) psit[k,j] = (psit.T * w*V) @ psit
    # (求和指标是 psit 的行 -> 必须转置; 曾因漏转置得到错误矩阵)
    v_mat = (psit.T * (w_gauss * v_at_nodes)) @ psit
    h_mat = hermite_ho_kinetic(n_ord, omega) + v_mat
    return h_mat, x_nodes
