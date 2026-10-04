"""
Smolyak sparse grids with nested Clenshaw-Curtis rules.

Python reference implementation of book section 18.4 and the proposed
`mod_sparse_grid` (book section 18.8):

    U_N^D f = sum_{N-D+1 <= |l|_1 <= N} (-1)^{N-|l|_1} C(D-1, N-|l|_1)
              [ x_ kappa U_{l_kappa}^(kappa) ] f,     l_kappa >= 1,

equivalently the difference form sum_{|l|_1 <= N} Delta_{l} with
Delta_l = U_l - U_{l-1}. For smooth integrands the grid grows like
O(N (log N)^(D-1)) instead of O(n^D); the growth depends on the smoothness
of f and on the 1D rule (book section 18.4 caveat).

Fortran counterpart (roadmap): `mod_sparse_grid` with `sparse_grid_t`,
`smolyak_build`, `sparse_grid_refine`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import List, Optional, Tuple

import numpy as np

__all__ = ["SparseGrid", "clenshaw_curtis", "smolyak_build",
           "smolyak_integrate", "smolyak_refine", "smolyak_multi_indices"]


@lru_cache(maxsize=128)
def _clenshaw_curtis_cached(n: int):
    """Cached CC computation; arrays are returned as immutable tuples."""
    k = np.arange(n)
    theta = math.pi * k / (n - 1)
    x = np.cos(theta)
    tj = np.cos(np.outer(np.arange(n), theta))          # (n, n): T_j(x_k)
    moments = np.zeros(n)
    moments[0] = 2.0
    for j in range(2, n, 2):
        moments[j] = 2.0 / (1.0 - j * j)
    w = np.linalg.solve(tj, moments)
    return tuple(x), tuple(w)


def clenshaw_curtis(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Nested Clenshaw-Curtis rule with n points on [-1, 1] (cached).

    Points x_k = cos(pi k / (n-1)) are nested across dyadic refinements
    (n = 2^(l-1) + 1 for level l >= 2, n = 1 for level 1), which is what makes
    the Smolyak grid economical. Weights follow from matching the moments of
    the Chebyshev polynomials T_j, j = 0..n-1.
    """
    n = int(n)
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return np.zeros(1), np.array([2.0])
    x, w = _clenshaw_curtis_cached(n)
    return np.asarray(x), np.asarray(w)


def _rule_size(level: int) -> int:
    """Number of 1D points at a Smolyak level (l = 1 -> 1 point)."""
    return 1 if level <= 1 else 2 ** (level - 1) + 1


def smolyak_multi_indices(d: int, n: int) -> List[Tuple[int, ...]]:
    """Multi-indices with l_k >= 1 and N-D+1 <= |l|_1 <= N (book section 18.4)."""
    out: List[Tuple[int, ...]] = []

    def rec(prefix, dims_left, budget_left):
        if dims_left == 0:
            s = sum(prefix)
            if n - d + 1 <= s <= n:
                out.append(tuple(prefix))
            return
        for val in range(1, budget_left - (dims_left - 1) + 1):
            rec(prefix + [val], dims_left - 1, budget_left - val)

    rec([], d, n)
    return out


@dataclass
class SparseGrid:
    """Smolyak sparse grid: points (n_pts, D) and quadrature weights."""
    d: int
    level: int
    points: np.ndarray
    weights: np.ndarray
    multi_indices: List[Tuple[int, ...]]

    def integrate(self, f_values: np.ndarray) -> float:
        """Quadrature of f sampled on grid.points (book section 18.4)."""
        return float(np.dot(self.weights, np.asarray(f_values, dtype=float)))


def smolyak_build(d: int, n: int) -> SparseGrid:
    """Build the Smolyak sparse grid U_N^D with nested Clenshaw-Curtis rules.

    Coefficients follow the book formula; coinciding points from nested
    1D rules are merged by summing their contributions.
    """
    if d < 1 or n < d:
        raise ValueError("need d >= 1 and N >= D")
    multi = smolyak_multi_indices(d, n)
    # cache 1D rules
    rules = {}
    for l in range(1, n + 1):
        m = _rule_size(l)
        if m not in rules:
            rules[m] = clenshaw_curtis(m)
    acc: dict = {}
    for l_multi in multi:
        coeff = (-1.0) ** (n - sum(l_multi)) * math.comb(d - 1, n - sum(l_multi))
        grids = [rules[_rule_size(l)] for l in l_multi]
        # tensor-product points
        shape = tuple(g[0].size for g in grids)
        flat_idx = np.indices(shape).reshape(d, -1).T
        xs = np.stack([grids[k][0][flat_idx[:, k]] for k in range(d)], axis=1)
        ws = np.prod(np.stack([grids[k][1][flat_idx[:, k]] for k in range(d)], axis=1), axis=1)
        for point, w in zip(xs, ws):
            key = tuple(np.round(point, 12))
            if key in acc:
                acc[key][1] += coeff * float(w)
            else:
                acc[key] = [point.copy(), coeff * float(w)]
    points = np.stack([v[0] for v in acc.values()])
    weights = np.array([v[1] for v in acc.values()])
    return SparseGrid(d=d, level=n, points=points, weights=weights,
                      multi_indices=multi)


def smolyak_integrate(grid: SparseGrid, f_values: np.ndarray) -> float:
    """Functional form of SparseGrid.integrate."""
    return grid.integrate(f_values)


def smolyak_refine(grid: SparseGrid) -> SparseGrid:
    """Refine to level N+1 (adaptive refinement driver, book section 18.4).

    Returns the refined grid; the difference between the two grids identifies
    the newly added points, whose contributions serve as hierarchical surpluses
    for adaptive refinement strategies.
    """
    return smolyak_build(grid.d, grid.level + 1)


def _delta_rule(level: int, rules: dict):
    """1D difference rule Delta_l = U_l - U_{l-1} on the U_l point set.

    Nested CC points make this a pure weight update: points of U_{l-1} keep
    w_l - w_{l-1}, new points keep w_l. U_0 = 0 gives Delta_1 = U_1.
    """
    x, w = rules[_rule_size(level)]
    if level == 1:
        return x, w.copy()
    x_prev, w_prev = rules[_rule_size(level - 1)]
    wmap = {round(float(xi), 12): float(wi) for xi, wi in zip(x_prev, w_prev)}
    dw = np.array([float(w[i]) - wmap.get(round(float(xi), 12), 0.0)
                   for i, xi in enumerate(x)])
    return x, dw


def _difference_form_grid(d: int, n: int) -> SparseGrid:
    """Reference construction via Delta_l (book section 18.4).

    U_N = sum_{|l|_1 <= N} Delta_{l1} x ... x Delta_{lD} with U_0 = 0. Used by
    the test suite to verify the coefficient form: both must produce identical
    points and weights.
    """
    from itertools import product
    rules = {}
    for l in range(1, n + 1):
        m = _rule_size(l)
        if m not in rules:
            rules[m] = clenshaw_curtis(m)
    delta_rules = [_delta_rule(l, rules) for l in range(1, n + 1)]
    acc: dict = {}

    def rec(prefix, dims_left, budget_left):
        if dims_left == 0:
            grids = [delta_rules[l - 1] for l in prefix]
            shape = [len(g[0]) for g in grids]
            for combo in product(*(range(s) for s in shape)):
                point = np.array([grids[k][0][combo[k]] for k in range(len(grids))])
                w = float(np.prod([grids[k][1][combo[k]] for k in range(len(grids))]))
                key = tuple(np.round(point, 12))
                if key in acc:
                    acc[key][1] += w
                else:
                    acc[key] = [point.copy(), w]
            return
        for val in range(1, budget_left - (dims_left - 1) + 1):
            rec(prefix + [val], dims_left - 1, budget_left - val)

    rec([], d, n)
    points = np.stack([v[0] for v in acc.values()])
    weights = np.array([v[1] for v in acc.values()])
    return SparseGrid(d=d, level=n, points=points, weights=weights, multi_indices=[])
