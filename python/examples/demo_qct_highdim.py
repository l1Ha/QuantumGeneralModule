#!/usr/bin/env python3
"""
GeneralModule Python Demonstration:
QCT reactive scattering (textbook Chapter 17) and high-dimensional quantum
methods (textbook Chapter 18) — SOP / tensor train / Smolyak / MCTDH.

Runtime: a few seconds on a laptop. Larger ensembles, bond dimensions and
grid levels should be run on a server or batch queue; scale n_traj, TT ranks
and Smolyak levels only after checking convergence at these settings.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

import pygenmod as pg


def section(title):
    print("=" * 66)
    print(f"  {title}")
    print("=" * 66)


def demo_qct():
    section("Chapter 17: QCT — H + H2 (LEPS) at E_coll = 1.5 eV")
    cfg = pg.QCTConfig(e_coll=1.5 * pg.EV2AU, b_max=3.0, r_start=9.0,
                       r_end=11.0, dt=3.0, n_traj=400, max_steps=8000,
                       v_initial=0, j_initial=0, seed=42)
    t0 = time.time()
    res = pg.run_qct_ensemble(cfg)
    print(f"  trajectories: {res.n_traj}   elapsed: {time.time()-t0:.2f} s")
    print(f"  reactive:     {res.n_reactive}")
    print(f"  sigma_rxn = {res.cross_section:.3f} +/- {res.stat_error:.3f} a0^2"
          f"   (Wilson 95%: [{res.wilson_lo:.3f}, {res.wilson_hi:.3f}])")
    strat = pg.stratified_cross_section(res.trajectories, cfg)
    print(f"  stratified estimate = {strat:.3f} a0^2")
    top = sorted(res.state_cross_section.items(), key=lambda kv: -kv[1])[:3]
    print("  top state-resolved channels:",
          {f"v'={v} j'={j}": round(s, 3) for (v, j), s in top})
    b, p = res.opacity_b, res.opacity_p
    nz = [(float(bi), float(pi)) for bi, pi in zip(b, p) if pi > 0]
    print("  opacity P(b) > 0 at b =",
          [(round(bi, 2), round(pi, 2)) for bi, pi in nz[:4]], "...")


def demo_highdim():
    section("Chapter 18: SOP / tensor train / Smolyak / MCTDH")
    rng = np.random.default_rng(7)
    n = 4
    terms = [pg.SOPTerm(coeff=float(rng.normal()),
                        factors=[rng.normal(size=(n, n)) for _ in range(3)])
             for _ in range(4)]
    sop = pg.SOPHamiltonian(dims=(n, n, n), terms=terms)
    v = rng.normal(size=(n, n, n))
    hv = pg.sop_apply(sop, v)
    print(f"  SOP action on a {n}^3 tensor: "
          f"output norm {np.linalg.norm(hv):.3f} (matrix-free)")

    pot = rng.normal(size=(5, 5, 5))
    potfit = pg.sop_from_potfit(pot, eps=1e-9)
    print(f"  POTFIT decomposition of a 5^3 potential: "
          f"{len(potfit.terms)} product terms")

    a = rng.normal(size=(4, 5, 6))
    tt = pg.tt_round(pg.tt_from_dense(a, eps=1e-10), eps=1e-10)
    rel = np.linalg.norm(pg.tt_to_dense(tt) - a) / np.linalg.norm(a)
    print(f"  TT-SVD + rounding: ranks {tt.ranks}, relative error {rel:.2e}")

    grid = pg.smolyak_build(3, 6)
    cs = np.array([1.0, 2.0, 3.0])
    exact = np.prod([2.0 * np.sin(c) / c for c in cs])
    val = grid.integrate(np.prod(np.cos(grid.points * cs), axis=-1))
    print(f"  Smolyak U_6^3: {len(grid.points)} points "
          f"(full tensor grid at 9 points/dim: 729), "
          f"integral error {abs(val - exact):.2e}")

    h1 = rng.normal(size=(5, 5)); h1 += h1.T
    h2 = rng.normal(size=(5, 5)); h2 += h2.T
    sop_m = pg.SOPHamiltonian(dims=(5, 5, 5), terms=[
        pg.SOPTerm(1.0, [h1, np.eye(5), np.eye(5)]),
        pg.SOPTerm(1.0, [np.eye(5), h2, np.eye(5)]),
        pg.SOPTerm(0.1, [h1, h2, np.eye(5)]),
    ])
    a0 = rng.normal(size=(5, 5, 5)) + 1j * rng.normal(size=(5, 5, 5))
    a0 /= np.linalg.norm(a0)
    st = pg.MCTDHState(A=a0.copy(),
                       spf=[np.eye(5, dtype=complex)] * 3)
    pg.mctdh_propagate(st, sop_m, pg.MCTDHConfig(dt=0.01, n_steps=100))
    print(f"  MCTDH (full SPF space, 100 RK4 steps): "
          f"norm = {st.norm():.12f}  (exact evolution reproduced)")

    spfs = []
    rr = np.random.default_rng(3)
    for _ in range(3):
        q, _ = np.linalg.qr(rr.normal(size=(5, 2)) + 1j * rr.normal(size=(5, 2)))
        spfs.append(q.T.copy())
    a_red = a0.copy()
    for k in range(3):
        a_red = np.moveaxis(np.tensordot(spfs[k].conj(), a_red,
                                         axes=([1], [k])), 0, k)
    a_red /= np.linalg.norm(a_red)
    st2 = pg.MCTDHState(A=a_red, spf=[s.copy() for s in spfs])
    e0 = float(np.real(np.sum(st2.A.conj() * _sop_apply_spf(st2, sop_m))))
    pg.mctdh_propagate(st2, sop_m, pg.MCTDHConfig(dt=0.01, n_steps=100))
    e1 = float(np.real(np.sum(st2.A.conj() * _sop_apply_spf(st2, sop_m))))
    print(f"  MCTDH (2 SPFs/mode): norm = {st2.norm():.10f}, "
          f"energy drift = {abs(e1-e0)/max(abs(e0), 1e-30):.2e}")


def _sop_apply_spf(state, sop):
    """H_A A in the SPF basis (eta contraction) — for the energy diagnostic."""
    import pygenmod.mctdh_core as M
    d = state.A.ndim
    out = np.zeros_like(state.A)
    for term in sop.terms:
        c = state.A
        for k in range(d):
            eta = M._spf_apply(state.spf[k], term.factors[k])
            c = M._apply_along_axis(eta, c, k)
        out += term.coeff * c
    return out


if __name__ == "__main__":
    demo_qct()
    print()
    demo_highdim()
    print("\nAll demos completed successfully.")
