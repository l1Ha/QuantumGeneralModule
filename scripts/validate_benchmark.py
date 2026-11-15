#!/usr/bin/env python3
"""Validate a benchmark JSON record and optionally compare with a baseline.

Usage:
    python3 scripts/validate_benchmark.py RESULTS.json [BASELINE.json] [--factor F]

Schema + correctness-gate validation always runs (exit 1 on violation).
When BASELINE.json is given, per-kernel timings are compared and any kernel
slower than ``F`` x baseline (default 1.5) is reported; ``--strict`` turns
regressions into a nonzero exit (off by default: shared CI runners are too
noisy for hard timing gates, so the comparison is advisory).

Records are produced by ``make bench`` (benchmarks/bench_kernels.f90).
"""
import json
import sys

KERNELS = (
    "fft_1d_ms", "fft_2d_ms", "diag_symmetric_ms", "diag_hermitian_ms",
    "krylov_100steps_ms", "chebyshev_100steps_ms", "spline3d_1e6_ms",
)


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def validate(doc):
    errs = []
    if doc.get("record_type") != "benchmark":
        errs.append("record_type must be 'benchmark'")
    if doc.get("module") != "GeneralModule":
        errs.append("module must be 'GeneralModule'")
    if doc.get("all_correctness_checks_passed") is not True:
        errs.append("correctness gate failed (all_correctness_checks_passed != true)")
    res = doc.get("results")
    if not isinstance(res, dict):
        errs.append("results object missing")
        return errs
    for key in KERNELS:
        v = res.get(key)
        if not isinstance(v, (int, float)) or v <= 0:
            errs.append(f"results.{key} missing or non-positive")
    scaling = res.get("fft2d_scaling")
    if not isinstance(scaling, list) or not scaling:
        errs.append("results.fft2d_scaling missing")
    else:
        for entry in scaling:
            if not isinstance(entry.get("threads"), int) or entry.get("fft_2d_ms", 0) <= 0:
                errs.append("fft2d_scaling entries must carry threads/fft_2d_ms")
                break
    return errs


def compare(cur, base, factor):
    rows = []
    bres = base.get("results", {})
    for key in KERNELS:
        if key in bres and key in cur.get("results", {}):
            ratio = cur["results"][key] / bres[key]
            rows.append((key, bres[key], cur["results"][key], ratio))
    return rows


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    results_path = argv[1]
    baseline_path = None
    factor = 1.5
    strict = "--strict" in argv
    if "--factor" in argv:
        factor = float(argv[argv.index("--factor") + 1])
    for a in argv[2:]:
        if not a.startswith("--"):
            baseline_path = a
            break

    doc = load(results_path)
    errs = validate(doc)
    if errs:
        for e in errs:
            print(f"[FAIL] schema: {e}")
        return 1
    print(f"[PASS] schema + correctness gate ({results_path})")
    res = doc["results"]
    print("       kernels: " + ", ".join(f"{k.replace('_ms','')}={res[k]:.2f}ms" for k in KERNELS))
    scaling = res["fft2d_scaling"]
    t1 = scaling[0]["fft_2d_ms"]
    best = scaling[-1]
    print(f"       fft_2d scaling: {scaling[0]['threads']}t {t1:.2f}ms -> "
          f"{best['threads']}t {best['fft_2d_ms']:.2f}ms ({t1 / best['fft_2d_ms']:.2f}x)")

    rc = 0
    if baseline_path:
        rows = compare(doc, load(baseline_path), factor)
        print(f"       vs baseline {baseline_path} (factor {factor}x):")
        for key, b, c, ratio in rows:
            flag = "  <-- REGRESSION" if ratio > factor else ""
            print(f"         {key:26s} {b:9.3f} -> {c:9.3f} ms  x{ratio:5.2f}{flag}")
            if ratio > factor and strict:
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))
