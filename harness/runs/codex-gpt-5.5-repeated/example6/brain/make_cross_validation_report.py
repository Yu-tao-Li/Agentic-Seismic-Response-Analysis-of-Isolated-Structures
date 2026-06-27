"""Create machine-readable cross-validation report for Example 6."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def arr(data: Any) -> np.ndarray:
    return np.asarray(data, dtype=float)


def rel_diff(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    denom = max(1.0e-30, abs(a), abs(b))
    return abs(a - b) / denom


def nrmse(a: Any, b: Any) -> float | None:
    if a is None or b is None:
        return None
    aa = arr(a)
    bb = arr(b)
    if aa.shape != bb.shape:
        return None
    mask = np.isfinite(aa) & np.isfinite(bb)
    if not np.any(mask):
        return None
    diff = aa[mask] - bb[mask]
    ref = aa[mask]
    denom = float(np.max(ref) - np.min(ref))
    if denom < 1.0e-30:
        denom = max(1.0, float(np.max(np.abs(ref))))
    return float(math.sqrt(float(np.mean(diff**2))) / denom)


def loop_area(result: dict[str, Any]) -> float | None:
    values = result.get("hysteresis_loop_area")
    if values is None:
        return None
    return float(np.sum(np.asarray(values, dtype=float)))


def convergence_summary(result: dict[str, Any]) -> dict[str, Any] | None:
    conv = result.get("convergence")
    if conv is None:
        return None
    out: dict[str, Any] = {}
    for key in [
        "relative_residual_tolerance",
        "test_type",
        "tolerance",
        "max_iterations_per_step",
        "failed_step_count",
        "fallback_step_count",
        "max_relative_residual",
        "max_iterations_observed",
        "reason",
    ]:
        if key in conv:
            out[key] = conv[key]
    history = conv.get("iteration_history")
    if isinstance(history, list) and history:
        hist = np.asarray(history, dtype=float)
        out["iteration_history_summary"] = {
            "count": int(hist.size),
            "min": float(np.nanmin(hist)),
            "max": float(np.nanmax(hist)),
            "mean": float(np.nanmean(hist)),
        }
    residual_history = conv.get("relative_residual_history")
    if isinstance(residual_history, list) and residual_history:
        hist = np.asarray(residual_history, dtype=float)
        out["relative_residual_history_summary"] = {
            "count": int(hist.size),
            "max": float(np.nanmax(hist)),
            "mean": float(np.nanmean(hist)),
        }
    return out


def peak(result: dict[str, Any], key: str) -> float | None:
    peaks = result.get("peak_responses")
    if peaks is None:
        return None
    value = peaks.get(key)
    if value is None:
        return None
    return float(value)


def scalar_metrics(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    keys = {
        "fundamental_period": (a.get("fundamental_period"), b.get("fundamental_period")),
        "peak_roof_displacement": (peak(a, "roof_displacement_abs_max"), peak(b, "roof_displacement_abs_max")),
        "maximum_interstory_drift_ratio": (peak(a, "max_interstory_drift_ratio_abs"), peak(b, "max_interstory_drift_ratio_abs")),
        "peak_isolation_displacement": (peak(a, "isolation_displacement_abs_max"), peak(b, "isolation_displacement_abs_max")),
        "peak_total_base_shear": (peak(a, "total_base_shear_abs_max"), peak(b, "total_base_shear_abs_max")),
        "total_isolator_hysteresis_loop_area": (loop_area(a), loop_area(b)),
    }
    out: dict[str, Any] = {}
    for name, (av, bv) in keys.items():
        out[name] = {"a": av, "b": bv, "relative_difference": rel_diff(av, bv)}
    return out


def compare_pair(name_a: str, a: dict[str, Any], name_b: str, b: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    comparable = a.get("status") == "OK" and b.get("status") == "OK"
    pair: dict[str, Any] = {
        "software_a": name_a,
        "software_b": name_b,
        "physically_comparable": comparable,
        "status_a": a.get("status"),
        "status_b": b.get("status"),
        "scalar_metrics": scalar_metrics(a, b),
        "curve_metrics": {
            "roof_displacement_nrmse": nrmse(a.get("roof_displacement"), b.get("roof_displacement")),
            "isolation_displacement_nrmse": nrmse(a.get("isolation_displacement"), b.get("isolation_displacement")),
            "base_shear_nrmse": nrmse(a.get("total_base_shear"), b.get("total_base_shear")),
        },
        "convergence": {
            "a": convergence_summary(a),
            "b": convergence_summary(b),
        },
    }
    if not comparable:
        pair["pass"] = False
        pair["reason"] = "At least one generated software result is missing required physical response output."
        return pair

    checks: dict[str, bool] = {}
    checks["fundamental_period"] = pair["scalar_metrics"]["fundamental_period"]["relative_difference"] <= thresholds["fundamental_period_relative"]
    checks["peak_roof_displacement"] = pair["scalar_metrics"]["peak_roof_displacement"]["relative_difference"] <= thresholds["peak_response_relative"]
    checks["maximum_interstory_drift_ratio"] = pair["scalar_metrics"]["maximum_interstory_drift_ratio"]["relative_difference"] <= thresholds["peak_response_relative"]
    checks["peak_isolation_displacement"] = pair["scalar_metrics"]["peak_isolation_displacement"]["relative_difference"] <= thresholds["peak_response_relative"]
    checks["peak_total_base_shear"] = pair["scalar_metrics"]["peak_total_base_shear"]["relative_difference"] <= thresholds["base_shear_peak_relative"]
    checks["total_isolator_hysteresis_loop_area"] = pair["scalar_metrics"]["total_isolator_hysteresis_loop_area"]["relative_difference"] <= thresholds["hysteresis_area_relative"]
    checks["roof_displacement_nrmse"] = pair["curve_metrics"]["roof_displacement_nrmse"] <= thresholds["curve_nrmse"]
    checks["isolation_displacement_nrmse"] = pair["curve_metrics"]["isolation_displacement_nrmse"] <= thresholds["curve_nrmse"]
    checks["base_shear_nrmse"] = pair["curve_metrics"]["base_shear_nrmse"] <= thresholds["base_shear_curve_nrmse"]
    checks["convergence"] = (
        int(a.get("convergence", {}).get("failed_step_count", 1) or 0) == 0
        and int(b.get("convergence", {}).get("failed_step_count", 1) or 0) == 0
    )
    pair["checks"] = checks
    pair["pass"] = all(checks.values())
    return pair


def input_consistency(results: dict[str, dict[str, Any]]) -> dict[str, Any]:
    values = {}
    for name, data in results.items():
        values[name] = {
            "status": data.get("status"),
            "input_file": data.get("input_file"),
            "time_step": data.get("time_step"),
            "num_input_samples": data.get("num_input_samples"),
            "raw_pga": data.get("raw_pga"),
            "normalization_scale": data.get("normalization_scale"),
            "normalized_pga": data.get("normalized_pga"),
        }
    ok = True
    keys = ["time_step", "num_input_samples", "raw_pga", "normalization_scale", "normalized_pga"]
    for key in keys:
        numeric = [v[key] for v in values.values() if v.get(key) is not None]
        if len(numeric) >= 2 and max(numeric) - min(numeric) > max(1.0e-10, 1.0e-10 * max(abs(x) for x in numeric)):
            ok = False
    return {"consistent": ok, "values": values}


def main() -> None:
    paths = {
        "matlab": ROOT / "verification" / "matlab" / "matlab_results.json",
        "openseespy": ROOT / "verification" / "openseespy" / "opensees_results.json",
        "abaqus": ROOT / "verification" / "abaqus" / "abaqus_results.json",
    }
    results = {name: load_json(path) for name, path in paths.items()}
    thresholds = {
        "fundamental_period_relative": 1.0e-6,
        "peak_response_relative": 5.0e-3,
        "base_shear_peak_relative": 5.0e-3,
        "curve_nrmse": 5.0e-3,
        "base_shear_curve_nrmse": 1.0e-2,
        "hysteresis_area_relative": 1.0e-2,
        "reasoning": (
            "MATLAB and OpenSeesPy use the same physical model but independent assemblers/solvers and different Bouc-Wen implementations. "
            "A 0.5% tolerance is used for nonlinear peak responses, 1% for base-shear curve and hysteresis area because resisting-force "
            "records are sensitive to element-level output timing, and a tight 1e-6 period tolerance is used for the initial elastic eigenproblem."
        ),
    }
    pairs = {
        "matlab_vs_openseespy": compare_pair("matlab", results["matlab"], "openseespy", results["openseespy"], thresholds),
        "matlab_vs_abaqus": compare_pair("matlab", results["matlab"], "abaqus", results["abaqus"], thresholds),
        "openseespy_vs_abaqus": compare_pair("openseespy", results["openseespy"], "abaqus", results["abaqus"], thresholds),
    }
    strict_three_software_pass = (
        all(results[name].get("status") == "OK" for name in results)
        and all(pair.get("pass") for pair in pairs.values())
    )
    if strict_three_software_pass:
        failure_reason = None
    elif any(results[name].get("status") != "OK" for name in results):
        bad = ", ".join(f"{name}={results[name].get('status')}" for name in results if results[name].get("status") != "OK")
        failure_reason = f"One or more software runs did not complete with required outputs: {bad}."
    else:
        failed_pairs = ", ".join(name for name, pair in pairs.items() if not pair.get("pass"))
        failure_reason = f"All software runs completed, but strict response agreement checks failed for: {failed_pairs}."

    report = {
        "generated_in_this_run": True,
        "source_outputs_used": {name: str(path.relative_to(ROOT)) for name, path in paths.items()},
        "input_consistency": input_consistency(results),
        "thresholds": thresholds,
        "pairs": pairs,
        "software_status": {name: data.get("status") for name, data in results.items()},
        "strict_three_software_pass": strict_three_software_pass,
        "final_result": "PASS" if strict_three_software_pass else "FAIL",
        "failure_reason": failure_reason,
    }
    (ROOT / "cross_validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
