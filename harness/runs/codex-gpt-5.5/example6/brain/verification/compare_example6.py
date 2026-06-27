#!/usr/bin/env python3
"""Cross-validate Example 6 response outputs.

The script reads available response_output.json files from MATLAB, OpenSeesPy,
and ABAQUS verification directories and writes cross_validation_report.json at
the run root.  PASS is allowed only when all three software outputs are present,
completed, and within the requested thresholds.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = SCRIPT_DIR.parent
REPORT_PATH = RUN_ROOT / "cross_validation_report.json"

SOFTWARE = {
    "matlab": RUN_ROOT / "verification" / "matlab" / "response_output.json",
    "openseespy": RUN_ROOT / "verification" / "openseespy" / "response_output.json",
    "abaqus": RUN_ROOT / "verification" / "abaqus" / "response_output.json",
}

THRESHOLDS = {
    "fundamental_period_relative_error": 0.03,
    "peak_roof_displacement_relative_error": 0.10,
    "peak_interstory_drift_relative_error": 0.10,
    "peak_isolation_displacement_relative_error": 0.10,
    "peak_base_shear_relative_error": 0.15,
    "roof_displacement_nrmse": 0.05,
    "isolation_displacement_nrmse": 0.05,
    "base_shear_nrmse": 0.10,
    "hysteresis_loop_area_relative_error": 0.15,
    "failed_step_rate": 0.05,
}

PEAK_KEYS = {
    "peak_roof_displacement_relative_error": "peak_abs_roof_displacement_m",
    "peak_interstory_drift_relative_error": "max_abs_story_drift_ratio",
    "peak_isolation_displacement_relative_error": "peak_abs_isolation_displacement_m",
    "peak_base_shear_relative_error": "peak_abs_total_base_shear_N",
    "hysteresis_loop_area_relative_error": "total_isolator_loop_area_abs_Nm",
}


def load_outputs() -> tuple[dict[str, Any], dict[str, Any]]:
    outputs: dict[str, Any] = {}
    status: dict[str, Any] = {}
    for name, path in SOFTWARE.items():
        if not path.exists():
            status[name] = {"status": "missing", "path": rel(path)}
            continue
        try:
            with path.open("r", encoding="utf-8") as f:
                payload = json.load(f)
            outputs[name] = payload
            status[name] = {
                "status": payload.get("status", "unknown"),
                "path": rel(path),
                "message": payload.get("message"),
            }
        except Exception as exc:
            status[name] = {
                "status": "unreadable",
                "path": rel(path),
                "message": f"{type(exc).__name__}: {exc}",
            }
    return outputs, status


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(RUN_ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def as_array(payload: dict[str, Any], key: str) -> np.ndarray:
    return np.asarray(payload.get(key, []), dtype=float).reshape(-1)


def story_matrix(payload: dict[str, Any]) -> np.ndarray:
    obj = payload.get("story_drift_ratios", {})
    if isinstance(obj, dict):
        return np.column_stack(
            [
                np.asarray(obj.get("story_1", []), dtype=float).reshape(-1),
                np.asarray(obj.get("story_2", []), dtype=float).reshape(-1),
                np.asarray(obj.get("story_3", []), dtype=float).reshape(-1),
            ]
        )
    return np.asarray(obj, dtype=float)


def safe_float(value: Any) -> float | None:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    return x


def rel_error(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    denom = max(abs(a), abs(b), 1.0e-30)
    return abs(a - b) / denom


def interpolate_to_common(a_payload: dict[str, Any], b_payload: dict[str, Any], key: str) -> tuple[np.ndarray, np.ndarray]:
    ta = as_array(a_payload, "time_s")
    tb = as_array(b_payload, "time_s")
    ya = as_array(a_payload, key)
    yb = as_array(b_payload, key)
    if ta.size == 0 or tb.size == 0 or ya.size == 0 or yb.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    n = min(ta.size, ya.size)
    ta = ta[:n]
    ya = ya[:n]
    n = min(tb.size, yb.size)
    tb = tb[:n]
    yb = yb[:n]
    t0 = max(float(ta[0]), float(tb[0]))
    t1 = min(float(ta[-1]), float(tb[-1]))
    if t1 < t0:
        return np.array([], dtype=float), np.array([], dtype=float)
    dt_candidates = []
    if ta.size > 1:
        dt_candidates.append(float(np.median(np.diff(ta))))
    if tb.size > 1:
        dt_candidates.append(float(np.median(np.diff(tb))))
    dt = min([x for x in dt_candidates if x > 0.0], default=0.01)
    common = np.arange(t0, t1 + 0.5 * dt, dt)
    return np.interp(common, ta, ya), np.interp(common, tb, yb)


def nrmse(a_payload: dict[str, Any], b_payload: dict[str, Any], key: str) -> float | None:
    ya, yb = interpolate_to_common(a_payload, b_payload, key)
    if ya.size == 0 or yb.size == 0:
        return None
    rmse = float(np.sqrt(np.mean((ya - yb) ** 2)))
    norm = max(float(np.max(np.abs(ya))), float(np.max(np.abs(yb))), 1.0e-30)
    return rmse / norm


def pairwise_metric(
    completed: dict[str, dict[str, Any]],
    metric_name: str,
    getter,
    threshold: float,
) -> dict[str, Any]:
    pairs: dict[str, Any] = {}
    values: list[float] = []
    for a, b in combinations(completed.keys(), 2):
        value = getter(completed[a], completed[b])
        pass_pair = value is not None and value <= threshold
        pairs[f"{a}_vs_{b}"] = {"value": value, "threshold": threshold, "pass": pass_pair}
        if value is not None:
            values.append(value)
    metric_pass = bool(pairs) and all(item["pass"] for item in pairs.values())
    return {
        "metric": metric_name,
        "threshold": threshold,
        "max_pairwise_value": max(values) if values else None,
        "pass": metric_pass,
        "pairs": pairs,
    }


def validate_inputs(completed: dict[str, dict[str, Any]]) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    for name, payload in completed.items():
        inp = payload.get("input", {})
        pga = safe_float(inp.get("normalized_peak_abs_m_per_s2"))
        target = safe_float(inp.get("target_pga_m_per_s2"))
        dt = safe_float(inp.get("dt_s"))
        sample_count = inp.get("sample_count")
        pga_ok = pga is not None and target is not None and abs(pga - target) <= max(1.0e-8, 1.0e-8 * abs(target))
        dt_ok = dt is not None and abs(dt - 0.01) <= 1.0e-12
        checks[name] = {
            "normalized_peak_abs_m_per_s2": pga,
            "target_pga_m_per_s2": target,
            "pga_normalization_pass": pga_ok,
            "dt_s": dt,
            "dt_pass": dt_ok,
            "sample_count": sample_count,
            "pass": pga_ok and dt_ok and sample_count is not None,
        }
    return {"pass": bool(checks) and all(x["pass"] for x in checks.values()), "software": checks}


def convergence_checks(completed: dict[str, dict[str, Any]]) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    for name, payload in completed.items():
        conv = payload.get("convergence", {})
        rate = safe_float(conv.get("failed_step_rate"))
        checks[name] = {
            "failed_step_rate": rate,
            "threshold": THRESHOLDS["failed_step_rate"],
            "pass": rate is not None and rate <= THRESHOLDS["failed_step_rate"],
            "failed_step_count": conv.get("failed_step_count"),
            "completed_steps": conv.get("completed_steps"),
        }
    return {"pass": bool(checks) and all(x["pass"] for x in checks.values()), "software": checks}


def build_metrics(completed: dict[str, dict[str, Any]]) -> dict[str, Any]:
    metrics: dict[str, Any] = {}
    metrics["fundamental_period_relative_error"] = pairwise_metric(
        completed,
        "fundamental_period_relative_error",
        lambda a, b: rel_error(safe_float(a.get("fundamental_period_s")), safe_float(b.get("fundamental_period_s"))),
        THRESHOLDS["fundamental_period_relative_error"],
    )

    for metric_name, peak_key in PEAK_KEYS.items():
        metrics[metric_name] = pairwise_metric(
            completed,
            metric_name,
            lambda a, b, key=peak_key: rel_error(
                safe_float(a.get("peaks", {}).get(key)),
                safe_float(b.get("peaks", {}).get(key)),
            ),
            THRESHOLDS[metric_name],
        )

    metrics["roof_displacement_nrmse"] = pairwise_metric(
        completed,
        "roof_displacement_nrmse",
        lambda a, b: nrmse(a, b, "roof_displacement_m"),
        THRESHOLDS["roof_displacement_nrmse"],
    )
    metrics["isolation_displacement_nrmse"] = pairwise_metric(
        completed,
        "isolation_displacement_nrmse",
        lambda a, b: nrmse(a, b, "isolation_displacement_m"),
        THRESHOLDS["isolation_displacement_nrmse"],
    )
    metrics["base_shear_nrmse"] = pairwise_metric(
        completed,
        "base_shear_nrmse",
        lambda a, b: nrmse(a, b, "total_base_shear_N"),
        THRESHOLDS["base_shear_nrmse"],
    )
    return metrics


def summarize_completed(outputs: dict[str, Any]) -> dict[str, dict[str, Any]]:
    completed: dict[str, dict[str, Any]] = {}
    for name, payload in outputs.items():
        if payload.get("status") == "completed":
            completed[name] = payload
    return completed


def main() -> int:
    outputs, status = load_outputs()
    completed = summarize_completed(outputs)
    required_complete = set(SOFTWARE) == set(completed)
    input_consistency = validate_inputs(completed) if completed else {"pass": False, "software": {}}
    convergence = convergence_checks(completed) if completed else {"pass": False, "software": {}}
    metrics = build_metrics(completed) if len(completed) >= 2 else {}
    metrics_pass = bool(metrics) and all(item["pass"] for item in metrics.values())

    blockers = []
    for name in SOFTWARE:
        if name not in completed:
            blockers.append(f"{name} response_output.json is not completed")
    if not input_consistency["pass"]:
        blockers.append("input consistency check failed or is incomplete")
    if not convergence["pass"]:
        blockers.append("convergence check failed or is incomplete")
    if not metrics_pass:
        blockers.append("cross-software metric thresholds failed or are incomplete")

    overall_pass = required_complete and input_consistency["pass"] and convergence["pass"] and metrics_pass
    verdict = "PASS" if overall_pass else "FAIL_INCOMPLETE_OR_OUT_OF_TOLERANCE"
    if overall_pass:
        conclusion = "MATLAB, OpenSeesPy, and ABAQUS completed Example 6 and satisfied all recorded cross-software agreement thresholds."
    else:
        missing_or_incomplete = [
            f"{name}: {status.get(name, {}).get('status', 'missing')}"
            for name in SOFTWARE
            if name not in completed
        ]
        detail = "; ".join(missing_or_incomplete) if missing_or_incomplete else "one or more threshold checks failed"
        conclusion = (
            "MATLAB and OpenSeesPy were checked where available, but full auditable MATLAB/OpenSeesPy/ABAQUS "
            f"agreement was not established because {detail}. The ABAQUS path must not be promoted to a pass "
            "unless it implements the same Bouc-Wen state evolution rather than a bilinear substitute."
        )

    report = {
        "schema_version": "example6-cross-validation-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": str(RUN_ROOT),
        "software_status": status,
        "required_software": list(SOFTWARE.keys()),
        "completed_software": list(completed.keys()),
        "thresholds": THRESHOLDS,
        "input_consistency": input_consistency,
        "metrics": metrics,
        "convergence": convergence,
        "overall_pass": overall_pass,
        "overall_verdict": verdict,
        "blockers": blockers,
        "conclusion": conclusion,
        "policy": "PASS is allowed only when MATLAB, OpenSeesPy, and ABAQUS all complete and every threshold check passes.",
    }
    with REPORT_PATH.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
