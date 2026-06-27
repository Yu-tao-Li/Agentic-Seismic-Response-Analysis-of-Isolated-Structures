"""Cross-validate generated MATLAB, OpenSeesPy, and ABAQUS response outputs."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
VERIFICATION = ROOT / "verification"
SOFTWARE_DIRS = {
    "matlab": VERIFICATION / "matlab",
    "openseespy": VERIFICATION / "openseespy",
    "abaqus": VERIFICATION / "abaqus",
}


def load_response(name: str) -> dict:
    path = SOFTWARE_DIRS[name] / "response.json"
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty response output: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def interp_to(time_ref: np.ndarray, response: dict, key: str) -> np.ndarray:
    time = np.asarray(response["time"], dtype=float)
    values = np.asarray(response[key], dtype=float)
    return np.interp(time_ref, time, values)


def peak_error(a: np.ndarray, b: np.ndarray) -> dict:
    pa = float(np.max(np.abs(a)))
    pb = float(np.max(np.abs(b)))
    denom = max(abs(pa), abs(pb), 1.0e-30)
    return {"a_peak": pa, "b_peak": pb, "relative_error": abs(pa - pb) / denom, "absolute_error": abs(pa - pb)}


def curve_error(a: np.ndarray, b: np.ndarray) -> dict:
    diff = a - b
    denom = max(float(np.sqrt(np.mean(a * a))), float(np.sqrt(np.mean(b * b))), 1.0e-30)
    return {
        "rms_absolute": float(np.sqrt(np.mean(diff * diff))),
        "rms_relative": float(np.sqrt(np.mean(diff * diff)) / denom),
        "max_absolute": float(np.max(np.abs(diff))),
    }


def first_yield_times(response: dict) -> list[float | None]:
    time = np.asarray(response["time"], dtype=float)
    yielding = np.asarray(response.get("yielding", []), dtype=bool)
    if yielding.ndim != 2 or yielding.shape[1] != 4:
        return [None, None, None, None]
    out: list[float | None] = []
    for i in range(4):
        idx = np.flatnonzero(yielding[:, i])
        out.append(float(time[idx[0]]) if idx.size else None)
    return out


def compare_pair(a_name: str, b_name: str, responses: dict[str, dict], thresholds: dict) -> dict:
    a = responses[a_name]
    b = responses[b_name]
    time_a = np.asarray(a["time"], dtype=float)
    time_b = np.asarray(b["time"], dtype=float)
    t0 = max(float(time_a[0]), float(time_b[0]))
    t1 = min(float(time_a[-1]), float(time_b[-1]))
    dt = max(float(np.median(np.diff(time_a))), float(np.median(np.diff(time_b))))
    time_ref = np.arange(t0, t1 + 0.5 * dt, dt)

    metrics = {}
    for key in ["top_story_displacement", "top_story_acceleration", "isolation_displacement", "base_shear"]:
        av = interp_to(time_ref, a, key)
        bv = interp_to(time_ref, b, key)
        metrics[key] = {"peak": peak_error(av, bv), "curve": curve_error(av, bv)}

    area_a = float(a["isolation_hysteresis"]["loop_area"])
    area_b = float(b["isolation_hysteresis"]["loop_area"])
    metrics["hysteresis_loop_area"] = {
        "a": area_a,
        "b": area_b,
        "relative_error": abs(area_a - area_b) / max(abs(area_a), abs(area_b), 1.0e-30),
        "absolute_error": abs(area_a - area_b),
    }
    y_a = first_yield_times(a)
    y_b = first_yield_times(b)
    yield_diffs = []
    for ya, yb in zip(y_a, y_b):
        if ya is None or yb is None:
            yield_diffs.append(0.0 if ya == yb else float("inf"))
        else:
            yield_diffs.append(abs(ya - yb))
    metrics["yielding_sequence"] = {
        "a_first_yield_times": y_a,
        "b_first_yield_times": y_b,
        "absolute_time_differences": yield_diffs,
    }

    checks = {
        "peak_top_displacement": metrics["top_story_displacement"]["peak"]["relative_error"] <= thresholds["peak_relative"],
        "peak_top_acceleration": metrics["top_story_acceleration"]["peak"]["relative_error"] <= thresholds["acceleration_peak_relative"],
        "curve_top_displacement": metrics["top_story_displacement"]["curve"]["rms_relative"] <= thresholds["curve_rms_relative"],
        "curve_top_acceleration": metrics["top_story_acceleration"]["curve"]["rms_relative"] <= thresholds["acceleration_curve_rms_relative"],
        "isolation_displacement": metrics["isolation_displacement"]["peak"]["relative_error"] <= thresholds["peak_relative"],
        "base_shear": metrics["base_shear"]["peak"]["relative_error"] <= thresholds["base_shear_peak_relative"],
        "hysteresis_loop_area": metrics["hysteresis_loop_area"]["relative_error"] <= thresholds["hysteresis_area_relative"],
        "yielding_sequence": all(d is not None and np.isfinite(d) and d <= thresholds["yield_time_abs_s"] for d in yield_diffs),
    }
    return {"time_window": [t0, t1], "metrics": metrics, "checks": checks, "pass": bool(all(checks.values()))}


def input_consistency(responses: dict[str, dict]) -> dict:
    keys = ["dt", "sample_count", "target_pga_mps2", "normalized_peak_mps2"]
    result = {}
    for key in keys:
        vals = {name: resp["input"].get(key) for name, resp in responses.items()}
        numeric = [float(v) for v in vals.values() if v is not None]
        consistent = len(numeric) == 3 and max(numeric) - min(numeric) <= max(1e-12, 1e-10 * max(abs(x) for x in numeric))
        result[key] = {"values": vals, "consistent": consistent}
    model_keys = ["masses_kg", "stiffness_N_per_m", "yield_N", "post_yield_ratio"]
    for key in model_keys:
        vals = {name: resp["model"].get(key) for name, resp in responses.items()}
        arrays = [np.asarray(v, dtype=float) for v in vals.values()]
        consistent = all(np.allclose(arrays[0], arr, rtol=0.0, atol=1e-9) for arr in arrays[1:])
        result[key] = {"values": vals, "consistent": bool(consistent)}
    return result


def main() -> None:
    responses = {name: load_response(name) for name in SOFTWARE_DIRS}
    thresholds = {
        "reasoning": (
            "MATLAB and OpenSeesPy use the same Newmark average-acceleration target model, while ABAQUS "
            "uses its direct-integration connector implementation and ODB frame interpolation. For nonlinear "
            "hysteretic response, 5% displacement/base-shear peak error, 10% acceleration peak error, "
            "5% displacement curve RMS error, 15% acceleration curve RMS error, 10% hysteresis area error, "
            "and 0.02 s first-yield timing tolerance were adopted before declaring pass/fail."
        ),
        "peak_relative": 0.05,
        "acceleration_peak_relative": 0.10,
        "curve_rms_relative": 0.05,
        "acceleration_curve_rms_relative": 0.15,
        "base_shear_peak_relative": 0.05,
        "hysteresis_area_relative": 0.10,
        "yield_time_abs_s": 0.02,
    }
    consistency = input_consistency(responses)
    pairwise = {}
    for a, b in combinations(SOFTWARE_DIRS.keys(), 2):
        pairwise[f"{a}_vs_{b}"] = compare_pair(a, b, responses, thresholds)

    output_files = {
        name: sorted(str(p.relative_to(ROOT)) for p in SOFTWARE_DIRS[name].glob("*") if p.is_file())
        for name in SOFTWARE_DIRS
    }
    generated_ok = {
        name: (SOFTWARE_DIRS[name] / "response.json").exists()
        and (SOFTWARE_DIRS[name] / "response.json").stat().st_size > 0
        and (SOFTWARE_DIRS[name] / "topStoDisIso2.txt").exists()
        and (SOFTWARE_DIRS[name] / "topStoAccIso2.txt").exists()
        for name in SOFTWARE_DIRS
    }
    convergence_ok = {
        name: int(responses[name].get("convergence", {}).get("failed_step_count", 0)) == 0
        for name in SOFTWARE_DIRS
    }
    final_pass = (
        all(generated_ok.values())
        and all(convergence_ok.values())
        and all(item["consistent"] for item in consistency.values())
        and all(pair["pass"] for pair in pairwise.values())
    )
    report = {
        "status": "PASS" if final_pass else "FAIL",
        "strict_three_software_pass": bool(final_pass),
        "generated_output_files": output_files,
        "input_consistency": consistency,
        "thresholds": thresholds,
        "pairwise": pairwise,
        "convergence_ok": convergence_ok,
        "generated_outputs_ok": generated_ok,
    }
    (ROOT / "cross_validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=True), encoding="utf-8")
    print(f"Cross-validation status: {report['status']}")


if __name__ == "__main__":
    main()
