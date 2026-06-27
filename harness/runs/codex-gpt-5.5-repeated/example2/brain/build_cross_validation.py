import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "cross_validation_report.json"

RESPONSE_FILES = {
    "MATLAB": ROOT / "verification" / "matlab" / "response.json",
    "OpenSeesPy": ROOT / "verification" / "openseespy" / "response.json",
    "ABAQUS": ROOT / "verification" / "abaqus" / "response.json",
}

THRESHOLDS = {
    "max_time_vector_difference_s": 1.0e-9,
    "displacement_peak_relative_error": 0.02,
    "acceleration_peak_relative_error": 0.05,
    "displacement_normalized_l2_error": 0.03,
    "acceleration_normalized_l2_error": 0.05,
    "minimum_trend_correlation": 0.995,
}

THRESHOLD_REASONING = (
    "All models are linear and use the same mass, stiffness, Rayleigh damping, time step, "
    "normalized record, and relative-coordinate inertial forcing. MATLAB and OpenSeesPy use "
    "average-acceleration Newmark integration directly, while ABAQUS/Standard direct dynamics "
    "may differ slightly in load interpolation and history output. The adopted limits therefore "
    "require essentially identical time vectors, <=2% peak displacement error, <=5% peak "
    "absolute-acceleration error, small normalized L2 curve errors, and >=0.995 correlation."
)


def load_response(path):
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(str(path))
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    required = [
        "time",
        "top_story_displacement",
        "top_story_acceleration",
        "input_pga_after_normalization",
        "dt",
        "num_samples",
    ]
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError("{} missing keys {}".format(path, missing))
    return data


def normalized_l2(a, b):
    diff = a - b
    denom = max(float(np.linalg.norm(a)), float(np.linalg.norm(b)), 1.0e-30)
    return float(np.linalg.norm(diff) / denom)


def peak_relative_error(a, b):
    pa = float(np.max(np.abs(a)))
    pb = float(np.max(np.abs(b)))
    denom = max(abs(pa), abs(pb), 1.0e-30)
    return abs(pa - pb) / denom


def corr(a, b):
    if a.size < 2 or b.size < 2:
        return 1.0 if np.allclose(a, b) else 0.0
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa < 1.0e-30 or sb < 1.0e-30:
        return 1.0 if np.allclose(a, b) else 0.0
    return float(np.corrcoef(a, b)[0, 1])


def main():
    report = {
        "status": "FAIL",
        "thresholds": THRESHOLDS,
        "threshold_reasoning": THRESHOLD_REASONING,
        "response_files": {name: str(path.relative_to(ROOT)) for name, path in RESPONSE_FILES.items()},
        "software": {},
        "pairwise": {},
        "missing_or_invalid_outputs": [],
    }

    responses = {}
    for name, path in RESPONSE_FILES.items():
        try:
            responses[name] = load_response(path)
        except Exception as exc:
            report["missing_or_invalid_outputs"].append({"software": name, "path": str(path.relative_to(ROOT)), "error": str(exc)})

    for name, data in responses.items():
        t = np.asarray(data["time"], dtype=float)
        d = np.asarray(data["top_story_displacement"], dtype=float)
        a = np.asarray(data["top_story_acceleration"], dtype=float)
        report["software"][name] = {
            "num_samples": int(data["num_samples"]),
            "dt": float(data["dt"]),
            "input_pga_after_normalization": float(data["input_pga_after_normalization"]),
            "time_start": float(t[0]) if t.size else None,
            "time_end": float(t[-1]) if t.size else None,
            "peak_abs_top_displacement_m": float(np.max(np.abs(d))) if d.size else None,
            "peak_abs_top_acceleration_m_per_s2": float(np.max(np.abs(a))) if a.size else None,
            "json_path": str(RESPONSE_FILES[name].relative_to(ROOT)),
        }

    all_pairwise_pass = True
    for a_name, b_name in combinations(RESPONSE_FILES.keys(), 2):
        key = "{}_vs_{}".format(a_name, b_name)
        if a_name not in responses or b_name not in responses:
            report["pairwise"][key] = {"status": "FAIL", "reason": "one or both response outputs missing"}
            all_pairwise_pass = False
            continue

        ra = responses[a_name]
        rb = responses[b_name]
        ta = np.asarray(ra["time"], dtype=float)
        tb = np.asarray(rb["time"], dtype=float)
        da = np.asarray(ra["top_story_displacement"], dtype=float)
        db = np.asarray(rb["top_story_displacement"], dtype=float)
        aa = np.asarray(ra["top_story_acceleration"], dtype=float)
        ab = np.asarray(rb["top_story_acceleration"], dtype=float)

        same_len = ta.size == tb.size == da.size == db.size == aa.size == ab.size
        if same_len:
            max_time_diff = float(np.max(np.abs(ta - tb)))
            disp_peak_err = peak_relative_error(da, db)
            acc_peak_err = peak_relative_error(aa, ab)
            disp_l2 = normalized_l2(da, db)
            acc_l2 = normalized_l2(aa, ab)
            disp_corr = corr(da, db)
            acc_corr = corr(aa, ab)
        else:
            max_time_diff = None
            disp_peak_err = None
            acc_peak_err = None
            disp_l2 = None
            acc_l2 = None
            disp_corr = None
            acc_corr = None

        checks = {
            "same_length": same_len,
            "time_vector": same_len and max_time_diff <= THRESHOLDS["max_time_vector_difference_s"],
            "displacement_peak": same_len and disp_peak_err <= THRESHOLDS["displacement_peak_relative_error"],
            "acceleration_peak": same_len and acc_peak_err <= THRESHOLDS["acceleration_peak_relative_error"],
            "displacement_curve": same_len and disp_l2 <= THRESHOLDS["displacement_normalized_l2_error"],
            "acceleration_curve": same_len and acc_l2 <= THRESHOLDS["acceleration_normalized_l2_error"],
            "displacement_trend": same_len and disp_corr >= THRESHOLDS["minimum_trend_correlation"],
            "acceleration_trend": same_len and acc_corr >= THRESHOLDS["minimum_trend_correlation"],
        }
        pair_pass = all(checks.values())
        if not pair_pass:
            all_pairwise_pass = False

        report["pairwise"][key] = {
            "status": "PASS" if pair_pass else "FAIL",
            "checks": checks,
            "max_time_vector_difference_s": max_time_diff,
            "displacement_peak_relative_error": disp_peak_err,
            "acceleration_peak_relative_error": acc_peak_err,
            "displacement_normalized_l2_error": disp_l2,
            "acceleration_normalized_l2_error": acc_l2,
            "displacement_correlation": disp_corr,
            "acceleration_correlation": acc_corr,
        }

    required_outputs_present = len(responses) == len(RESPONSE_FILES) and not report["missing_or_invalid_outputs"]
    report["strict_three_software_pass"] = bool(required_outputs_present and all_pairwise_pass)
    report["status"] = "PASS" if report["strict_three_software_pass"] else "FAIL"

    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "strict_three_software_pass": report["strict_three_software_pass"],
        "missing_or_invalid_outputs": report["missing_or_invalid_outputs"],
    }, indent=2))


if __name__ == "__main__":
    main()
