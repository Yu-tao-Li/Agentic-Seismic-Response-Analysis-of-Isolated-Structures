import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RESPONSES = {
    "matlab": ROOT / "verification" / "matlab" / "response_matlab.json",
    "openseespy": ROOT / "verification" / "openseespy" / "response_openseespy.json",
    "abaqus": ROOT / "verification" / "abaqus" / "response_abaqus.json",
}


def load_json(path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def sanitize_for_json(value):
    if isinstance(value, dict):
        return {k: sanitize_for_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [sanitize_for_json(v) for v in value]
    if isinstance(value, tuple):
        return [sanitize_for_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return sanitize_for_json(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def finite_array(values):
    return np.asarray(values, dtype=float)


def clean_nan(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return v


def rel_error(a, b, floor=1.0e-12):
    denom = max(abs(a), abs(b), floor)
    return abs(a - b) / denom


def equivalent_model_value(a, b, rel_tol=1.0e-10, abs_tol=1.0e-8):
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        if not isinstance(a, (list, tuple)) or not isinstance(b, (list, tuple)) or len(a) != len(b):
            return False
        return all(equivalent_model_value(x, y, rel_tol, abs_tol) for x, y in zip(a, b))
    try:
        af = float(a)
        bf = float(b)
    except (TypeError, ValueError):
        return a == b
    if not math.isfinite(af) or not math.isfinite(bf):
        return (not math.isfinite(af)) and (not math.isfinite(bf))
    return abs(af - bf) <= max(abs_tol, rel_tol * max(abs(af), abs(bf), 1.0))


def curve_metrics(a_time, a_values, b_time, b_values):
    a_time = finite_array(a_time)
    b_time = finite_array(b_time)
    a_values = finite_array(a_values)
    b_values = finite_array(b_values)
    start = max(float(a_time[0]), float(b_time[0]))
    end = min(float(a_time[-1]), float(b_time[-1]))
    dt = max(np.median(np.diff(a_time)), np.median(np.diff(b_time)))
    common = np.arange(start, end + 0.25 * dt, dt)
    av = np.interp(common, a_time, a_values)
    bv = np.interp(common, b_time, b_values)
    diff = av - bv
    scale = max(float(np.max(np.abs(av))), float(np.max(np.abs(bv))), 1.0e-12)
    return {
        "sample_count": int(common.size),
        "rms_abs": float(np.sqrt(np.mean(diff**2))),
        "rms_norm_by_pair_peak": float(np.sqrt(np.mean(diff**2)) / scale),
        "max_abs": float(np.max(np.abs(diff))),
        "max_norm_by_pair_peak": float(np.max(np.abs(diff)) / scale),
    }


def peak_pair_metrics(ra, rb):
    keys = [
        "peak_top_displacement_m",
        "peak_top_acceleration_mps2",
        "peak_isolation_displacement_m",
        "peak_base_shear_N",
    ]
    out = {}
    for key in keys:
        a = clean_nan(ra["summary"].get(key))
        b = clean_nan(rb["summary"].get(key))
        out[key] = {
            "a": a,
            "b": b,
            "abs_error": None if a is None or b is None else abs(a - b),
            "relative_error": None if a is None or b is None else rel_error(a, b),
        }
    return out


def first_yield(response):
    values = response.get("summary", {}).get("first_yield_time_s")
    if values is None:
        values = response.get("yielding", {}).get("first_yield_time", [])
    return [clean_nan(v) for v in values]


def yield_metrics(ra, rb):
    a = first_yield(ra)
    b = first_yield(rb)
    max_len = max(len(a), len(b))
    a += [None] * (max_len - len(a))
    b += [None] * (max_len - len(b))
    per_story = []
    all_match = True
    max_time_error = 0.0
    for i, (ta, tb) in enumerate(zip(a, b), start=2):
        yielded_match = (ta is None) == (tb is None)
        if ta is not None and tb is not None:
            terr = abs(ta - tb)
            max_time_error = max(max_time_error, terr)
        else:
            terr = None
        if not yielded_match:
            all_match = False
        per_story.append({
            "story": i,
            "a_first_yield_time_s": ta,
            "b_first_yield_time_s": tb,
            "yielded_match": yielded_match,
            "time_error_s": terr,
        })
    return {
        "per_story": per_story,
        "all_yielded_flags_match": all_match,
        "max_first_yield_time_error_s": max_time_error,
    }


def input_consistency(responses):
    fields = ["raw_count", "dt", "scale_to_mps2", "normalized_pga_mps2", "target_pga_g"]
    model_fields = ["masses_kg", "story_stiffness_N_per_m", "yield_forces_N", "post_yield_ratios", "isolation_dashpot_Ns_per_m"]
    details = {}
    ok = True
    for field in fields:
        values = {name: responses[name]["input"].get(field) for name in responses}
        details[field] = values
        finite = [float(v) for v in values.values() if v is not None]
        if len(finite) != len(values) or max(finite) - min(finite) > max(1.0e-10, 1.0e-10 * max(abs(x) for x in finite)):
            ok = False
    for field in model_fields:
        values = {name: responses[name]["model"].get(field) for name in responses}
        details[field] = values
        ref = next(iter(values.values()))
        if any(not equivalent_model_value(ref, v) for v in values.values()):
            ok = False
    rayleigh = {
        name: {
            "alpha_mass": responses[name]["model"]["rayleigh"].get("alpha_mass"),
            "beta_initial_stiffness": responses[name]["model"]["rayleigh"].get("beta_initial_stiffness"),
        }
        for name in responses
    }
    details["rayleigh_alpha_beta"] = rayleigh
    alphas = [float(v["alpha_mass"]) for v in rayleigh.values()]
    betas = [float(v["beta_initial_stiffness"]) for v in rayleigh.values()]
    if max(alphas) - min(alphas) > 1.0e-10 or max(betas) - min(betas) > 1.0e-12:
        ok = False
    return {"consistent": ok, "details": details}


def file_manifest():
    items = {}
    for folder in ["matlab", "openseespy", "abaqus"]:
        base = ROOT / "verification" / folder
        files = {}
        for path in sorted(base.iterdir()):
            if path.is_file():
                files[path.name] = {"bytes": path.stat().st_size, "relative_path": str(path.relative_to(ROOT)).replace("\\", "/")}
        items[folder] = files
    return items


def main():
    responses = {name: load_json(path) for name, path in RESPONSES.items()}
    thresholds = {
        "input_consistency_required": True,
        "peak_relative_error_limit": 0.01,
        "curve_rms_norm_limit": 0.02,
        "curve_max_norm_limit": 0.06,
        "isolation_peak_relative_error_limit": 0.01,
        "first_yield_time_error_limit_s": 0.02,
        "convergence_failed_step_limit": 0,
        "reasoning": (
            "MATLAB and OpenSeesPy use average-acceleration Newmark and should match nearly exactly. "
            "ABAQUS/Standard is configured with HHT alpha = 0.0, the Newmark average-acceleration limit. "
            "The retained 1% peak and 2% normalized RMS curve limits are tight enough to catch modeling errors "
            "while allowing small solver implementation and ODB interpolation differences."
        ),
    }

    pairwise = {}
    for a, b in combinations(responses.keys(), 2):
        ra = responses[a]
        rb = responses[b]
        curves = {}
        for field in ["top_story_displacement", "top_story_acceleration", "isolation_displacement", "base_shear"]:
            curves[field] = curve_metrics(ra["time"], ra[field], rb["time"], rb[field])
        peak = peak_pair_metrics(ra, rb)
        yld = yield_metrics(ra, rb)
        conv = {
            a: int(ra["convergence"].get("failed_step_count", -1)),
            b: int(rb["convergence"].get("failed_step_count", -1)),
        }
        pass_checks = {
            "peak_errors_within_limit": all(
                v["relative_error"] is not None and v["relative_error"] <= thresholds["peak_relative_error_limit"]
                for v in peak.values()
            ),
            "curve_rms_within_limit": all(v["rms_norm_by_pair_peak"] <= thresholds["curve_rms_norm_limit"] for v in curves.values()),
            "curve_max_within_limit": all(v["max_norm_by_pair_peak"] <= thresholds["curve_max_norm_limit"] for v in curves.values()),
            "isolation_peak_within_limit": (
                peak["peak_isolation_displacement_m"]["relative_error"] is not None
                and peak["peak_isolation_displacement_m"]["relative_error"] <= thresholds["isolation_peak_relative_error_limit"]
            ),
            "yield_sequence_within_limit": (
                yld["all_yielded_flags_match"]
                and yld["max_first_yield_time_error_s"] <= thresholds["first_yield_time_error_limit_s"]
            ),
            "convergence_within_limit": all(v <= thresholds["convergence_failed_step_limit"] for v in conv.values()),
        }
        pairwise[f"{a}_vs_{b}"] = {
            "peak_response_error": peak,
            "curve_error": curves,
            "yielding_sequence": yld,
            "convergence": conv,
            "pass_checks": pass_checks,
            "pair_pass": all(pass_checks.values()),
        }

    software_status = {
        name: {
            "response_file": str(RESPONSES[name].relative_to(ROOT)).replace("\\", "/"),
            "response_file_bytes": RESPONSES[name].stat().st_size,
            "failed_step_count": int(resp.get("convergence", {}).get("failed_step_count", -1)),
            "time_points": len(resp.get("time", [])),
            "topStoDisIso1_bytes": (RESPONSES[name].parent / "topStoDisIso1.txt").stat().st_size,
            "topStoAccIso1_bytes": (RESPONSES[name].parent / "topStoAccIso1.txt").stat().st_size,
        }
        for name, resp in responses.items()
    }

    input_ok = input_consistency(responses)
    final_pass = (
        input_ok["consistent"]
        and all(v["pair_pass"] for v in pairwise.values())
        and all(s["response_file_bytes"] > 0 and s["topStoDisIso1_bytes"] > 0 and s["topStoAccIso1_bytes"] > 0 for s in software_status.values())
    )

    report = {
        "generated_at": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "basis": "Only generated outputs under verification/matlab, verification/openseespy, and verification/abaqus were compared.",
        "thresholds": thresholds,
        "software_status": software_status,
        "input_consistency": input_ok,
        "summaries": {name: responses[name]["summary"] for name in responses},
        "pairwise_metrics": pairwise,
        "file_manifest": file_manifest(),
        "final_result": {
            "strict_three_software_pass": final_pass,
            "reason": "All required generated outputs are present and pairwise metrics meet adopted thresholds." if final_pass else "One or more required outputs or thresholds failed.",
        },
    }

    report = sanitize_for_json(report)
    (ROOT / "cross_validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")

    lines = []
    lines.append("# Cross-Validation Report")
    lines.append("")
    lines.append("## Scope")
    lines.append("Generated and executed three independent analysis paths for the four-story isolated nonlinear shear-building model using the supplied `Northridge_01_NO_968.txt` record, `dt = 0.01 s`, and PGA normalized to `0.40 g`.")
    lines.append("")
    lines.append("## Model and Numerical Assumptions")
    lines.append("- Masses, stiffnesses, yield forces, and post-yield ratios match the fixed prompt values in all generated inputs.")
    lines.append("- Damping uses `zeta = 0.05` Rayleigh coefficients from initial modes 1 and 4, implemented as `alpha*M + beta*K_initial`; the isolation dashpot `c = 1000e3 N*s/m` is added explicitly.")
    lines.append("- MATLAB uses average-acceleration Newmark with Newton iteration and return mapping for ideal elastic-plastic stories.")
    lines.append("- OpenSeesPy uses zeroLength springs with Elastic/ElasticPP materials and explicit viscous elements for the same damping matrix.")
    lines.append("- ABAQUS uses lumped mass elements, T2D2 truss spring elements with elastic-perfectly-plastic material behavior, DASHPOT2 elements for damping, and `*DYNAMIC, ALPHA=0.0` to match the Newmark average-acceleration limit. A negligible density `1.0e-12 kg/m^3` was added to truss materials because ABAQUS/Standard requires at least one active material density in general dynamic steps; inertia remains governed by the specified lumped masses.")
    lines.append("")
    lines.append("## Execution Notes and Corrections")
    lines.append("- MATLAB completed with zero failed steps.")
    lines.append("- The default Python 3.13 environment contained an OpenSeesPy package linked to `python312.dll`, causing an import failure. The analysis was rerun successfully under `<python_env>/python.exe` with OpenSeesPy 3.8.0.")
    lines.append("- ABAQUS first failed preprocessing due to zero active material density. After adding negligible truss density and setting `ALPHA=0.0` for Newmark-equivalent integration, ABAQUS completed successfully with 1899 increments and zero cutbacks. The wrapper was updated to verify success from the `.sta` file before postprocessing.")
    lines.append("")
    lines.append("## Consistency Criteria")
    lines.append(f"- Peak response relative error limit: {thresholds['peak_relative_error_limit']:.1%}.")
    lines.append(f"- Normalized RMS curve error limit: {thresholds['curve_rms_norm_limit']:.1%}; normalized max curve error limit: {thresholds['curve_max_norm_limit']:.1%}.")
    lines.append(f"- Isolation peak relative error limit: {thresholds['isolation_peak_relative_error_limit']:.1%}.")
    lines.append(f"- First-yield time tolerance: {thresholds['first_yield_time_error_limit_s']} s; failed step count must be zero.")
    lines.append(f"- Reasoning: {thresholds['reasoning']}")
    lines.append("")
    lines.append("## Generated Outputs Used")
    for name, status in software_status.items():
        lines.append(f"- {name}: `{status['response_file']}`, `verification/{name}/topStoDisIso1.txt`, `verification/{name}/topStoAccIso1.txt`; time points = {status['time_points']}, failed steps = {status['failed_step_count']}.")
    lines.append("")
    lines.append("## Peak Response Summary")
    lines.append("| Software | Peak top disp (m) | Peak top acc (m/s^2) | Peak iso disp (m) | Peak base shear (N) | First yield times s2/s3/s4 |")
    lines.append("|---|---:|---:|---:|---:|---|")
    for name, resp in responses.items():
        s = resp["summary"]
        fy = first_yield(resp)
        fy_text = "/".join("none" if v is None else f"{v:.3f}" for v in fy)
        lines.append(
            f"| {name} | {s['peak_top_displacement_m']:.8f} | {s['peak_top_acceleration_mps2']:.8f} | "
            f"{s['peak_isolation_displacement_m']:.8f} | {s['peak_base_shear_N']:.3f} | {fy_text} |"
        )
    lines.append("")
    lines.append("## Pairwise Metrics")
    lines.append("| Pair | Max peak rel error | Max RMS curve norm | Max curve norm | Yield max time err (s) | Pair pass |")
    lines.append("|---|---:|---:|---:|---:|---|")
    for pair, metrics in pairwise.items():
        max_peak = max(v["relative_error"] for v in metrics["peak_response_error"].values())
        max_rms = max(v["rms_norm_by_pair_peak"] for v in metrics["curve_error"].values())
        max_curve = max(v["max_norm_by_pair_peak"] for v in metrics["curve_error"].values())
        yerr = metrics["yielding_sequence"]["max_first_yield_time_error_s"]
        lines.append(f"| {pair} | {max_peak:.6e} | {max_rms:.6e} | {max_curve:.6e} | {yerr:.6e} | {metrics['pair_pass']} |")
    lines.append("")
    lines.append("## Final Verdict")
    lines.append(f"Strict three-software pass: **{final_pass}**.")
    lines.append(report["final_result"]["reason"])
    lines.append("")
    lines.append("Detailed machine-readable metrics are in `cross_validation_report.json`.")
    (ROOT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
