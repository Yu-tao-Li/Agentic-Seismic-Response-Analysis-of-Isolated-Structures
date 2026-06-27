"""Cross-validate generated MATLAB, OpenSeesPy, and ABAQUS modal outputs."""

from __future__ import annotations

import json
import math
import os
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOFTWARE_OUTPUTS = {
    "matlab": ROOT / "verification" / "matlab" / "modal_results.json",
    "openseespy": ROOT / "verification" / "openseespy" / "modal_results.json",
    "abaqus": ROOT / "verification" / "abaqus" / "modal_results.json",
}

THRESHOLDS = {
    "circular_frequency_relative_tolerance": 5.0e-5,
    "period_relative_tolerance": 5.0e-5,
    "mode_shape_component_abs_tolerance": 1.0e-4,
    "mode_shape_mac_minimum": 0.999999,
}

THRESHOLD_REASONING = (
    "All three implementations represent the same three-DOF linear undamped "
    "eigenproblem in SI units. The exact model is small and well conditioned, "
    "so MATLAB and OpenSeesPy should agree to near roundoff. In this ABAQUS "
    "run the ODB API exposes modal frequency as cycles/time rounded to four "
    "decimal places, giving a first-mode relative quantization scale of about "
    "3.2e-5. The adopted relative frequency and period tolerance of 5e-5 and "
    "mode-shape component tolerance of 1e-4 account for that documented ABAQUS "
    "output precision while remaining strict enough to catch model or unit "
    "inconsistencies."
)


def load_output(name: str, path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"{name} output is missing: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{name} output is empty: {path}")
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    for key in ("circular_frequencies_rad_per_s", "periods_s", "normalized_mode_shapes"):
        if key not in data:
            raise KeyError(f"{name} output missing required key {key}")
    if len(data["circular_frequencies_rad_per_s"]) != 3:
        raise ValueError(f"{name} output does not contain 3 circular frequencies.")
    if len(data["periods_s"]) != 3:
        raise ValueError(f"{name} output does not contain 3 periods.")
    shapes = data["normalized_mode_shapes"]
    if len(shapes) != 3 or any(len(row) != 3 for row in shapes):
        raise ValueError(f"{name} normalized_mode_shapes must be a 3x3 floor-by-mode matrix.")
    return data


def as_float_list(values) -> list[float]:
    return [float(value) for value in values]


def shape_column(data: dict, mode_index: int) -> list[float]:
    return [float(row[mode_index]) for row in data["normalized_mode_shapes"]]


def normalize_for_compare(vector: list[float]) -> list[float]:
    scale = max(abs(value) for value in vector)
    if scale <= 0.0:
        raise ValueError("Cannot compare zero mode shape.")
    return [value / scale for value in vector]


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm2(a: list[float]) -> float:
    return math.sqrt(dot(a, a))


def compare_shapes(a: list[float], b: list[float]) -> dict:
    aa = normalize_for_compare(a)
    bb = normalize_for_compare(b)
    if dot(aa, bb) < 0.0:
        bb = [-value for value in bb]
    max_abs_component_difference = max(abs(x - y) for x, y in zip(aa, bb))
    denominator = dot(aa, aa) * dot(bb, bb)
    mac = (dot(aa, bb) ** 2) / denominator if denominator > 0.0 else 0.0
    return {
        "max_abs_component_difference": max_abs_component_difference,
        "mac": mac,
        "sign_corrected_reference": aa,
        "sign_corrected_candidate": bb,
    }


def compare_pair(name_a: str, data_a: dict, name_b: str, data_b: dict) -> dict:
    omega_a = as_float_list(data_a["circular_frequencies_rad_per_s"])
    omega_b = as_float_list(data_b["circular_frequencies_rad_per_s"])
    period_a = as_float_list(data_a["periods_s"])
    period_b = as_float_list(data_b["periods_s"])

    frequency_relative_differences = [
        abs(a - b) / max(abs(a), abs(b), 1.0e-300) for a, b in zip(omega_a, omega_b)
    ]
    period_relative_differences = [
        abs(a - b) / max(abs(a), abs(b), 1.0e-300) for a, b in zip(period_a, period_b)
    ]
    shape_comparisons = [
        compare_shapes(shape_column(data_a, mode_index), shape_column(data_b, mode_index))
        for mode_index in range(3)
    ]
    max_shape_component_difference = max(
        item["max_abs_component_difference"] for item in shape_comparisons
    )
    min_mac = min(item["mac"] for item in shape_comparisons)
    pair_pass = (
        max(frequency_relative_differences) <= THRESHOLDS["circular_frequency_relative_tolerance"]
        and max(period_relative_differences) <= THRESHOLDS["period_relative_tolerance"]
        and max_shape_component_difference <= THRESHOLDS["mode_shape_component_abs_tolerance"]
        and min_mac >= THRESHOLDS["mode_shape_mac_minimum"]
    )
    return {
        "pair": [name_a, name_b],
        "circular_frequency_relative_differences": frequency_relative_differences,
        "max_circular_frequency_relative_difference": max(frequency_relative_differences),
        "period_relative_differences": period_relative_differences,
        "max_period_relative_difference": max(period_relative_differences),
        "mode_shape_comparisons": shape_comparisons,
        "max_mode_shape_component_abs_difference": max_shape_component_difference,
        "min_mode_shape_mac": min_mac,
        "pass": pair_pass,
    }


def main() -> int:
    os.chdir(ROOT)
    loaded = {}
    errors = {}
    for name, path in SOFTWARE_OUTPUTS.items():
        try:
            loaded[name] = load_output(name, path)
        except Exception as exc:
            errors[name] = str(exc)

    pairwise = []
    if not errors and len(loaded) == 3:
        for name_a, name_b in combinations(SOFTWARE_OUTPUTS.keys(), 2):
            pairwise.append(compare_pair(name_a, loaded[name_a], name_b, loaded[name_b]))

    strict_pass = not errors and len(pairwise) == 3 and all(item["pass"] for item in pairwise)
    report = {
        "status": "PASS" if strict_pass else "FAIL",
        "strict_three_software_pass": strict_pass,
        "thresholds": THRESHOLDS,
        "threshold_reasoning_recorded_before_pass_fail": THRESHOLD_REASONING,
        "software_result_files": {name: str(path) for name, path in SOFTWARE_OUTPUTS.items()},
        "software_output_load_errors": errors,
        "pairwise_comparisons": pairwise,
    }
    with (ROOT / "cross_validation_report.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
        f.write("\n")

    write_markdown_report(report, loaded)
    print(json.dumps({"status": report["status"], "strict_three_software_pass": strict_pass}, indent=2))
    return 0 if strict_pass else 2


def write_markdown_report(report: dict, loaded: dict) -> None:
    lines = []
    lines.append("# Three-Story Shear Building Modal Verification")
    lines.append("")
    lines.append("## Model")
    lines.append("")
    lines.append("The verified model is a three-story linear shear building in SI units.")
    lines.append("Each floor has a 1000 kg lumped mass and each story has a 500000 N/m lateral stiffness.")
    lines.append("The floor DOF order is floor 1, floor 2, floor 3.")
    lines.append("")
    lines.append("The MATLAB mass and stiffness matrices are:")
    lines.append("")
    lines.append("```text")
    lines.append("M = 1000 * eye(3)")
    lines.append("K = 500000 * [ 2 -1  0; -1  2 -1; 0 -1 1 ]")
    lines.append("```")
    lines.append("")
    lines.append("Mode shapes in all JSON files are 3x3 matrices with rows as floors and columns as modes. Each mode is normalized to max(abs(component)) = 1 and signed so the roof component is positive.")
    lines.append("")
    lines.append("## Generated Outputs Used")
    lines.append("")
    for name, path in report["software_result_files"].items():
        lines.append(f"- {name}: `{Path(path).as_posix()}`")
    lines.append("- cross validation: `cross_validation_report.json`")
    lines.append("")
    lines.append("## Software Implementations")
    lines.append("")
    lines.append("- MATLAB builds the 3x3 mass and stiffness matrices directly and solves `K phi = lambda M phi`.")
    lines.append("- OpenSeesPy uses a 1D chain of floor nodes with lumped masses and zeroLength elastic story springs.")
    lines.append("- ABAQUS uses MASS elements at the floor nodes and SPRING2 elements in global U1, followed by Lanczos frequency extraction and ODB post-processing.")
    lines.append("")
    lines.append("## Consistency Criteria")
    lines.append("")
    lines.append(report["threshold_reasoning_recorded_before_pass_fail"])
    lines.append("")
    lines.append("Adopted thresholds:")
    for key, value in report["thresholds"].items():
        lines.append(f"- {key}: `{value}`")
    lines.append("")
    lines.append("## Errors Encountered and Corrections")
    lines.append("")
    lines.append("- The default Python environment was Python 3.13 while the installed `openseespywin` binary imports `python312.dll`. This caused an OpenSeesPy DLL import failure during environment probing. The actual OpenSeesPy analysis was run with the matching Python 3.12 environment at `<python_env>\\python.exe`; the probe is saved in `verification/openseespy/environment_probe.log`.")
    lines.append("- ABAQUS completed successfully and produced an ODB. The `.dat` file reports warnings that requested displacement normalization was replaced by mass normalization under the default SIM architecture and that some constrained DOFs were inactive; the post-processor renormalizes extracted U1 mode shapes, and the participating floor DOF remains the intended U1 shear-building DOF.")
    lines.append("- ABAQUS ODB frame frequencies are exposed as four-decimal cycles/time values in this run. The final tolerance was selected to account for that generated-output precision before computing the pass/fail result.")
    if report["software_output_load_errors"]:
        for name, error in report["software_output_load_errors"].items():
            lines.append(f"- {name} result load error: {error}")
    lines.append("")
    lines.append("## Modal Results")
    lines.append("")
    for name, data in loaded.items():
        lines.append(f"### {name}")
        lines.append("")
        lines.append("| Mode | Circular frequency (rad/s) | Period (s) | Floor 1 | Floor 2 | Floor 3 |")
        lines.append("|---:|---:|---:|---:|---:|---:|")
        omegas = data["circular_frequencies_rad_per_s"]
        periods = data["periods_s"]
        shapes = data["normalized_mode_shapes"]
        for mode_index in range(3):
            lines.append(
                f"| {mode_index + 1} | {float(omegas[mode_index]):.12g} | "
                f"{float(periods[mode_index]):.12g} | "
                f"{float(shapes[0][mode_index]):.12g} | "
                f"{float(shapes[1][mode_index]):.12g} | "
                f"{float(shapes[2][mode_index]):.12g} |"
            )
        lines.append("")
    lines.append("## Pairwise Metrics")
    lines.append("")
    if report["pairwise_comparisons"]:
        lines.append("| Pair | Max freq rel diff | Max period rel diff | Max mode component diff | Min MAC | Pass |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for item in report["pairwise_comparisons"]:
            pair_name = " vs ".join(item["pair"])
            lines.append(
                f"| {pair_name} | "
                f"{item['max_circular_frequency_relative_difference']:.6g} | "
                f"{item['max_period_relative_difference']:.6g} | "
                f"{item['max_mode_shape_component_abs_difference']:.6g} | "
                f"{item['min_mode_shape_mac']:.12g} | "
                f"{item['pass']} |"
            )
    else:
        lines.append("Pairwise metrics were not computed because one or more required software outputs were missing or invalid.")
    lines.append("")
    lines.append("## Final Agreement")
    lines.append("")
    if report["strict_three_software_pass"]:
        lines.append("Final status: strict three-software PASS. MATLAB, OpenSeesPy, and ABAQUS generated outputs are present and agree within the adopted thresholds.")
    else:
        lines.append("Final status: FAIL. A strict three-software pass is not declared because at least one required output is missing, invalid, or outside the adopted thresholds.")
    lines.append("")

    with (ROOT / "report.md").open("w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
