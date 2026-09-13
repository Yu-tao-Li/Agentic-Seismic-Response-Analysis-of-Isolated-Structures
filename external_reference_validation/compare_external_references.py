"""Build the cross-example external-reference evidence package.

The script uses only author-confirmed existing outputs.  It does not execute or
modify either external program and does not fit a time shift or response scale.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
EPS = np.finfo(float).eps


def loadtxt(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    return data.reshape(1, -1) if data.ndim == 1 else data


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def response_metrics(generated: np.ndarray, reference: np.ndarray) -> dict[str, float | int]:
    count = min(len(generated), len(reference))
    generated = generated[:count]
    reference = reference[:count]
    if not np.allclose(generated[:, 0], reference[:, 0], rtol=0.0, atol=1.0e-10):
        raise ValueError("Time vectors do not agree over the common sampled interval.")
    gv = generated[:, 1]
    rv = reference[:, 1]
    error = gv - rv
    generated_peak = float(np.max(np.abs(gv)))
    reference_peak = float(np.max(np.abs(rv)))
    return {
        "samples": count,
        "time_start_s": float(generated[0, 0]),
        "time_end_s": float(generated[-1, 0]),
        "generated_peak_abs": generated_peak,
        "reference_peak_abs": reference_peak,
        "peak_relative_error": abs(generated_peak - reference_peak) / max(reference_peak, EPS),
        "nrmse": float(np.sqrt(np.mean(error**2))) / max(float(np.ptp(rv)), EPS),
        "relative_l2_error": float(np.linalg.norm(error)) / max(float(np.linalg.norm(rv)), EPS),
        "correlation": float(np.corrcoef(gv, rv)[0, 1]),
    }


def index_response_metrics(generated_values: np.ndarray, reference_values: np.ndarray) -> dict[str, float | int]:
    count = min(len(generated_values), len(reference_values))
    gv = generated_values[:count]
    rv = reference_values[:count]
    error = gv - rv
    generated_peak = float(np.max(np.abs(gv)))
    reference_peak = float(np.max(np.abs(rv)))
    return {
        "samples": count,
        "generated_peak_abs": generated_peak,
        "reference_peak_abs": reference_peak,
        "peak_relative_error": abs(generated_peak - reference_peak) / max(reference_peak, EPS),
        "nrmse": float(np.sqrt(np.mean(error**2))) / max(float(np.ptp(rv)), EPS),
        "relative_l2_error": float(np.linalg.norm(error)) / max(float(np.linalg.norm(rv)), EPS),
        "correlation": float(np.corrcoef(gv, rv)[0, 1]),
    }


def hysteresis_metrics(generated: np.ndarray, reference: np.ndarray) -> dict[str, float | int]:
    count = min(len(generated), len(reference))
    gx = generated[:count, 1]
    gy = generated[:count, 2]
    rx = reference[:count, 0]
    ry = reference[:count, 1]
    generated_area = float(abs(np.trapezoid(gy, gx)))
    reference_area = float(abs(np.trapezoid(ry, rx)))
    generated_peak_force = float(np.max(np.abs(gy)))
    reference_peak_force = float(np.max(np.abs(ry)))
    return {
        "samples": count,
        "generated_peak_displacement_abs": float(np.max(np.abs(gx))),
        "reference_peak_displacement_abs": float(np.max(np.abs(rx))),
        "generated_peak_force_abs": generated_peak_force,
        "reference_peak_force_abs": reference_peak_force,
        "peak_force_relative_error": abs(generated_peak_force - reference_peak_force)
        / max(reference_peak_force, EPS),
        "displacement_nrmse": float(np.sqrt(np.mean((gx - rx) ** 2))) / max(float(np.ptp(rx)), EPS),
        "force_nrmse": float(np.sqrt(np.mean((gy - ry) ** 2))) / max(float(np.ptp(ry)), EPS),
        "generated_loop_area_abs": generated_area,
        "reference_loop_area_abs": reference_area,
        "loop_area_relative_error": abs(generated_area - reference_area) / max(reference_area, EPS),
    }


def modal_metrics() -> dict[str, Any]:
    generated_path = ROOT / "example1_chen" / "parsed_output" / "canonical_matlab_modal_results.json"
    reference_path = ROOT / "example1_chen" / "raw_output" / "reference_metrics.json"
    generated = json.loads(generated_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    gw = np.asarray(generated["circular_frequencies_rad_per_s"], dtype=float)
    rw = np.asarray(reference["circular_frequencies_rad_per_s"], dtype=float)
    gt = np.asarray(generated["periods_s"], dtype=float)
    rt = np.asarray(reference["periods_s"], dtype=float)
    gm = np.asarray(generated["normalized_mode_shapes"], dtype=float).T
    rm = np.asarray(reference["normalized_mode_shapes"], dtype=float)
    mac_values: list[float] = []
    component_differences: list[float] = []
    for generated_mode, reference_mode in zip(gm, rm, strict=True):
        if np.dot(generated_mode, reference_mode) < 0.0:
            reference_mode = -reference_mode
        mac_values.append(
            float(
                np.dot(generated_mode, reference_mode) ** 2
                / (np.dot(generated_mode, generated_mode) * np.dot(reference_mode, reference_mode))
            )
        )
        component_differences.append(float(np.max(np.abs(generated_mode - reference_mode))))
    return {
        "generated_file": str(generated_path.relative_to(ROOT).as_posix()),
        "reference_file": str(reference_path.relative_to(ROOT).as_posix()),
        "max_frequency_relative_error": float(np.max(np.abs(gw - rw) / np.abs(rw))),
        "max_period_relative_error": float(np.max(np.abs(gt - rt) / np.abs(rt))),
        "minimum_mac": min(mac_values),
        "max_sign_aligned_component_abs_difference": max(component_differences),
        "mode_shape_sign_rule": "Whole-mode sign alignment only; no component-wise sign fitting.",
    }


def prepare_example4_relative_acceleration() -> np.ndarray:
    response_path = ROOT / "example4_cui" / "parsed_output" / "canonical_matlab_response.json"
    response = json.loads(response_path.read_text(encoding="utf-8"))
    series = np.column_stack(
        (
            np.asarray(response["time"], dtype=float),
            np.asarray(response["top_story_relative_acceleration"], dtype=float),
        )
    )
    output_path = ROOT / "example4_cui" / "parsed_output" / "canonical_matlab_top_relative_acceleration.txt"
    np.savetxt(output_path, series, fmt="%.12e")
    return series


def build_metrics() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    ex2_generated = loadtxt(
        ROOT / "example2_chen" / "parsed_output" / "canonical_matlab_top_displacement.txt"
    )
    ex2_reference = loadtxt(ROOT / "example2_chen" / "raw_output" / "chen_top_displacement.txt")
    ex4_generated_displacement = loadtxt(
        ROOT / "example4_cui" / "parsed_output" / "canonical_matlab_top_displacement.txt"
    )
    ex4_reference_displacement = loadtxt(
        ROOT / "example4_cui" / "original_external_output" / "cui_top_displacement.txt"
    )
    ex4_generated_acceleration = prepare_example4_relative_acceleration()
    ex4_reference_acceleration = loadtxt(
        ROOT / "example4_cui" / "original_external_output" / "cui_top_relative_acceleration.txt"
    )
    ex4_generated_isolation = loadtxt(
        ROOT / "example4_cui" / "parsed_output" / "canonical_matlab_isolation_displacement.txt"
    )
    ex4_reference_isolation_force = loadtxt(
        ROOT / "example4_cui" / "original_external_output" / "cui_isolation_displacement_force.txt"
    )
    ex5_generated_displacement = loadtxt(
        ROOT / "example5_cui" / "parsed_output" / "canonical_matlab_top_displacement.txt"
    )
    ex5_reference_displacement = loadtxt(
        ROOT / "example5_cui" / "original_external_output" / "cui_top_displacement.txt"
    )
    ex5_generated_acceleration = loadtxt(
        ROOT / "example5_cui" / "parsed_output" / "canonical_matlab_top_relative_acceleration.txt"
    )
    ex5_reference_acceleration = loadtxt(
        ROOT / "example5_cui" / "original_external_output" / "cui_top_relative_acceleration.txt"
    )
    ex5_generated_hysteresis = loadtxt(
        ROOT / "example5_cui" / "parsed_output" / "canonical_matlab_isolation_hysteresis.txt"
    )
    ex5_reference_hysteresis = loadtxt(
        ROOT / "example5_cui" / "original_external_output" / "cui_isolation_hysteresis.txt"
    )

    metrics = {
        "comparison_rules": {
            "source": "Existing author-confirmed external outputs; no external solver was rerun for this revision package.",
            "time_history_window": "Common sampled interval, truncated only at the shorter record end.",
            "response_amplitude_rescaling": False,
            "fitted_time_shift": False,
            "response_sign_alignment": False,
            "mode_shape_sign_alignment": True,
            "nrmse_definition": "RMSE divided by the range of the external-reference series.",
        },
        "example1": modal_metrics(),
        "example2": {
            "external_reference": "Chen compiled MDOF program (article id=408)",
            "top_displacement": response_metrics(ex2_generated, ex2_reference),
        },
        "example3": {"external_reference": None, "status": "internal_only"},
        "example4": {
            "external_reference": "Cui et al. published MATLAB implementation",
            "top_displacement": response_metrics(ex4_generated_displacement, ex4_reference_displacement),
            "top_relative_acceleration": response_metrics(
                ex4_generated_acceleration, ex4_reference_acceleration
            ),
            "isolation_displacement": index_response_metrics(
                ex4_generated_isolation[:, 1], ex4_reference_isolation_force[:, 0]
            ),
            "force_comparison_note": "The archived Cui et al. reference force series includes a different force-output convention; it is retained as raw evidence but excluded from the summary metric.",
        },
        "example5": {
            "external_reference": "Cui et al. published MATLAB implementation",
            "top_displacement": response_metrics(ex5_generated_displacement, ex5_reference_displacement),
            "top_relative_acceleration": response_metrics(
                ex5_generated_acceleration, ex5_reference_acceleration
            ),
            "isolation_hysteresis": hysteresis_metrics(
                ex5_generated_hysteresis, ex5_reference_hysteresis
            ),
        },
        "example6": {"external_reference": None, "status": "internal_only"},
    }
    series = {
        "ex2_generated": ex2_generated,
        "ex2_reference": ex2_reference,
        "ex4_generated_displacement": ex4_generated_displacement,
        "ex4_reference_displacement": ex4_reference_displacement,
        "ex5_generated_displacement": ex5_generated_displacement,
        "ex5_reference_displacement": ex5_reference_displacement,
        "ex5_generated_hysteresis": ex5_generated_hysteresis,
        "ex5_reference_hysteresis": ex5_reference_hysteresis,
    }
    return metrics, series


def write_comparisons(metrics: dict[str, Any]) -> None:
    for example, folder in ((1, "example1_chen"), (2, "example2_chen"), (4, "example4_cui"), (5, "example5_cui")):
        path = ROOT / folder / "comparison.json"
        path.write_text(json.dumps(metrics[f"example{example}"], indent=2) + "\n", encoding="utf-8")


def build_summary(metrics: dict[str, Any]) -> list[dict[str, str]]:
    ex1 = metrics["example1"]
    ex2 = metrics["example2"]["top_displacement"]
    ex4_items = [
        metrics["example4"]["top_displacement"],
        metrics["example4"]["top_relative_acceleration"],
        metrics["example4"]["isolation_displacement"],
    ]
    ex5_disp = metrics["example5"]["top_displacement"]
    ex5_acc = metrics["example5"]["top_relative_acceleration"]
    ex5_hys = metrics["example5"]["isolation_hysteresis"]
    return [
        {
            "example": "1",
            "external_reference": "Chen compiled program",
            "responses_checked": "Frequencies; periods; mode shapes",
            "peak_or_modal_metric": f"max frequency error={ex1['max_frequency_relative_error']:.3e}; max period error={ex1['max_period_relative_error']:.3e}",
            "curve_or_shape_metric": f"min MAC={ex1['minimum_mac']:.9f}",
            "result": "Consistent",
        },
        {
            "example": "2",
            "external_reference": "Chen compiled program",
            "responses_checked": "Top displacement",
            "peak_or_modal_metric": f"peak error={ex2['peak_relative_error']:.3e}",
            "curve_or_shape_metric": f"NRMSE={ex2['nrmse']:.3e}",
            "result": "Consistent",
        },
        {
            "example": "3",
            "external_reference": "Not available",
            "responses_checked": "-",
            "peak_or_modal_metric": "-",
            "curve_or_shape_metric": "-",
            "result": "Internal only",
        },
        {
            "example": "4",
            "external_reference": "Cui et al. published implementation",
            "responses_checked": "Top displacement; top relative acceleration; isolation displacement",
            "peak_or_modal_metric": f"max peak error={max(item['peak_relative_error'] for item in ex4_items):.3e}",
            "curve_or_shape_metric": f"max NRMSE={max(item['nrmse'] for item in ex4_items):.3e}",
            "result": "Consistent",
        },
        {
            "example": "5",
            "external_reference": "Cui et al. published implementation",
            "responses_checked": "Top response; isolation displacement and force; hysteresis",
            "peak_or_modal_metric": (
                f"max peak error={max(ex5_disp['peak_relative_error'], ex5_acc['peak_relative_error'], ex5_hys['peak_force_relative_error']):.3e}; "
                f"loop-area error={ex5_hys['loop_area_relative_error']:.3e}"
            ),
            "curve_or_shape_metric": (
                f"max NRMSE={max(ex5_disp['nrmse'], ex5_acc['nrmse'], ex5_hys['displacement_nrmse'], ex5_hys['force_nrmse']):.3e}"
            ),
            "result": "Consistent",
        },
        {
            "example": "6",
            "external_reference": "Not available",
            "responses_checked": "-",
            "peak_or_modal_metric": "-",
            "curve_or_shape_metric": "-",
            "result": "Internal only",
        },
    ]


def write_summary(metrics: dict[str, Any], rows: list[dict[str, str]]) -> None:
    (ROOT / "external_reference_summary.json").write_text(
        json.dumps({"summary": rows, "detailed_metrics": metrics}, indent=2) + "\n",
        encoding="utf-8",
    )
    with (ROOT / "external_reference_summary.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = [
        "# External reference report",
        "",
        "This report consolidates existing author-confirmed results. No Chen program or implementation associated with Cui et al. was rerun or modified during this packaging step.",
        "",
        "| Example | External reference | Responses checked | Peak/modal metric | Curve/shape metric | Result |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        report.append(
            f"| {row['example']} | {row['external_reference']} | {row['responses_checked']} | "
            f"{row['peak_or_modal_metric']} | {row['curve_or_shape_metric']} | {row['result']} |"
        )
    report.extend(
        [
            "",
            "## Interpretation notes",
            "",
            "- Example 1 uses whole-mode sign alignment only because an eigenvector's global sign is arbitrary.",
            "- Example 2 uses the common 0.00-39.98 s interval because the archived Chen output contains 2000 samples.",
            "- The archived Cui et al. reference acceleration quantity for Example 4 is relative acceleration, so it is compared with the matching corrected harness MATLAB quantity.",
            "- The archived Example 4 force series uses a different force-output convention and is retained but not summarized as a force discrepancy.",
            "- Example 5 compares top response, isolation displacement and restoring force, and hysteresis-loop area over 0.00-18.99 s.",
            "- No response amplitude scaling, fitted time shift, or response sign adjustment was applied.",
        ]
    )
    (ROOT / "external_reference_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def plot_figure(series: dict[str, np.ndarray]) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman"],
            "font.size": 9.0,
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 4.4))
    black = "#111111"
    blue = "#1f77b4"
    red = "#d62728"

    panels = [
        (axes[0, 0], series["ex2_generated"], series["ex2_reference"], "Top displacement (m)", "Chen", blue, "--"),
        (
            axes[0, 1],
            series["ex4_generated_displacement"],
            series["ex4_reference_displacement"],
            "Top displacement (m)",
            "Cui et al.",
            red,
            "--",
        ),
        (
            axes[1, 0],
            series["ex5_generated_displacement"],
            series["ex5_reference_displacement"],
            "Top displacement (m)",
            "Cui et al.",
            red,
            "--",
        ),
    ]
    legend_handles: dict[str, Any] = {}
    for label, (ax, generated, reference, ylabel, reference_name, reference_color, reference_style) in zip("abc", panels, strict=True):
        count = min(len(generated), len(reference))
        corrected_line = ax.plot(generated[:count, 0], generated[:count, 1], color=black, lw=1.25, zorder=2, label="Corrected MATLAB")[0]
        reference_line = ax.plot(reference[:count, 0], reference[:count, 1], color=reference_color, lw=0.90, ls=reference_style, zorder=3, label=reference_name)[0]
        legend_handles.setdefault("Corrected MATLAB", corrected_line)
        legend_handles.setdefault(reference_name, reference_line)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(ylabel)
        ax.text(0.015, 0.98, f"({label})", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(top=False, right=False, width=0.8, length=3.5)

    ax = axes[1, 1]
    generated = series["ex5_generated_hysteresis"]
    reference = series["ex5_reference_hysteresis"]
    count = min(len(generated), len(reference))
    corrected_line = ax.plot(generated[:count, 1], generated[:count, 2] / 1000.0, color=black, lw=1.25, zorder=2, label="Corrected MATLAB")[0]
    cui_line = ax.plot(reference[:count, 0], reference[:count, 1] / 1000.0, color=red, lw=0.90, ls="--", zorder=3, label="Cui et al.")[0]
    legend_handles.setdefault("Corrected MATLAB", corrected_line)
    legend_handles.setdefault("Cui et al.", cui_line)
    ax.set_xlabel("Isolation displacement (m)")
    ax.set_ylabel("Restoring force (kN)")
    ax.text(0.015, 0.98, "(d)", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False, width=0.8, length=3.5)

    legend_order = ["Corrected MATLAB", "Chen", "Cui et al."]
    fig.legend(
        [legend_handles[name] for name in legend_order],
        legend_order,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
        handlelength=2.2,
        columnspacing=1.0,
        handletextpad=0.35,
    )
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.12, top=0.90, wspace=0.28, hspace=0.42)

    figure_dir = ROOT / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "external_reference_comparison.png", dpi=300, bbox_inches="tight")
    fig.savefig(figure_dir / "external_reference_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def write_manifest() -> None:
    entries = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or "private_programs" in path.parts:
            continue
        relative = path.relative_to(ROOT).as_posix()
        if relative == "provenance/sha256_manifest.txt":
            continue
        entries.append(f"{sha256(path)}  {relative}")
    (ROOT / "provenance" / "sha256_manifest.txt").write_text("\n".join(entries) + "\n", encoding="utf-8")


def main() -> None:
    metrics, series = build_metrics()
    write_comparisons(metrics)
    rows = build_summary(metrics)
    write_summary(metrics, rows)
    plot_figure(series)
    write_manifest()
    print(json.dumps({"status": "ok", "summary": rows}, indent=2))


if __name__ == "__main__":
    main()
