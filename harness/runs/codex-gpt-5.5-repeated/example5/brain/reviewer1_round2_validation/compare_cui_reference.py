"""Compare the corrected Example 5 result with Cui et al.'s external reference."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[5]
CUI_DIR = (
    REPO_ROOT
    / "harness"
    / "runs"
    / "codex-gpt-5.5"
    / "example5"
    / "brain"
    / "reference"
    / "origin_extracted"
    / "cui"
)


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
    if not np.allclose(generated[:, 0], reference[:, 0], rtol=0.0, atol=1.0e-12):
        raise ValueError("Generated and Cui reference time vectors do not match.")

    generated_values = generated[:, 1]
    reference_values = reference[:, 1]
    error = generated_values - reference_values
    reference_peak = float(np.max(np.abs(reference_values)))
    generated_peak = float(np.max(np.abs(generated_values)))
    response_range = float(np.ptp(reference_values))
    reference_norm = float(np.linalg.norm(reference_values))

    return {
        "samples": count,
        "generated_peak_abs": generated_peak,
        "reference_peak_abs": reference_peak,
        "peak_relative_error": abs(generated_peak - reference_peak)
        / max(reference_peak, np.finfo(float).eps),
        "nrmse": float(np.sqrt(np.mean(error**2)))
        / max(response_range, np.finfo(float).eps),
        "relative_l2_error": float(np.linalg.norm(error))
        / max(reference_norm, np.finfo(float).eps),
        "correlation": float(np.corrcoef(generated_values, reference_values)[0, 1]),
    }


def hysteresis_metrics(generated: np.ndarray, reference: np.ndarray) -> dict[str, float | int]:
    count = min(len(generated), len(reference))
    generated_disp = generated[:count, 1]
    generated_force = generated[:count, 2]
    reference_disp = reference[:count, 0]
    reference_force = reference[:count, 1]

    generated_area = float(abs(np.trapezoid(generated_force, generated_disp)))
    reference_area = float(abs(np.trapezoid(reference_force, reference_disp)))
    disp_error = generated_disp - reference_disp
    force_error = generated_force - reference_force

    return {
        "samples": count,
        "generated_peak_displacement_abs": float(np.max(np.abs(generated_disp))),
        "reference_peak_displacement_abs": float(np.max(np.abs(reference_disp))),
        "generated_peak_force_abs": float(np.max(np.abs(generated_force))),
        "reference_peak_force_abs": float(np.max(np.abs(reference_force))),
        "displacement_nrmse": float(np.sqrt(np.mean(disp_error**2)))
        / max(float(np.ptp(reference_disp)), np.finfo(float).eps),
        "force_nrmse": float(np.sqrt(np.mean(force_error**2)))
        / max(float(np.ptp(reference_force)), np.finfo(float).eps),
        "generated_loop_area_abs": generated_area,
        "reference_loop_area_abs": reference_area,
        "loop_area_relative_error": abs(generated_area - reference_area)
        / max(reference_area, np.finfo(float).eps),
    }


def main() -> None:
    generated_paths = {
        "top_displacement": HERE / "topStoDisIso2.txt",
        "top_acceleration": HERE / "topStoAccIso2.txt",
        "isolation_hysteresis": HERE / "hysteresis_layer1.txt",
    }
    reference_paths = {
        "top_displacement": CUI_DIR / "topStoDisIso2.txt",
        "top_acceleration": CUI_DIR / "topStoAccIso2.txt",
        "isolation_hysteresis": CUI_DIR / "hysteresis_layer1.txt",
    }
    for path in [*generated_paths.values(), *reference_paths.values()]:
        if not path.is_file():
            raise FileNotFoundError(path)

    displacement = response_metrics(
        np.loadtxt(generated_paths["top_displacement"]),
        np.loadtxt(reference_paths["top_displacement"]),
    )
    acceleration = response_metrics(
        np.loadtxt(generated_paths["top_acceleration"]),
        np.loadtxt(reference_paths["top_acceleration"]),
    )
    hysteresis = hysteresis_metrics(
        np.loadtxt(generated_paths["isolation_hysteresis"]),
        np.loadtxt(reference_paths["isolation_hysteresis"]),
    )

    result = {
        "external_reference": {
            "source": "Cui, Shen, and Yang (2022), Structural Seismic Response Analysis Programming and Application, China Architecture & Building Press",
            "origin": "Reference output generated outside the coding-agent harness and extracted from the archived Origin project",
            "comparison_window": "0.00-18.99 s, matching the 1900-sample generated record",
            "files": {
                name: {
                    "repository_relative_path": path.relative_to(REPO_ROOT).as_posix(),
                    "sha256": sha256(path),
                }
                for name, path in reference_paths.items()
            },
        },
        "metrics": {
            "top_displacement": displacement,
            "top_acceleration": acceleration,
            "isolation_hysteresis": hysteresis,
        },
    }
    (HERE / "cui_external_reference_comparison.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    rows = [
        ("top displacement peak relative error", displacement["peak_relative_error"]),
        ("top displacement NRMSE", displacement["nrmse"]),
        ("top displacement relative L2 error", displacement["relative_l2_error"]),
        ("top acceleration peak relative error", acceleration["peak_relative_error"]),
        ("top acceleration NRMSE", acceleration["nrmse"]),
        ("top acceleration relative L2 error", acceleration["relative_l2_error"]),
        ("hysteresis displacement NRMSE", hysteresis["displacement_nrmse"]),
        ("hysteresis force NRMSE", hysteresis["force_nrmse"]),
        ("hysteresis loop-area relative error", hysteresis["loop_area_relative_error"]),
    ]
    with (HERE / "cui_external_reference_comparison.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)

    report_lines = [
        "# External Independent Reference Comparison",
        "",
        "The corrected Example 5 MATLAB implementation was compared against the archived output from Cui, Shen, and Yang (2022). The Cui output predates and was generated outside the coding-agent harness. SHA-256 hashes identify the exact deposited reference files used in this comparison.",
        "",
        "The comparison uses the common 0.00-18.99 s window (1900 samples at 0.01 s). No sign alignment, amplitude rescaling, or fitted time shift was applied.",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
    ]
    report_lines.extend(f"| {name} | {value:.6g} |" for name, value in rows)
    report_lines.extend(
        [
            "",
            "Reference provenance and file hashes are recorded in `cui_external_reference_comparison.json`.",
        ]
    )
    (HERE / "cui_external_reference_comparison.md").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
