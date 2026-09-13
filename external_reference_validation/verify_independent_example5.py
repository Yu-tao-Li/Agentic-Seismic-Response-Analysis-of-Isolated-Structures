"""Verify the independently executed published-algorithm reference for Example 5.

The public files contain response histories only.  The copyrighted book listing
and page scans are not redistributed.  This script compares the deposited
independent series with the corrected harness MATLAB output using the response
definitions reported in the revised manuscript.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
EXAMPLE = ROOT / "example5_cui"
REFERENCE = EXAMPLE / "independent_published_algorithm_output"
HARNESS = EXAMPLE / "parsed_output"
OUTPUT = EXAMPLE / "independent_run_comparison.json"
EPS = np.finfo(float).eps


def load(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    return data.reshape(1, -1) if data.ndim == 1 else data


def compare_history(reference: np.ndarray, harness: np.ndarray) -> dict[str, float | int]:
    count = min(len(reference), len(harness))
    reference = reference[:count]
    harness = harness[:count]
    if not np.allclose(reference[:, 0], harness[:, 0], rtol=0.0, atol=1.0e-10):
        raise ValueError("The reference and harness time vectors do not agree.")
    rv = reference[:, 1]
    hv = harness[:, 1]
    error = hv - rv
    reference_peak = float(np.max(np.abs(rv)))
    harness_peak = float(np.max(np.abs(hv)))
    return {
        "samples": count,
        "time_start_s": round(float(reference[0, 0]), 10),
        "time_end_s": round(float(reference[-1, 0]), 10),
        "reference_peak_abs": reference_peak,
        "harness_peak_abs": harness_peak,
        "maximum_absolute_difference": float(np.max(np.abs(error))),
        "peak_relative_error": abs(harness_peak - reference_peak) / max(reference_peak, EPS),
        "nrmse": float(np.sqrt(np.mean(error**2))) / max(float(np.ptp(rv)), EPS),
    }


def compare_hysteresis(reference: np.ndarray, harness: np.ndarray) -> dict[str, float | int]:
    count = min(len(reference), len(harness))
    reference = reference[:count]
    harness = harness[:count]
    reference_x, reference_y = reference[:, 0], reference[:, 1]
    harness_x, harness_y = harness[:, 1], harness[:, 2]
    reference_area = float(abs(np.trapezoid(reference_y, reference_x)))
    harness_area = float(abs(np.trapezoid(harness_y, harness_x)))
    return {
        "samples": count,
        "reference_peak_displacement_abs": float(np.max(np.abs(reference_x))),
        "harness_peak_displacement_abs": float(np.max(np.abs(harness_x))),
        "reference_peak_force_abs": float(np.max(np.abs(reference_y))),
        "harness_peak_force_abs": float(np.max(np.abs(harness_y))),
        "maximum_absolute_displacement_difference": float(np.max(np.abs(harness_x - reference_x))),
        "maximum_absolute_force_difference": float(np.max(np.abs(harness_y - reference_y))),
        "displacement_nrmse": float(np.sqrt(np.mean((harness_x - reference_x) ** 2)))
        / max(float(np.ptp(reference_x)), EPS),
        "force_nrmse": float(np.sqrt(np.mean((harness_y - reference_y) ** 2)))
        / max(float(np.ptp(reference_y)), EPS),
        "reference_loop_area_abs": reference_area,
        "harness_loop_area_abs": harness_area,
        "loop_area_relative_error": abs(harness_area - reference_area) / max(reference_area, EPS),
    }


def main() -> None:
    displacement = compare_history(
        load(REFERENCE / "top_displacement_m.txt"),
        load(HARNESS / "canonical_matlab_top_displacement.txt"),
    )
    acceleration = compare_history(
        load(REFERENCE / "top_relative_acceleration_mps2.txt"),
        load(HARNESS / "canonical_matlab_top_relative_acceleration.txt"),
    )
    hysteresis = compare_hysteresis(
        load(REFERENCE / "isolator_hysteresis_m_N.txt"),
        load(HARNESS / "canonical_matlab_isolation_hysteresis.txt"),
    )
    result = {
        "external_reference": "Cui et al. published algorithm (independent run)",
        "source": {
            "citation": "J. Cui, X. Shen, and M. Yang, Structural Seismic Response Analysis Programming and Application, China Architecture & Building Press, 2022 (in Chinese)",
            "chapter": 11,
            "appendices": [1, 6],
            "execution_environment": "MATLAB R2025b",
            "source_listing_distributed": False,
        },
        "input_convention": {
            "ground_motion": "example5_cui/input/Northridge_01_NO_968.txt",
            "pga_scaling": "Normalized to 0.40 g to match the benchmark input convention.",
            "response_history_amplitude_rescaling": False,
            "time_shift": False,
            "sign_adjustment": False,
        },
        "comparison_target": "Corrected harness MATLAB output used in the revised manuscript.",
        "comparison_window": "0.00--18.99 s (1900 samples)",
        "reference_series": {
            "top_displacement": "example5_cui/independent_published_algorithm_output/top_displacement_m.txt",
            "top_relative_acceleration": "example5_cui/independent_published_algorithm_output/top_relative_acceleration_mps2.txt",
            "isolator_hysteresis": "example5_cui/independent_published_algorithm_output/isolator_hysteresis_m_N.txt",
        },
        "harness_series": {
            "top_displacement": "example5_cui/parsed_output/canonical_matlab_top_displacement.txt",
            "top_relative_acceleration": "example5_cui/parsed_output/canonical_matlab_top_relative_acceleration.txt",
            "isolator_hysteresis": "example5_cui/parsed_output/canonical_matlab_isolation_hysteresis.txt",
        },
        "metrics": {
            "top_displacement": displacement,
            "top_relative_acceleration": acceleration,
            "isolator_hysteresis": hysteresis,
            "maximum_peak_relative_error": max(
                displacement["peak_relative_error"],
                acceleration["peak_relative_error"],
                abs(hysteresis["harness_peak_force_abs"] - hysteresis["reference_peak_force_abs"])
                / max(hysteresis["reference_peak_force_abs"], EPS),
            ),
            "maximum_nrmse": max(
                displacement["nrmse"],
                acceleration["nrmse"],
                hysteresis["displacement_nrmse"],
                hysteresis["force_nrmse"],
            ),
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["metrics"], indent=2))


if __name__ == "__main__":
    main()
