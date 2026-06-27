"""
Generate cross_validation_report.json by comparing modal results from
MATLAB, OpenSeesPy, and ABAQUS.

Consistency thresholds (based on modal-analysis practice):
  - Frequency: relative difference < 1% (well-posed linear eigenvalue problem
    should yield near-identical results across software)
  - Period: relative difference < 1%
  - Mode shape: MAC (Modal Assurance Criterion) > 0.99 for each mode pair,
    AND max absolute element-wise difference < 0.05 after sign correction
"""

import json
import math
import os
import sys

THRESH_FREQ_REL_DIFF = 0.01    # 1% relative difference
THRESH_PERIOD_REL_DIFF = 0.01  # 1% relative difference
THRESH_MAC = 0.99              # Modal Assurance Criterion
THRESH_SHAPE_MAX_DIFF = 0.05   # max absolute element-wise difference

SOFTWARE = ["MATLAB", "OpenSeesPy", "ABAQUS"]
FILES = {
    "MATLAB":     "verification/matlab/modal_results.json",
    "OpenSeesPy": "verification/openseespy/modal_results.json",
    "ABAQUS":     "verification/abaqus/modal_results.json",
}


def load_results(path):
    with open(path, 'r') as f:
        return json.load(f)


def compute_mac(phi1, phi2):
    """Modal Assurance Criterion between two mode shape vectors."""
    dot12 = sum(a * b for a, b in zip(phi1, phi2))
    dot11 = sum(a * a for a in phi1)
    dot22 = sum(b * b for b in phi2)
    denom = math.sqrt(dot11 * dot22)
    if denom < 1e-30:
        return 0.0
    return (dot12 ** 2) / (dot11 * dot22)


def sign_correct(phi_ref, phi):
    """Flip sign of phi if it has negative dot product with phi_ref."""
    dot = sum(a * b for a, b in zip(phi_ref, phi))
    if dot < 0:
        return [-x for x in phi]
    return phi


def rel_diff(a, b):
    if abs(a) < 1e-30 and abs(b) < 1e-30:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def main():
    # Load all results
    results = {}
    for sw in SOFTWARE:
        path = FILES[sw]
        if not os.path.exists(path):
            print(f"ERROR: {sw} results file not found at {path}")
            sys.exit(1)
        results[sw] = load_results(path)
        n_modes = len(results[sw]["circular_frequencies_rad_per_s"])
        print(f"Loaded {sw}: {n_modes} modes")

    nModes = len(results["MATLAB"]["circular_frequencies_rad_per_s"])

    report = {
        "thresholds": {
            "frequency_relative_difference": THRESH_FREQ_REL_DIFF,
            "period_relative_difference": THRESH_PERIOD_REL_DIFF,
            "mac_minimum": THRESH_MAC,
            "mode_shape_max_absolute_difference": THRESH_SHAPE_MAX_DIFF,
            "rationale": (
                "For a well-posed linear eigenvalue problem (3-DOF shear building), "
                "all three software packages should produce essentially identical results. "
                "1% frequency/period tolerance accounts for minor numerical precision differences. "
                "MAC > 0.99 and max element-wise difference < 0.05 confirm mode shape agreement."
            )
        },
        "software_comparison": {},
        "overall_pass": True
    }

    # Compare each pair of software
    pairs = [("MATLAB", "OpenSeesPy"), ("MATLAB", "ABAQUS"), ("OpenSeesPy", "ABAQUS")]

    for sw1, sw2 in pairs:
        pair_key = f"{sw1}_vs_{sw2}"
        pair_report = {
            "frequency_comparison": [],
            "period_comparison": [],
            "mode_shape_comparison": [],
            "pair_pass": True
        }

        for mode in range(nModes):
            # Frequency comparison
            w1 = results[sw1]["circular_frequencies_rad_per_s"][mode]
            w2 = results[sw2]["circular_frequencies_rad_per_s"][mode]
            freq_diff = rel_diff(w1, w2)
            freq_ok = freq_diff < THRESH_FREQ_REL_DIFF

            pair_report["frequency_comparison"].append({
                "mode": mode + 1,
                f"{sw1}_rad_s": w1,
                f"{sw2}_rad_s": w2,
                "relative_difference": freq_diff,
                "pass": freq_ok
            })

            # Period comparison
            T1 = results[sw1]["periods_s"][mode]
            T2 = results[sw2]["periods_s"][mode]
            per_diff = rel_diff(T1, T2)
            per_ok = per_diff < THRESH_PERIOD_REL_DIFF

            pair_report["period_comparison"].append({
                "mode": mode + 1,
                f"{sw1}_s": T1,
                f"{sw2}_s": T2,
                "relative_difference": per_diff,
                "pass": per_ok
            })

            # Mode shape comparison
            phi1 = results[sw1]["normalized_mode_shapes"][mode]
            phi2 = results[sw2]["normalized_mode_shapes"][mode]

            # Sign-correct phi2 relative to phi1
            phi2_corrected = sign_correct(phi1, phi2)

            mac = compute_mac(phi1, phi2_corrected)
            max_diff = max(abs(a - b) for a, b in zip(phi1, phi2_corrected))
            shape_ok = mac >= THRESH_MAC and max_diff < THRESH_SHAPE_MAX_DIFF

            pair_report["mode_shape_comparison"].append({
                "mode": mode + 1,
                f"{sw1}_shape": phi1,
                f"{sw2}_shape_corrected": phi2_corrected,
                "mac": mac,
                "max_absolute_difference": max_diff,
                "pass": shape_ok
            })

            if not (freq_ok and per_ok and shape_ok):
                pair_report["pair_pass"] = False

        report["software_comparison"][pair_key] = pair_report
        if not pair_report["pair_pass"]:
            report["overall_pass"] = False

    # Write report
    with open("cross_validation_report.json", "w") as f:
        json.dump(report, f, indent=2)

    # Print summary
    print("\n=== Cross-Validation Summary ===")
    for pair_key, pair_report in report["software_comparison"].items():
        status = "PASS" if pair_report["pair_pass"] else "FAIL"
        print(f"\n{pair_key}: {status}")
        for fc in pair_report["frequency_comparison"]:
            print(f"  Mode {fc['mode']}: freq rel_diff = {fc['relative_difference']:.2e} ({'OK' if fc['pass'] else 'FAIL'})")
        for sc in pair_report["mode_shape_comparison"]:
            print(f"  Mode {sc['mode']}: MAC = {sc['mac']:.6f}, max_diff = {sc['max_absolute_difference']:.2e} ({'OK' if sc['pass'] else 'FAIL'})")

    overall = "PASS" if report["overall_pass"] else "FAIL"
    print(f"\nOverall: {overall}")
    print("\nReport written to cross_validation_report.json")


if __name__ == "__main__":
    main()
