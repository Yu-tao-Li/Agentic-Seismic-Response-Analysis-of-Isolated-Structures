"""Build the provenance and SHA-256 table for the external-reference series."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "provenance" / "external_reference_provenance_sha256.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


ROWS = [
    {
        "example": "1",
        "series": "modal frequencies, periods, and mode shapes",
        "reference_role": "archived external reference",
        "external_reference_file": "example1_chen/raw_output/reference_metrics.json",
        "harness_output_file": "example1_chen/parsed_output/canonical_matlab_modal_results.json",
        "provenance": "Chen compiled MDOF program; author-confirmed archived output.",
    },
    {
        "example": "2",
        "series": "top displacement",
        "reference_role": "archived external reference",
        "external_reference_file": "example2_chen/raw_output/chen_top_displacement.txt",
        "harness_output_file": "example2_chen/parsed_output/canonical_matlab_top_displacement.txt",
        "provenance": "Chen compiled MDOF program; author-confirmed archived output.",
    },
    {
        "example": "4",
        "series": "top displacement",
        "reference_role": "archived external reference",
        "external_reference_file": "example4_cui/original_external_output/cui_top_displacement.txt",
        "harness_output_file": "example4_cui/parsed_output/canonical_matlab_top_displacement.txt",
        "provenance": "Cui et al. published MATLAB implementation; archived output retained as an immutable comparison input.",
    },
    {
        "example": "4",
        "series": "top relative acceleration",
        "reference_role": "archived external reference",
        "external_reference_file": "example4_cui/original_external_output/cui_top_relative_acceleration.txt",
        "harness_output_file": "example4_cui/parsed_output/canonical_matlab_top_relative_acceleration.txt",
        "provenance": "Cui et al. published MATLAB implementation; matched relative-acceleration definition.",
    },
    {
        "example": "4",
        "series": "isolation displacement",
        "reference_role": "archived external reference",
        "external_reference_file": "example4_cui/original_external_output/cui_isolation_displacement_force.txt",
        "harness_output_file": "example4_cui/parsed_output/canonical_matlab_isolation_displacement.txt",
        "provenance": "Cui et al. published MATLAB implementation; displacement column of the archived displacement-force output.",
    },
    {
        "example": "5",
        "series": "top displacement",
        "reference_role": "historical deposited comparison record",
        "external_reference_file": "example5_cui/original_external_output/cui_top_displacement.txt",
        "harness_output_file": "example5_cui/parsed_output/canonical_matlab_top_displacement.txt",
        "provenance": "Previously deposited archived series; retained unchanged as an integrity record.",
    },
    {
        "example": "5",
        "series": "top relative acceleration",
        "reference_role": "historical deposited comparison record",
        "external_reference_file": "example5_cui/original_external_output/cui_top_relative_acceleration.txt",
        "harness_output_file": "example5_cui/parsed_output/canonical_matlab_top_relative_acceleration.txt",
        "provenance": "Previously deposited archived series; retained unchanged as an integrity record.",
    },
    {
        "example": "5",
        "series": "isolator hysteresis displacement and force",
        "reference_role": "historical deposited comparison record",
        "external_reference_file": "example5_cui/original_external_output/cui_isolation_hysteresis.txt",
        "harness_output_file": "example5_cui/parsed_output/canonical_matlab_isolation_hysteresis.txt",
        "provenance": "Previously deposited archived series; retained unchanged as an integrity record.",
    },
    {
        "example": "5",
        "series": "top displacement",
        "reference_role": "independently executed published algorithm",
        "external_reference_file": "example5_cui/independent_published_algorithm_output/top_displacement_m.txt",
        "harness_output_file": "example5_cui/parsed_output/canonical_matlab_top_displacement.txt",
        "provenance": "Cui et al. Chapter 11 algorithm retranscribed and executed independently in MATLAB R2025b; PGA normalized to 0.40 g.",
    },
    {
        "example": "5",
        "series": "top relative acceleration",
        "reference_role": "independently executed published algorithm",
        "external_reference_file": "example5_cui/independent_published_algorithm_output/top_relative_acceleration_mps2.txt",
        "harness_output_file": "example5_cui/parsed_output/canonical_matlab_top_relative_acceleration.txt",
        "provenance": "Cui et al. Chapter 11 algorithm retranscribed and executed independently in MATLAB R2025b; PGA normalized to 0.40 g.",
    },
    {
        "example": "5",
        "series": "isolator hysteresis displacement and force",
        "reference_role": "independently executed published algorithm",
        "external_reference_file": "example5_cui/independent_published_algorithm_output/isolator_hysteresis_m_N.txt",
        "harness_output_file": "example5_cui/parsed_output/canonical_matlab_isolation_hysteresis.txt",
        "provenance": "Cui et al. Chapter 11 algorithm retranscribed and executed independently in MATLAB R2025b; PGA normalized to 0.40 g.",
    },
]


def main() -> None:
    fieldnames = [
        "example",
        "series",
        "reference_role",
        "external_reference_file",
        "external_reference_sha256",
        "harness_output_file",
        "harness_output_sha256",
        "provenance",
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in ROWS:
            external = ROOT / row["external_reference_file"]
            harness = ROOT / row["harness_output_file"]
            if not external.is_file() or not harness.is_file():
                raise FileNotFoundError(f"Missing manifest input: {external} or {harness}")
            writer.writerow(
                {
                    **row,
                    "external_reference_sha256": sha256(external),
                    "harness_output_sha256": sha256(harness),
                }
            )
    print(OUTPUT)


if __name__ == "__main__":
    main()
