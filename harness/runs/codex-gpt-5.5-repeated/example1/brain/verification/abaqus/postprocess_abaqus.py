"""Extract modal frequencies and mode shapes from the Abaqus ODB."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

from odbAccess import openOdb


SCRIPT_DIR = Path(__file__).resolve().parent
FLOOR_NODE_LABELS = [2, 3, 4]


def normalize_mode(vector):
    scale = max(abs(value) for value in vector)
    if scale <= 0.0:
        raise ValueError("Zero mode shape encountered.")
    normalized = [float(value / scale) for value in vector]
    if normalized[-1] < 0.0:
        normalized = [-value for value in normalized]
    return normalized


def write_text_output(path, circular_frequencies, periods, mode_shapes_by_floor):
    with path.open("w", encoding="utf-8") as f:
        f.write("Three-story shear building modal analysis - ABAQUS\n")
        f.write("Units: kg, N, m, s\n")
        f.write("Story mass: 1000 kg\n")
        f.write("Story stiffness: 500000 N/m\n\n")
        f.write("Model: MASS elements at floor nodes and SPRING2 elements in global U1.\n\n")
        f.write("Mode results:\n")
        f.write(
            "Mode, circular_frequency_rad_per_s, period_s, "
            "normalized_floor_1, normalized_floor_2, normalized_floor_3\n"
        )
        for mode_id in range(3):
            f.write(
                f"{mode_id + 1}, "
                f"{circular_frequencies[mode_id]:.15g}, "
                f"{periods[mode_id]:.15g}, "
                f"{mode_shapes_by_floor[0][mode_id]:.15g}, "
                f"{mode_shapes_by_floor[1][mode_id]:.15g}, "
                f"{mode_shapes_by_floor[2][mode_id]:.15g}\n"
            )


def extract_results(odb_path):
    odb = openOdb(str(odb_path), readOnly=True)
    try:
        step = odb.steps["FREQUENCY"]
        mode_records = []
        for frame in step.frames:
            mode_number = getattr(frame, "mode", None)
            frequency_hz = getattr(frame, "frequency", None)
            if mode_number is None or mode_number <= 0:
                continue
            if frequency_hz is None or frequency_hz <= 0.0:
                continue
            if "U" not in frame.fieldOutputs:
                raise RuntimeError(f"Frame for mode {mode_number} has no U field output.")
            u_field = frame.fieldOutputs["U"]
            values_by_node = {}
            for value in u_field.values:
                if value.nodeLabel in FLOOR_NODE_LABELS:
                    values_by_node[value.nodeLabel] = float(value.data[0])
            missing = [node for node in FLOOR_NODE_LABELS if node not in values_by_node]
            if missing:
                raise RuntimeError(f"Missing U output for floor node labels {missing} in mode {mode_number}.")
            vector = [values_by_node[node] for node in FLOOR_NODE_LABELS]
            mode_records.append(
                {
                    "mode": int(mode_number),
                    "frequency_hz": float(frequency_hz),
                    "circular_frequency": 2.0 * math.pi * float(frequency_hz),
                    "period": 1.0 / float(frequency_hz),
                    "mode_shape": normalize_mode(vector),
                }
            )
        if len(mode_records) != 3:
            raise RuntimeError(f"Expected 3 positive frequency modes, found {len(mode_records)}.")
        mode_records.sort(key=lambda item: item["circular_frequency"])
        circular_frequencies = [record["circular_frequency"] for record in mode_records]
        periods = [record["period"] for record in mode_records]
        mode_vectors = [record["mode_shape"] for record in mode_records]
        mode_shapes_by_floor = [[mode_vectors[mode][floor] for mode in range(3)] for floor in range(3)]
        frequencies_hz = [record["frequency_hz"] for record in mode_records]
        return circular_frequencies, periods, mode_shapes_by_floor, frequencies_hz, mode_records
    finally:
        odb.close()


def main():
    odb_path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else SCRIPT_DIR / "three_story_shear_spring_mass.odb"
    circular_frequencies, periods, mode_shapes_by_floor, frequencies_hz, mode_records = extract_results(odb_path)

    results = {
        "software": "ABAQUS",
        "execution": {
            "source_input_file": str(SCRIPT_DIR / "three_story_shear_spring_mass.inp"),
            "postprocess_file": str(SCRIPT_DIR / "postprocess_abaqus.py"),
            "odb_file": str(odb_path),
        },
        "model": {
            "description": "Three-story shear building with equal lumped floor masses and equal story stiffnesses.",
            "units": "SI: kg, N, m, s",
            "story_mass_kg": 1000.0,
            "story_stiffness_N_per_m": 500000.0,
            "floor_node_labels": FLOOR_NODE_LABELS,
            "mass_element_type": "MASS",
            "spring_element_type": "SPRING2 in global U1",
        },
        "mode_shape_convention": (
            "Rows are floors 1..3; columns are modes 1..3; each mode is normalized "
            "to max(abs(component)) = 1 and signed so the roof component is positive."
        ),
        "frequencies_hz": frequencies_hz,
        "circular_frequencies_rad_per_s": circular_frequencies,
        "periods_s": periods,
        "normalized_mode_shapes": mode_shapes_by_floor,
        "raw_abaqus_mode_records": mode_records,
    }

    with (SCRIPT_DIR / "modal_results.json").open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        f.write("\n")
    write_text_output(SCRIPT_DIR / "modal_results.txt", circular_frequencies, periods, mode_shapes_by_floor)
    print("ABAQUS ODB post-processing completed.")
    print(SCRIPT_DIR / "modal_results.json")
    return 0


if __name__ == "__main__":
    os.chdir(SCRIPT_DIR)
    raise SystemExit(main())
