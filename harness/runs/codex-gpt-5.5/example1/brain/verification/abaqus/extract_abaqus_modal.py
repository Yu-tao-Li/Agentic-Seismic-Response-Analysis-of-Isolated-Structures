from __future__ import annotations

import json
import math
from pathlib import Path

from odbAccess import openOdb


OUT = Path(__file__).resolve().parent
ODB = OUT / "three_story_modal.odb"


def normalize(mode):
    peak = max(abs(v) for v in mode)
    row = [v / peak for v in mode]
    if row[-1] < 0:
        row = [-v for v in row]
    return row


def main() -> None:
    odb = openOdb(path=str(ODB), readOnly=True)
    try:
        step = odb.steps["Frequency"]
        frequencies_hz = []
        periods = []
        circular = []
        modes = []

        for index in range(1, 4):
            frame = step.frames[index]
            freq_hz = float(frame.frequency)
            frequencies_hz.append(freq_hz)
            circular.append(2.0 * math.pi * freq_hz)
            periods.append(1.0 / freq_hz)

            displacement_field = frame.fieldOutputs["U"]
            by_label = {
                value.nodeLabel: float(value.data[1])
                for value in displacement_field.values
                if value.nodeLabel in (2, 3, 4)
            }
            modes.append(normalize([by_label[2], by_label[3], by_label[4]]))

        results = {
            "software": "ABAQUS",
            "model": "SPRING2/MASS equivalent three-story shear building",
            "frequencies_hz": frequencies_hz,
            "circular_frequencies_rad_per_s": circular,
            "periods_s": periods,
            "normalized_mode_shapes": modes,
        }

        (OUT / "modal_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    finally:
        odb.close()


if __name__ == "__main__":
    main()
