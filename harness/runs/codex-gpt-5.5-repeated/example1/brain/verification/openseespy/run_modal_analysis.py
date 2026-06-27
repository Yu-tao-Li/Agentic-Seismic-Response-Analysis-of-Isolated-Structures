"""Modal analysis of a three-story shear building with OpenSeesPy.

Rows in mode shape output are floors 1..3; columns are modes 1..3.
The model uses one horizontal translational DOF per floor, lumped floor
masses, and zeroLength elastic story springs in series.
"""

from __future__ import annotations

import json
import math
import os
import platform
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def normalize_mode(vector: list[float]) -> list[float]:
    scale = max(abs(value) for value in vector)
    if scale <= 0.0:
        raise ValueError("Zero mode shape encountered.")
    normalized = [float(value / scale) for value in vector]
    if normalized[-1] < 0.0:
        normalized = [-value for value in normalized]
    return normalized


def write_text_output(
    path: Path,
    story_mass: float,
    story_stiffness: float,
    circular_frequencies: list[float],
    periods: list[float],
    mode_shapes_by_floor: list[list[float]],
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("Three-story shear building modal analysis - OpenSeesPy\n")
        f.write("Units: kg, N, m, s\n")
        f.write(f"Story mass: {story_mass:.12g} kg\n")
        f.write(f"Story stiffness: {story_stiffness:.12g} N/m\n\n")
        f.write("Model: 1D floor nodes with lumped masses and zeroLength Elastic story springs.\n\n")
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


def main() -> int:
    import openseespy.opensees as ops

    story_mass = 1000.0
    story_stiffness = 500000.0

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)

    base_node = 0
    floor_nodes = [1, 2, 3]
    ops.node(base_node, 0.0)
    ops.fix(base_node, 1)

    for node in floor_nodes:
        ops.node(node, 0.0)
        ops.mass(node, story_mass)

    mat_tag = 1
    ops.uniaxialMaterial("Elastic", mat_tag, story_stiffness)
    ops.element("zeroLength", 1, 0, 1, "-mat", mat_tag, "-dir", 1)
    ops.element("zeroLength", 2, 1, 2, "-mat", mat_tag, "-dir", 1)
    ops.element("zeroLength", 3, 2, 3, "-mat", mat_tag, "-dir", 1)

    raw_eigenvalues = [float(value) for value in ops.eigen("-fullGenLapack", 3)]
    ordered_modes = sorted(enumerate(raw_eigenvalues, start=1), key=lambda item: item[1])

    circular_frequencies: list[float] = []
    periods: list[float] = []
    mode_vectors: list[list[float]] = []
    for original_mode, eigenvalue in ordered_modes:
        if eigenvalue <= 0.0:
            raise ValueError(f"Non-positive eigenvalue for mode {original_mode}: {eigenvalue}")
        circular_frequency = math.sqrt(eigenvalue)
        circular_frequencies.append(circular_frequency)
        periods.append(2.0 * math.pi / circular_frequency)
        vector = [float(ops.nodeEigenvector(node, original_mode, 1)) for node in floor_nodes]
        mode_vectors.append(normalize_mode(vector))

    mode_shapes_by_floor = [[mode_vectors[mode][floor] for mode in range(3)] for floor in range(3)]

    results = {
        "software": "OpenSeesPy",
        "execution": {
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "opensees_version": getattr(ops, "version", lambda: "unknown")(),
            "source_file": str(SCRIPT_DIR / "run_modal_analysis.py"),
        },
        "model": {
            "description": "Three-story shear building with equal lumped floor masses and equal story stiffnesses.",
            "units": "SI: kg, N, m, s",
            "story_mass_kg": story_mass,
            "story_stiffness_N_per_m": story_stiffness,
            "base_node": base_node,
            "floor_nodes": floor_nodes,
            "element_type": "zeroLength with Elastic uniaxial material",
        },
        "mode_shape_convention": (
            "Rows are floors 1..3; columns are modes 1..3; each mode is normalized "
            "to max(abs(component)) = 1 and signed so the roof component is positive."
        ),
        "circular_frequencies_rad_per_s": circular_frequencies,
        "periods_s": periods,
        "normalized_mode_shapes": mode_shapes_by_floor,
    }

    with (SCRIPT_DIR / "modal_results.json").open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        f.write("\n")

    write_text_output(
        SCRIPT_DIR / "modal_results.txt",
        story_mass,
        story_stiffness,
        circular_frequencies,
        periods,
        mode_shapes_by_floor,
    )

    ops.wipe()
    print("OpenSeesPy modal analysis completed.")
    print(SCRIPT_DIR / "modal_results.json")
    return 0


if __name__ == "__main__":
    os.chdir(SCRIPT_DIR)
    raise SystemExit(main())
