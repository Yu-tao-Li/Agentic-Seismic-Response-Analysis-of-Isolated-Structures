from __future__ import annotations

import json
import math
from pathlib import Path

import openseespy.opensees as ops


OUT = Path(__file__).resolve().parent


def normalize_modes(raw_modes: list[list[float]]) -> list[list[float]]:
    normalized = []
    for mode in raw_modes:
        peak = max(abs(v) for v in mode)
        if peak == 0:
            normalized.append(mode)
            continue
        row = [v / peak for v in mode]
        if row[-1] < 0:
            row = [-v for v in row]
        normalized.append(row)
    return normalized


def main() -> None:
    mass = 1000.0
    stiffness = 500000.0

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)

    # Node 0 is the fixed ground. Nodes 1-3 are the three floor DOFs.
    for node in range(4):
        ops.node(node, 0.0)
    ops.fix(0, 1)

    for node in range(1, 4):
        ops.mass(node, mass)

    ops.uniaxialMaterial("Elastic", 1, stiffness)
    for ele in range(1, 4):
        ops.element("zeroLength", ele, ele - 1, ele, "-mat", 1, "-dir", 1)

    eigenvalues = ops.eigen("-fullGenLapack", 3)
    frequencies = [math.sqrt(value) for value in eigenvalues]
    periods = [2.0 * math.pi / value for value in frequencies]

    raw_modes = []
    for mode in range(1, 4):
        raw_modes.append([ops.nodeEigenvector(node, mode, 1) for node in range(1, 4)])

    results = {
        "software": "OpenSeesPy",
        "model": "1D three-DOF shear building with zeroLength elastic story springs",
        "circular_frequencies_rad_per_s": frequencies,
        "periods_s": periods,
        "normalized_mode_shapes": normalize_modes(raw_modes),
        "raw_eigenvalues": list(eigenvalues),
        "raw_mode_shapes": raw_modes,
    }

    (OUT / "modal_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    ops.wipe()


if __name__ == "__main__":
    main()
