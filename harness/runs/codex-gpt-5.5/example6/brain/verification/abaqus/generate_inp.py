#!/usr/bin/env python3
"""Generate the Abaqus input deck for Example 6.

The generated model follows the prompt literally: a 3-story, 2-bay elastic
beam-column frame with three direct Bouc-Wen UEL isolators at the base.
Equivalent horizontal inertia loads are used instead of base acceleration so
that the response coordinates match the MATLAB relative-coordinate benchmark.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = SCRIPT_DIR.parents[1]
INPUT_PATH = RUN_ROOT / "input_data" / "Northridge_01_NO_968.txt"
INP_PATH = SCRIPT_DIR / "example6_model.inp"
MOTION_PATH = SCRIPT_DIR / "normalized_ground_motion.txt"
METADATA_PATH = SCRIPT_DIR / "abaqus_model_metadata.json"


PARAMS: dict[str, Any] = {
    "g": 9.81,
    "target_pga_g": 0.40,
    "dt": 0.01,
    "stories": 3,
    "bays": 2,
    "bay_width_m": 6.0,
    "story_height_m": 3.6,
    "E_Pa": 2.06e11,
    "nu": 0.30,
    "rho_kg_per_m3": 7850.0,
    "column_area_m2": 0.020,
    "column_I_m4": 8.0e-4,
    "beam_area_m2": 0.015,
    "beam_I_m4": 4.5e-4,
    "floor_mass_kg": 2.5e5,
    "isolator_k0_N_per_m": 2.0e6,
    "isolator_characteristic_yield_force_N": 1.0e4,
    "isolator_characteristic_yield_displacement_m": 0.005,
    "boucwen_alpha": 0.05,
    "boucwen_n": 2.0,
    "boucwen_beta": 2.0e4,
    "boucwen_gamma": 2.0e4,
    "boucwen_Ao": 1.0,
    "boucwen_deltaA": 0.0,
    "boucwen_deltaNu": 0.0,
    "boucwen_deltaEta": 0.0,
    "isolator_kv_N_per_m": 1.0e10,
    "damping_ratio": 0.05,
    "rayleigh_alpha_mass": 0.18649697508516591,
    "rayleigh_beta_stiffness": 0.0070671988247071145,
    "abaqus_density_placeholder_kg_per_m3": 1.0e-9,
    "newmark_beta": 0.25,
    "newmark_gamma": 0.50,
}


def normalize_motion() -> tuple[np.ndarray, dict[str, Any]]:
    raw = np.loadtxt(INPUT_PATH, dtype=float).reshape(-1)
    if raw.size == 0:
        raise ValueError(f"Ground-motion file is empty: {INPUT_PATH}")
    raw_peak = float(np.max(np.abs(raw)))
    if raw_peak <= 0.0:
        raise ValueError("Ground-motion peak acceleration is zero.")
    target_pga = PARAMS["target_pga_g"] * PARAMS["g"]
    scale = target_pga / raw_peak
    ag = raw * scale
    meta = {
        "dt_s": PARAMS["dt"],
        "sample_count": int(raw.size),
        "raw_peak_abs": raw_peak,
        "target_pga_m_per_s2": target_pga,
        "target_pga_g": PARAMS["target_pga_g"],
        "normalization_scale": scale,
        "normalized_peak_abs_m_per_s2": float(np.max(np.abs(ag))),
        "source_file": "input_data/Northridge_01_NO_968.txt",
    }
    return ag, meta


def rect_dims(area: float, inertia: float) -> tuple[float, float]:
    """Return Abaqus rectangular beam-section dimensions matching A and I."""

    height = math.sqrt(12.0 * inertia / area)
    width = area / height
    return width, height


def node_label(level: int, col: int) -> int:
    if level == 0:
        return col + 1
    return level * 10 + col + 1


def write_amp(lines: list[str], ag: np.ndarray) -> None:
    lines.append("*AMPLITUDE, NAME=GM_ACCEL, DEFINITION=TABULAR, TIME=TOTAL TIME, SMOOTH=0.0")
    pairs = []
    dt = PARAMS["dt"]
    for i, value in enumerate(ag):
        pairs.append(f"{i * dt:.8f}, {value:.12e}")
    for i in range(0, len(pairs), 4):
        lines.append(", ".join(pairs[i : i + 4]))


def build_inp() -> dict[str, Any]:
    ag, input_meta = normalize_motion()
    np.savetxt(MOTION_PATH, ag, fmt="%.12e")

    p = PARAMS
    col_b, col_h = rect_dims(p["column_area_m2"], p["column_I_m4"])
    beam_b, beam_h = rect_dims(p["beam_area_m2"], p["beam_I_m4"])
    m_node = p["floor_mass_kg"] / 3.0
    alpha_m = p["rayleigh_alpha_mass"]
    beta_k = p["rayleigh_beta_stiffness"]
    floor_dashpot_c = alpha_m * m_node
    total_time = (ag.size - 1) * p["dt"]

    lines: list[str] = []
    lines.extend(
        [
            "*HEADING",
            "Example 6: 3-story 2-bay elastic frame with Bouc-Wen UEL isolators",
            "*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO",
            "*NODE",
        ]
    )

    for col in range(3):
        x = col * p["bay_width_m"]
        lines.append(f"{101 + col}, {x:.8f}, 0.0")
    for level in range(4):
        y = level * p["story_height_m"]
        for col in range(3):
            x = col * p["bay_width_m"]
            lines.append(f"{node_label(level, col)}, {x:.8f}, {y:.8f}")

    lines.extend(
        [
            "*NSET, NSET=GROUND",
            "101, 102, 103",
            "*NSET, NSET=BASE_TOP",
            "1, 2, 3",
            "*NSET, NSET=FLOOR1",
            "11, 12, 13",
            "*NSET, NSET=FLOOR2",
            "21, 22, 23",
            "*NSET, NSET=ROOF",
            "31, 32, 33",
            "*NSET, NSET=MASS_NODES",
            "11, 12, 13, 21, 22, 23, 31, 32, 33",
            "*ELEMENT, TYPE=B23, ELSET=COLUMNS",
        ]
    )
    etag = 1
    for col in range(3):
        for level in range(3):
            lines.append(f"{etag}, {node_label(level, col)}, {node_label(level + 1, col)}")
            etag += 1
    lines.append("*ELEMENT, TYPE=B23, ELSET=BEAMS")
    for level in range(1, 4):
        for col in range(2):
            lines.append(f"{etag}, {node_label(level, col)}, {node_label(level, col + 1)}")
            etag += 1

    lines.extend(
        [
            "*USER ELEMENT, TYPE=U100, NODES=2, COORDINATES=2, PROPERTIES=8, VARIABLES=6",
            "1, 2",
            "*ELEMENT, TYPE=U100, ELSET=ISO_H",
            "1001, 101, 1",
            "1002, 102, 2",
            "1003, 103, 3",
            "*UEL PROPERTY, ELSET=ISO_H",
            (
                f"{p['isolator_k0_N_per_m']:.12e}, {p['boucwen_alpha']:.12e}, "
                f"{p['boucwen_n']:.12e}, {p['boucwen_gamma']:.12e}, "
                f"{p['boucwen_beta']:.12e}, {p['boucwen_Ao']:.12e}, "
                f"{p['isolator_kv_N_per_m']:.12e}, {beta_k:.12e}"
            ),
            "*ELEMENT, TYPE=MASS, ELSET=FLOOR_MASS",
        ]
    )
    mass_tag = 2001
    for level in range(1, 4):
        for col in range(3):
            lines.append(f"{mass_tag}, {node_label(level, col)}")
            mass_tag += 1
    lines.extend(
        [
            "*MASS, ELSET=FLOOR_MASS, TYPE=ANISOTROPIC",
            f"{m_node:.12e}, 0.0, 0.0",
            "*ELEMENT, TYPE=DASHPOT1, ELSET=MASS_PROP_DASHPOTS",
        ]
    )
    dashpot_tag = 3001
    for level in range(1, 4):
        for col in range(3):
            lines.append(f"{dashpot_tag}, {node_label(level, col)}")
            dashpot_tag += 1
    lines.extend(
        [
            "*DASHPOT, ELSET=MASS_PROP_DASHPOTS",
            "1,",
            f"{floor_dashpot_c:.12e},",
            "*MATERIAL, NAME=STEEL_RAYLEIGH",
            "*DENSITY",
            f"{p['abaqus_density_placeholder_kg_per_m3']:.12e}",
            "*ELASTIC",
            f"{p['E_Pa']:.12e}, {p['nu']:.12e}",
            f"*DAMPING, ALPHA=0.0, BETA={beta_k:.12e}",
            "*BEAM SECTION, SECTION=RECT, MATERIAL=STEEL_RAYLEIGH, ELSET=COLUMNS",
            f"{col_b:.12e}, {col_h:.12e}",
            "0.0, 0.0, -1.0",
            "*BEAM SECTION, SECTION=RECT, MATERIAL=STEEL_RAYLEIGH, ELSET=BEAMS",
            f"{beam_b:.12e}, {beam_h:.12e}",
            "0.0, 0.0, -1.0",
        ]
    )

    for level in range(1, 4):
        master = node_label(level, 0)
        for slave_col in (1, 2):
            slave = node_label(level, slave_col)
            lines.extend(
                [
                    "*EQUATION",
                    "2",
                    f"{slave}, 1, 1.0, {master}, 1, -1.0",
                ]
            )

    lines.extend(
        [
            "*BOUNDARY",
            "GROUND, 1, 2",
        ]
    )
    write_amp(lines, ag)
    lines.extend(
        [
            "*STEP, NAME=DYNAMIC_GM, NLGEOM=NO, INC=200000",
            "*DYNAMIC, ALPHA=0.0",
            f"{0.5 * p['dt']:.8f}, {total_time:.8f}, 1.0e-7, {p['dt']:.8f}",
            "*CLOAD, AMPLITUDE=GM_ACCEL",
        ]
    )
    for level in range(1, 4):
        lines.append(f"{node_label(level, 0)}, 1, {-p['floor_mass_kg']:.12e}")
    lines.extend(
        [
            "*OUTPUT, FIELD, FREQUENCY=1",
            "*NODE OUTPUT, NSET=BASE_TOP",
            "U",
            "*NODE OUTPUT, NSET=FLOOR1",
            "U",
            "*NODE OUTPUT, NSET=FLOOR2",
            "U",
            "*NODE OUTPUT, NSET=ROOF",
            "U",
            "*NODE OUTPUT, NSET=GROUND",
            "U, RF",
            "*END STEP",
        ]
    )

    INP_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")
    metadata = {
        "input": input_meta,
        "model": PARAMS,
        "files": {
            "input_deck": INP_PATH.name,
            "uel": "boucwen_uel.for",
            "normalized_ground_motion": MOTION_PATH.name,
        },
        "sections": {
            "column_rect_width_m": col_b,
            "column_rect_height_m": col_h,
            "beam_rect_width_m": beam_b,
            "beam_rect_height_m": beam_h,
            "section_note": "Rectangular B23 sections are back-calculated to match the specified A and out-of-plane I.",
        },
        "damping": {
            "rayleigh_alpha_mass": alpha_m,
            "rayleigh_beta_stiffness": beta_k,
            "dashpot_coefficient_per_mass_node_Ns_per_m": floor_dashpot_c,
            "uel_uses_betaK_initial_stiffness_velocity_term": True,
        },
        "time_integration_note": (
            "Abaqus uses automatic time incrementation with maximum increment equal to the "
            "ground-motion sampling interval 0.01 s; smaller internal cutbacks are allowed "
            "only for nonlinear convergence."
        ),
        "mass_model_note": (
            "Abaqus/Standard requires a nonzero active material density in a general dynamic step. "
            "A tiny numerical placeholder density is therefore assigned to B23 frame elements; the "
            "benchmark inertia remains the prescribed horizontal anisotropic floor MASS elements."
        ),
        "abaqus_boucwen_implementation": {
            "direct_boucwen_state_evolution": True,
            "substituted_bilinear_model": False,
            "implementation": "Abaqus/Standard UEL in boucwen_uel.for",
        },
    }
    METADATA_PATH.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main() -> int:
    metadata = build_inp()
    print(json.dumps({"generated": str(INP_PATH), "metadata": metadata["input"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
