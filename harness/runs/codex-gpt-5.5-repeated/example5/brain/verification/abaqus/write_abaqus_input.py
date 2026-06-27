"""Generate the Abaqus connector model input file for the isolated 4-DOF building."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
GM_FILE = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"
INP_FILE = SCRIPT_DIR / "isolated_shear_abaqus.inp"


def assemble_story_matrix(story_stiffness: np.ndarray) -> np.ndarray:
    k_global = np.zeros((4, 4), dtype=float)
    for i, kval in enumerate(story_stiffness):
        e = np.zeros(4, dtype=float)
        e[i] = 1.0
        if i > 0:
            e[i - 1] = -1.0
        k_global += kval * np.outer(e, e)
    return k_global


def modal_rayleigh(masses: np.ndarray, stiffness: np.ndarray, damping_ratio: float) -> tuple[float, float, np.ndarray]:
    m_mat = np.diag(masses)
    k_mat = assemble_story_matrix(stiffness)
    eigvals = np.linalg.eigvals(np.linalg.solve(m_mat, k_mat)).real
    omega = np.sqrt(np.sort(eigvals[eigvals > 0.0]))
    mat = np.array([[1.0 / (2.0 * omega[0]), omega[0] / 2.0],
                    [1.0 / (2.0 * omega[3]), omega[3] / 2.0]])
    alpha_m, beta_k = np.linalg.solve(mat, np.array([damping_ratio, damping_ratio]))
    return float(alpha_m), float(beta_k), omega


def fmt(value: float) -> str:
    return f"{value:.12e}"


def connector_behavior(name: str, stiffness: float | None = None, yield_force: float | None = None,
                       post_ratio: float = 0.0, damping: float | None = None) -> list[str]:
    lines = [f"*CONNECTOR BEHAVIOR, NAME={name}"]
    if stiffness is not None:
        lines.extend(["*CONNECTOR ELASTICITY, COMPONENT=1", fmt(stiffness)])
    if yield_force is not None:
        if stiffness is None:
            raise ValueError("Plastic connector behavior requires an elastic stiffness.")
        hardening = 0.0 if post_ratio == 0.0 else post_ratio * stiffness / (1.0 - post_ratio)
        lines.extend([
            "*CONNECTOR PLASTICITY, COMPONENT=1",
            "*CONNECTOR HARDENING, TYPE=KINEMATIC, DEFINITION=PARAMETERS",
            f"{fmt(yield_force)}, {fmt(hardening)}, 0.0",
        ])
    if damping is not None and abs(damping) > 0.0:
        lines.extend(["*CONNECTOR DAMPING, COMPONENT=1", fmt(damping)])
    return lines


def main() -> None:
    dt = 0.01
    g0 = 9.80665
    target_pga = 0.40 * g0
    masses = np.array([2.50e5, 2.70e5, 2.70e5, 1.80e5], dtype=float)
    stiffness = np.array([50e6, 245e6, 195e6, 98e6], dtype=float)
    fy = np.array([500e3, 1225e3, 975e3, 490e3], dtype=float)
    b = np.array([0.05, 0.0, 0.0, 0.0], dtype=float)
    damping_ratio = 0.05

    raw = np.loadtxt(GM_FILE, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw)))
    if raw_peak <= 0.0:
        raise RuntimeError("Ground motion file has zero peak acceleration.")
    ag = raw * (target_pga / raw_peak)
    n = ag.size
    time = np.arange(n, dtype=float) * dt
    total_time = float(time[-1])

    alpha_m, beta_k, omega = modal_rayleigh(masses, stiffness, damping_ratio)
    story_damping = beta_k * stiffness
    mass_damping = alpha_m * masses

    lines: list[str] = []
    lines.extend([
        "*HEADING",
        "Four-DOF isolated shear building; SI units; relative-coordinate earthquake loading.",
        "*PREPRINT, ECHO=NO, MODEL=YES, HISTORY=NO, CONTACT=NO",
        "*NODE, NSET=ALLNODES",
        "1000, 0., 0., 0.",
        "1, 0., 0., 0.",
        "2, 0., 0., 0.",
        "3, 0., 0., 0.",
        "4, 0., 0., 0.",
        "*NSET, NSET=MASSNODES",
        "1, 2, 3, 4",
        "*NSET, NSET=BASE",
        "1000",
    ])
    for i, mass in enumerate(masses, start=1):
        lines.extend([
            f"*ELEMENT, TYPE=MASS, ELSET=MASS{i}",
            f"{100+i}, {i}",
            f"*MASS, ELSET=MASS{i}",
            fmt(float(mass)),
        ])

    lines.extend([
        "*ELEMENT, TYPE=CONN3D2, ELSET=E_STORY1",
        "1, 1000, 1",
        "*ELEMENT, TYPE=CONN3D2, ELSET=E_STORY2",
        "2, 1, 2",
        "*ELEMENT, TYPE=CONN3D2, ELSET=E_STORY3",
        "3, 2, 3",
        "*ELEMENT, TYPE=CONN3D2, ELSET=E_STORY4",
        "4, 3, 4",
        "*ELSET, ELSET=STORY_CONNECTORS",
        "1, 2, 3, 4",
    ])
    for i in range(1, 5):
        lines.extend([
            f"*ELEMENT, TYPE=CONN3D2, ELSET=E_MASSDAMP{i}",
            f"{20+i}, 1000, {i}",
        ])

    for i in range(4):
        lines.extend(connector_behavior(
            f"B_STORY{i+1}",
            stiffness=float(stiffness[i]),
            yield_force=float(fy[i]),
            post_ratio=float(b[i]),
            damping=float(story_damping[i]),
        ))
    for i in range(4):
        lines.extend(connector_behavior(f"B_MASSDAMP{i+1}", damping=float(mass_damping[i])))

    for i in range(1, 5):
        lines.extend([
            f"*CONNECTOR SECTION, ELSET=E_STORY{i}, BEHAVIOR=B_STORY{i}",
            "CARTESIAN",
        ])
    for i in range(1, 5):
        lines.extend([
            f"*CONNECTOR SECTION, ELSET=E_MASSDAMP{i}, BEHAVIOR=B_MASSDAMP{i}",
            "CARTESIAN",
        ])

    lines.extend([
        "*BOUNDARY",
        "BASE, 1, 6, 0.",
        "MASSNODES, 2, 6, 0.",
        "*AMPLITUDE, NAME=AG, TIME=TOTAL TIME",
    ])
    for tt, aa in zip(time, ag):
        lines.append(f"{fmt(float(tt))}, {fmt(float(aa))}")

    lines.extend([
        "*STEP, NAME=EARTHQUAKE, NLGEOM=NO, INC=100000",
        "*DYNAMIC, DIRECT",
        f"{fmt(dt)}, {fmt(total_time)}, {fmt(1.0e-10)}, {fmt(dt)}",
        "*CLOAD, AMPLITUDE=AG",
    ])
    for i, mass in enumerate(masses, start=1):
        lines.append(f"{i}, 1, {fmt(-float(mass))}")
    lines.extend([
        "*OUTPUT, FIELD, TIME INTERVAL=0.01",
        "*NODE OUTPUT, NSET=MASSNODES",
        "U, V, A",
        "*ELEMENT OUTPUT, ELSET=STORY_CONNECTORS",
        "CU, CUE, CEF, CTF",
        "*OUTPUT, HISTORY, TIME INTERVAL=0.01",
        "*NODE OUTPUT, NSET=MASSNODES",
        "U1, A1",
        "*ELEMENT OUTPUT, ELSET=STORY_CONNECTORS",
        "CU1, CEF1, CTF1",
        "*END STEP",
    ])

    INP_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    metadata = {
        "input_file": str(INP_FILE),
        "ground_motion_file": str(GM_FILE),
        "dt": dt,
        "sample_count": int(n),
        "raw_peak": raw_peak,
        "target_pga_mps2": target_pga,
        "normalized_peak_mps2": float(np.max(np.abs(ag))),
        "rayleigh_alphaM": alpha_m,
        "rayleigh_betaKinit": beta_k,
        "elastic_omega_rad_s": omega.tolist(),
        "story_connector_damping": story_damping.tolist(),
        "mass_connector_damping": mass_damping.tolist(),
    }
    (SCRIPT_DIR / "abaqus_input_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Wrote {INP_FILE}")
    print(f"Samples: {n}, total time: {total_time:.6f}, alphaM: {alpha_m:.12g}, betaKinit: {beta_k:.12g}")


if __name__ == "__main__":
    main()
