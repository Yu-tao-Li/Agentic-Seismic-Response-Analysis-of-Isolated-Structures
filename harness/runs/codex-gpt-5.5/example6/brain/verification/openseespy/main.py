#!/usr/bin/env python3
"""OpenSeesPy benchmark for Example 6.

The model is a 2D three-story, two-bay elastic steel frame with nonlinear
Bouc-Wen zeroLength isolators below the three column lines. The script writes
response_output.json plus CSV histories in this directory.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = SCRIPT_DIR.parents[1]
INPUT_PATH = RUN_ROOT / "input_data" / "Northridge_01_NO_968.txt"
OUTPUT_PATH = SCRIPT_DIR / "response_output.json"
HISTORY_CSV = SCRIPT_DIR / "time_histories.csv"
HYSTERESIS_CSV = SCRIPT_DIR / "isolator_hysteresis.csv"
NORMALIZED_GM = SCRIPT_DIR / "normalized_ground_motion.txt"


PARAMS = {
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
    "isolator_model": "BoucWen",
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
    "newmark_beta": 0.25,
    "newmark_gamma": 0.50,
}


def node_tag(level: int, col: int) -> int:
    if level == 0:
        return col + 1
    return level * 10 + col + 1


def ground_tag(col: int) -> int:
    return 100 + col + 1


def iso_ele_tag(col: int) -> int:
    return 1000 + col + 1


def read_and_normalize_motion() -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    if not INPUT_PATH.exists():
        write_unavailable_output(
            "missing_input",
            f"Required ground-motion file is missing: {INPUT_PATH}",
        )
        raise SystemExit(2)

    raw = np.loadtxt(INPUT_PATH, dtype=float).reshape(-1)
    if raw.size == 0:
        write_unavailable_output("empty_input", f"Ground-motion file is empty: {INPUT_PATH}")
        raise SystemExit(2)

    raw_peak = float(np.max(np.abs(raw)))
    if raw_peak <= 0.0:
        write_unavailable_output("invalid_input", "Ground-motion peak acceleration is zero.")
        raise SystemExit(2)

    target_pga = PARAMS["target_pga_g"] * PARAMS["g"]
    scale = target_pga / raw_peak
    ag = raw * scale
    np.savetxt(NORMALIZED_GM, ag, fmt="%.12e")
    time = np.arange(raw.size, dtype=float) * PARAMS["dt"]
    return raw, ag, {
        "dt_s": PARAMS["dt"],
        "sample_count": int(raw.size),
        "raw_peak_abs": raw_peak,
        "target_pga_m_per_s2": target_pga,
        "target_pga_g": PARAMS["target_pga_g"],
        "normalization_scale": scale,
        "normalized_peak_abs_m_per_s2": float(np.max(np.abs(ag))),
        "source_file": str(INPUT_PATH.relative_to(RUN_ROOT)),
    }


def import_opensees():
    try:
        import openseespy.opensees as ops  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on external install
        write_unavailable_output(
            "openseespy_unavailable",
            f"OpenSeesPy could not be imported: {type(exc).__name__}: {exc}",
        )
        raise SystemExit(2)
    return ops


def build_opensees_model(ops: Any) -> None:
    p = PARAMS
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    for col in range(3):
        x = col * p["bay_width_m"]
        ops.node(ground_tag(col), x, 0.0)
        ops.fix(ground_tag(col), 1, 1, 1)

    for level in range(4):
        y = level * p["story_height_m"]
        for col in range(3):
            x = col * p["bay_width_m"]
            ops.node(node_tag(level, col), x, y)

    m_node = p["floor_mass_kg"] / 3.0
    for level in range(1, 4):
        for col in range(3):
            ops.mass(node_tag(level, col), m_node, 0.0, 0.0)

    for level in range(1, 4):
        master = node_tag(level, 0)
        for col in (1, 2):
            ops.equalDOF(master, node_tag(level, col), 1)

    ops.geomTransf("Linear", 1)
    ele = 1
    for col in range(3):
        for level in range(3):
            ops.element(
                "elasticBeamColumn",
                ele,
                node_tag(level, col),
                node_tag(level + 1, col),
                p["column_area_m2"],
                p["E_Pa"],
                p["column_I_m4"],
                1,
            )
            ele += 1

    for level in range(1, 4):
        for col in range(2):
            ops.element(
                "elasticBeamColumn",
                ele,
                node_tag(level, col),
                node_tag(level, col + 1),
                p["beam_area_m2"],
                p["E_Pa"],
                p["beam_I_m4"],
                1,
            )
            ele += 1

    ops.uniaxialMaterial(
        "BoucWen",
        1,
        p["boucwen_alpha"],
        p["isolator_k0_N_per_m"],
        p["boucwen_n"],
        p["boucwen_gamma"],
        p["boucwen_beta"],
        p["boucwen_Ao"],
        p["boucwen_deltaA"],
        p["boucwen_deltaNu"],
        p["boucwen_deltaEta"],
    )
    ops.uniaxialMaterial("Elastic", 2, p["isolator_kv_N_per_m"])
    for col in range(3):
        ops.element(
            "zeroLength",
            iso_ele_tag(col),
            ground_tag(col),
            node_tag(0, col),
            "-mat",
            1,
            2,
            "-dir",
            1,
            2,
            "-doRayleigh",
            1,
        )


def frame_element_stiffness(E: float, A: float, I: float, xi: float, yi: float, xj: float, yj: float) -> np.ndarray:
    L = math.hypot(xj - xi, yj - yi)
    c = (xj - xi) / L
    s = (yj - yi) / L
    k = np.array(
        [
            [A * E / L, 0, 0, -A * E / L, 0, 0],
            [0, 12 * E * I / L**3, 6 * E * I / L**2, 0, -12 * E * I / L**3, 6 * E * I / L**2],
            [0, 6 * E * I / L**2, 4 * E * I / L, 0, -6 * E * I / L**2, 2 * E * I / L],
            [-A * E / L, 0, 0, A * E / L, 0, 0],
            [0, -12 * E * I / L**3, -6 * E * I / L**2, 0, 12 * E * I / L**3, -6 * E * I / L**2],
            [0, 6 * E * I / L**2, 2 * E * I / L, 0, -6 * E * I / L**2, 4 * E * I / L],
        ],
        dtype=float,
    )
    T = np.array(
        [
            [c, s, 0, 0, 0, 0],
            [-s, c, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, s, 0],
            [0, 0, 0, -s, c, 0],
            [0, 0, 0, 0, 0, 1],
        ],
        dtype=float,
    )
    return T.T @ k @ T


def reference_modal_properties() -> tuple[float, list[float], float, float]:
    """Independent condensation check used if OpenSees eigen extraction fails."""
    p = PARAMS
    ncols, nlev = 3, 4
    nn = ncols * nlev
    ndof = nn * 3

    def nid(level: int, col: int) -> int:
        return level * ncols + col

    def dofs(n: int) -> list[int]:
        return [3 * n, 3 * n + 1, 3 * n + 2]

    coords = [(col * p["bay_width_m"], level * p["story_height_m"]) for level in range(nlev) for col in range(ncols)]
    K = np.zeros((ndof, ndof), dtype=float)
    M = np.zeros((ndof, ndof), dtype=float)

    for col in range(ncols):
        for level in range(3):
            ni, nj = nid(level, col), nid(level + 1, col)
            ke = frame_element_stiffness(p["E_Pa"], p["column_area_m2"], p["column_I_m4"], *coords[ni], *coords[nj])
            ed = dofs(ni) + dofs(nj)
            K[np.ix_(ed, ed)] += ke

    for level in range(1, 4):
        for col in range(2):
            ni, nj = nid(level, col), nid(level, col + 1)
            ke = frame_element_stiffness(p["E_Pa"], p["beam_area_m2"], p["beam_I_m4"], *coords[ni], *coords[nj])
            ed = dofs(ni) + dofs(nj)
            K[np.ix_(ed, ed)] += ke

    for col in range(3):
        n = nid(0, col)
        K[3 * n, 3 * n] += p["isolator_k0_N_per_m"]
        K[3 * n + 1, 3 * n + 1] += p["isolator_kv_N_per_m"]

    for level in range(1, 4):
        for col in range(3):
            M[3 * nid(level, col), 3 * nid(level, col)] += p["floor_mass_kg"] / 3.0

    red_map: dict[tuple[Any, ...], int] = {}
    red_ids = []
    for level in range(nlev):
        for col in range(ncols):
            n = nid(level, col)
            for comp in range(3):
                gd = 3 * n + comp
                if comp == 0 and level >= 1:
                    key = ("floor_ux", level)
                else:
                    key = ("dof", gd)
                if key not in red_map:
                    red_map[key] = len(red_map)
                red_ids.append(red_map[key])

    T = np.zeros((ndof, len(red_map)), dtype=float)
    for gd, rd in enumerate(red_ids):
        T[gd, rd] = 1.0
    Kr = T.T @ K @ T
    Mr = T.T @ M @ T
    dyn = np.flatnonzero(np.diag(Mr) > 1e-9)
    stat = np.array([i for i in range(Kr.shape[0]) if i not in set(dyn)], dtype=int)
    Kdd = Kr[np.ix_(dyn, dyn)]
    Kds = Kr[np.ix_(dyn, stat)]
    Ksd = Kr[np.ix_(stat, dyn)]
    Kss = Kr[np.ix_(stat, stat)]
    Kcond = Kdd - Kds @ np.linalg.solve(Kss, Ksd)
    Mcond = Mr[np.ix_(dyn, dyn)]
    eigvals = np.linalg.eigvals(np.linalg.solve(Mcond, Kcond)).real
    eigvals = np.sort(eigvals[eigvals > 1e-8])
    omegas = np.sqrt(eigvals)
    periods = (2.0 * math.pi / omegas).tolist()
    alpha = 2.0 * p["damping_ratio"] * omegas[0] * omegas[1] / (omegas[0] + omegas[1])
    beta = 2.0 * p["damping_ratio"] / (omegas[0] + omegas[1])
    return periods[0], periods, float(alpha), float(beta)


def opensees_modal_properties(ops: Any) -> tuple[float, list[float], float, float, str]:
    try:
        lambdas = ops.eigen("-fullGenLapack", 2)
        omegas = [math.sqrt(float(x)) for x in lambdas if float(x) > 0.0]
        if len(omegas) < 2:
            raise RuntimeError(f"OpenSees returned fewer than two positive eigenvalues: {lambdas}")
        periods = [2.0 * math.pi / w for w in omegas]
        source = "OpenSees eigen('-fullGenLapack', 2)"
    except Exception as exc:
        first_period, periods, alpha, beta = reference_modal_properties()
        return first_period, periods, alpha, beta, f"reference condensation fallback because OpenSees eigen failed: {exc}"

    p = PARAMS
    alpha = 2.0 * p["damping_ratio"] * omegas[0] * omegas[1] / (omegas[0] + omegas[1])
    beta = 2.0 * p["damping_ratio"] / (omegas[0] + omegas[1])
    return periods[0], periods, float(alpha), float(beta), source


def value_from_ele_response(response: Any, default: float = 0.0) -> float:
    if response is None:
        return default
    if isinstance(response, (int, float)):
        return float(response)
    if len(response) == 0:
        return default
    return float(response[0])


def boucwen_z_from_force(force: float, deformation: float) -> float:
    p = PARAMS
    alpha = p["boucwen_alpha"]
    k0 = p["isolator_k0_N_per_m"]
    return (force - alpha * k0 * deformation) / ((1.0 - alpha) * k0)


def collect_response(ops: Any) -> dict[str, Any]:
    floor = [float(ops.nodeDisp(node_tag(level, 0), 1)) for level in (1, 2, 3)]
    base = [float(ops.nodeDisp(node_tag(0, col), 1)) for col in range(3)]
    iso_disp = []
    iso_force = []
    iso_z = []
    for col in range(3):
        ele = iso_ele_tag(col)
        strain = value_from_ele_response(ops.eleResponse(ele, "material", 1, "strain"), base[col])
        stress = value_from_ele_response(ops.eleResponse(ele, "material", 1, "stress"), 0.0)
        iso_disp.append(strain)
        iso_force.append(stress)
        iso_z.append(boucwen_z_from_force(stress, strain))
    h = PARAMS["story_height_m"]
    base_mean = float(np.mean(base))
    return {
        "roof": floor[2],
        "story_drifts": [(floor[0] - base_mean) / h, (floor[1] - floor[0]) / h, (floor[2] - floor[1]) / h],
        "isolation": float(np.mean(iso_disp)),
        "base_shear": float(np.sum(iso_force)),
        "iso_disp": iso_disp,
        "iso_force": iso_force,
        "iso_z": iso_z,
        "floor_disp": floor,
        "base_top_disp": base,
    }


def loop_area(force: np.ndarray, disp: np.ndarray) -> tuple[float, float]:
    if force.size < 2:
        return 0.0, 0.0
    inc_work = 0.5 * (force[1:] + force[:-1]) * np.diff(disp)
    return float(np.sum(inc_work)), float(np.sum(np.abs(inc_work)))


def peaks_from_histories(
    roof: np.ndarray,
    story: np.ndarray,
    iso: np.ndarray,
    base_shear: np.ndarray,
    iso_force: np.ndarray,
    iso_disp: np.ndarray,
) -> dict[str, Any]:
    areas = []
    signed = []
    for i in range(3):
        s, a = loop_area(iso_force[:, i], iso_disp[:, i])
        signed.append(s)
        areas.append(a)
    return {
        "peak_abs_roof_displacement_m": float(np.max(np.abs(roof))),
        "max_abs_story_drift_ratio": float(np.max(np.abs(story))),
        "peak_abs_isolation_displacement_m": float(np.max(np.abs(iso))),
        "peak_abs_total_base_shear_N": float(np.max(np.abs(base_shear))),
        "isolator_loop_area_signed_Nm": signed,
        "isolator_loop_area_abs_Nm": areas,
        "total_isolator_loop_area_signed_Nm": float(np.sum(signed)),
        "total_isolator_loop_area_abs_Nm": float(np.sum(areas)),
    }


def write_csvs(
    time: np.ndarray,
    ag: np.ndarray,
    roof: np.ndarray,
    story: np.ndarray,
    iso: np.ndarray,
    base_shear: np.ndarray,
    iso_disp: np.ndarray,
    iso_force: np.ndarray,
    iso_z: np.ndarray,
) -> None:
    with HISTORY_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "time_s",
                "ground_accel_m_per_s2",
                "roof_displacement_m",
                "story_1_drift_ratio",
                "story_2_drift_ratio",
                "story_3_drift_ratio",
                "isolation_displacement_m",
                "total_base_shear_N",
            ]
        )
        for i in range(time.size):
            writer.writerow([time[i], ag[i], roof[i], story[i, 0], story[i, 1], story[i, 2], iso[i], base_shear[i]])

    with HYSTERESIS_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "time_s",
                "iso1_disp_m",
                "iso1_force_N",
                "iso1_z_m",
                "iso2_disp_m",
                "iso2_force_N",
                "iso2_z_m",
                "iso3_disp_m",
                "iso3_force_N",
                "iso3_z_m",
            ]
        )
        for i in range(time.size):
            writer.writerow(
                [
                    time[i],
                    iso_disp[i, 0],
                    iso_force[i, 0],
                    iso_z[i, 0],
                    iso_disp[i, 1],
                    iso_force[i, 1],
                    iso_z[i, 1],
                    iso_disp[i, 2],
                    iso_force[i, 2],
                    iso_z[i, 2],
                ]
            )


def write_output(payload: dict[str, Any]) -> None:
    with OUTPUT_PATH.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def write_unavailable_output(status: str, message: str) -> None:
    payload = {
        "schema_version": "example6-response-v1",
        "software": "OpenSeesPy",
        "status": status,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "message": message,
        "assumptions": [
            "No engineering constants were invented; this diagnostic output is written only because an external runtime/input is unavailable.",
        ],
    }
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    write_output(payload)


def run() -> int:
    _, ag, input_meta = read_and_normalize_motion()
    ops = import_opensees()

    build_opensees_model(ops)
    fundamental_period, modal_periods, alpha_m, beta_k_init, period_source = opensees_modal_properties(ops)
    ops.rayleigh(alpha_m, 0.0, beta_k_init, 0.0)

    # Use value-backed PathSeries instead of file-backed PathSeries.  OpenSees
    # can fail to reopen files when the working path contains non-ASCII
    # characters, while the in-memory values follow the same normalized record.
    ops.timeSeries("Path", 1, "-dt", PARAMS["dt"], "-values", *[float(x) for x in ag.tolist()], "-factor", 1.0)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)

    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("UmfPack")
        linear_solver = "UmfPack"
    except Exception:
        ops.system("BandGeneral")
        linear_solver = "BandGeneral"
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", PARAMS["newmark_gamma"], PARAMS["newmark_beta"])
    ops.analysis("Transient")

    n = ag.size
    time = np.arange(n, dtype=float) * PARAMS["dt"]
    roof = np.zeros(n, dtype=float)
    story = np.zeros((n, 3), dtype=float)
    iso = np.zeros(n, dtype=float)
    base_shear = np.zeros(n, dtype=float)
    iso_disp = np.zeros((n, 3), dtype=float)
    iso_force = np.zeros((n, 3), dtype=float)
    iso_z = np.zeros((n, 3), dtype=float)
    iterations = np.zeros(n, dtype=int)
    retry_steps: list[int] = []
    failed_steps: list[int] = []

    r0 = collect_response(ops)
    roof[0] = r0["roof"]
    story[0, :] = r0["story_drifts"]
    iso[0] = r0["isolation"]
    base_shear[0] = r0["base_shear"]
    iso_disp[0, :] = r0["iso_disp"]
    iso_force[0, :] = r0["iso_force"]
    iso_z[0, :] = r0["iso_z"]

    for i in range(1, n):
        ok = ops.analyze(1, PARAMS["dt"])
        if ok != 0:
            retry_steps.append(i)
            ops.test("NormDispIncr", 1.0e-9, 100, 0)
            ops.algorithm("KrylovNewton")
            ok = ops.analyze(10, PARAMS["dt"] / 10.0)
            ops.test("NormDispIncr", 1.0e-10, 50, 0)
            ops.algorithm("Newton")
        if ok != 0:
            failed_steps.append(i)
            break

        try:
            iterations[i] = int(ops.testIter())
        except Exception:
            iterations[i] = -1
        r = collect_response(ops)
        roof[i] = r["roof"]
        story[i, :] = r["story_drifts"]
        iso[i] = r["isolation"]
        base_shear[i] = r["base_shear"]
        iso_disp[i, :] = r["iso_disp"]
        iso_force[i, :] = r["iso_force"]
        iso_z[i, :] = r["iso_z"]

    last = failed_steps[0] if failed_steps else n
    if last < n:
        time = time[:last]
        ag = ag[:last]
        roof = roof[:last]
        story = story[:last, :]
        iso = iso[:last]
        base_shear = base_shear[:last]
        iso_disp = iso_disp[:last, :]
        iso_force = iso_force[:last, :]
        iso_z = iso_z[:last, :]
        iterations = iterations[:last]

    write_csvs(time, ag, roof, story, iso, base_shear, iso_disp, iso_force, iso_z)
    peaks = peaks_from_histories(roof, story, iso, base_shear, iso_force, iso_disp)
    failed_count = len(failed_steps)
    attempted = max(0, n - 1)

    payload = {
        "schema_version": "example6-response-v1",
        "software": "OpenSeesPy",
        "status": "completed" if failed_count == 0 else "failed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "units": {"length": "m", "force": "N", "mass": "kg", "time": "s"},
        "input": input_meta,
        "model": PARAMS,
        "damping": {
            "type": "Rayleigh",
            "target_modes": [1, 2],
            "target_damping_ratio": PARAMS["damping_ratio"],
            "alpha_mass": alpha_m,
            "beta_stiffness_initial": beta_k_init,
            "period_source": period_source,
            "zeroLength_doRayleigh": True,
        },
        "fundamental_period_s": fundamental_period,
        "modal_periods_s": modal_periods[:6],
        "time_s": time.tolist(),
        "ground_accel_m_per_s2": ag.tolist(),
        "roof_displacement_m": roof.tolist(),
        "story_drift_ratios": {
            "story_1": story[:, 0].tolist(),
            "story_2": story[:, 1].tolist(),
            "story_3": story[:, 2].tolist(),
        },
        "isolation_displacement_m": iso.tolist(),
        "total_base_shear_N": base_shear.tolist(),
        "isolator_hysteresis": {
            f"isolator_{i + 1}": {
                "deformation_m": iso_disp[:, i].tolist(),
                "force_N": iso_force[:, i].tolist(),
                "boucwen_z_m": iso_z[:, i].tolist(),
                "loop_area_signed_Nm": peaks["isolator_loop_area_signed_Nm"][i],
                "loop_area_abs_Nm": peaks["isolator_loop_area_abs_Nm"][i],
            }
            for i in range(3)
        },
        "peaks": peaks,
        "convergence": {
            "analysis": "OpenSees Transient Newmark gamma=0.5 beta=0.25",
            "test": "NormDispIncr",
            "linear_solver": linear_solver,
            "attempted_steps": attempted,
            "completed_steps": int(time.size - 1),
            "retry_step_indices": retry_steps,
            "failed_step_indices": failed_steps,
            "failed_step_count": failed_count,
            "failed_step_rate": float(failed_count / max(1, attempted)),
            "iterations_per_completed_step": iterations.tolist(),
        },
        "assumptions": [
            "Floor horizontal diaphragms are enforced with equalDOF constraints at levels 1 through 3.",
            "Floor mass is assigned only to horizontal translational DOFs, one third of each floor mass at each floor node.",
            "Horizontal isolators use OpenSees uniaxialMaterial BoucWen with F = alpha*k0*u + (1-alpha)*k0*z; z is back-calculated from the reported material force for audit output.",
            "Story-1 drift uses the average top-of-isolator displacement as the lower reference; isolation displacement is the mean of the three isolator deformations.",
            "Total base shear is the sum of the three horizontal uniaxial-material forces, positive in the internal resisting-force convention.",
            "Gravity load and P-Delta effects are intentionally excluded.",
        ],
        "output_files": {
            "json": str(OUTPUT_PATH.relative_to(RUN_ROOT)),
            "time_histories_csv": str(HISTORY_CSV.relative_to(RUN_ROOT)),
            "isolator_hysteresis_csv": str(HYSTERESIS_CSV.relative_to(RUN_ROOT)),
            "normalized_ground_motion": str(NORMALIZED_GM.relative_to(RUN_ROOT)),
        },
    }
    write_output(payload)
    return 0 if failed_count == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
