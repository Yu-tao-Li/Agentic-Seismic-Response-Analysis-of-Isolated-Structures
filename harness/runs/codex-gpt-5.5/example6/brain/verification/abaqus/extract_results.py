#!/usr/bin/env python3
"""Extract and package Abaqus Example 6 results.

Run with `abaqus python extract_results.py` after `example6_model.odb` is
created. The ODB supplies displacements; isolator forces are recomputed from
the same Bouc-Wen backward-Euler state law implemented in `boucwen_uel.for`,
because Abaqus does not reliably expose custom UEL force components as standard
ODB stress output.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = SCRIPT_DIR.parents[1]
INPUT_PATH = RUN_ROOT / "input_data" / "Northridge_01_NO_968.txt"
OUTPUT_PATH = SCRIPT_DIR / "response_output.json"
HISTORY_CSV = SCRIPT_DIR / "abaqus_history.csv"
HYSTERESIS_CSV = SCRIPT_DIR / "isolator_hysteresis.csv"
ODB_PATH = SCRIPT_DIR / "example6_model.odb"
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


def frame_element_stiffness(E: float, area: float, inertia: float, xi: float, yi: float, xj: float, yj: float) -> np.ndarray:
    length = math.hypot(xj - xi, yj - yi)
    c = (xj - xi) / length
    s = (yj - yi) / length
    k = np.array(
        [
            [area * E / length, 0, 0, -area * E / length, 0, 0],
            [0, 12 * E * inertia / length**3, 6 * E * inertia / length**2, 0, -12 * E * inertia / length**3, 6 * E * inertia / length**2],
            [0, 6 * E * inertia / length**2, 4 * E * inertia / length, 0, -6 * E * inertia / length**2, 2 * E * inertia / length],
            [-area * E / length, 0, 0, area * E / length, 0, 0],
            [0, -12 * E * inertia / length**3, -6 * E * inertia / length**2, 0, 12 * E * inertia / length**3, -6 * E * inertia / length**2],
            [0, 6 * E * inertia / length**2, 2 * E * inertia / length, 0, -6 * E * inertia / length**2, 4 * E * inertia / length],
        ],
        dtype=float,
    )
    trans = np.array(
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
    return trans.T @ k @ trans


def reference_modal_properties() -> tuple[float, list[float], float, float]:
    p = PARAMS
    ncols, nlev = 3, 4
    nn = ncols * nlev
    ndof = nn * 3

    def nid(level: int, col: int) -> int:
        return level * ncols + col

    def dofs(n: int) -> list[int]:
        return [3 * n, 3 * n + 1, 3 * n + 2]

    coords = [(col * p["bay_width_m"], level * p["story_height_m"]) for level in range(nlev) for col in range(ncols)]
    stiffness = np.zeros((ndof, ndof), dtype=float)
    mass = np.zeros((ndof, ndof), dtype=float)

    for col in range(ncols):
        for level in range(3):
            ni, nj = nid(level, col), nid(level + 1, col)
            ke = frame_element_stiffness(p["E_Pa"], p["column_area_m2"], p["column_I_m4"], *coords[ni], *coords[nj])
            ed = dofs(ni) + dofs(nj)
            stiffness[np.ix_(ed, ed)] += ke

    for level in range(1, 4):
        for col in range(2):
            ni, nj = nid(level, col), nid(level, col + 1)
            ke = frame_element_stiffness(p["E_Pa"], p["beam_area_m2"], p["beam_I_m4"], *coords[ni], *coords[nj])
            ed = dofs(ni) + dofs(nj)
            stiffness[np.ix_(ed, ed)] += ke

    for col in range(3):
        n = nid(0, col)
        stiffness[3 * n, 3 * n] += p["isolator_k0_N_per_m"]
        stiffness[3 * n + 1, 3 * n + 1] += p["isolator_kv_N_per_m"]

    for level in range(1, 4):
        for col in range(3):
            mass[3 * nid(level, col), 3 * nid(level, col)] += p["floor_mass_kg"] / 3.0

    red_map: dict[tuple[Any, ...], int] = {}
    red_ids = []
    for level in range(nlev):
        for col in range(ncols):
            n = nid(level, col)
            for comp in range(3):
                gd = 3 * n + comp
                key = ("floor_ux", level) if comp == 0 and level >= 1 else ("dof", gd)
                red_map.setdefault(key, len(red_map))
                red_ids.append(red_map[key])

    transform = np.zeros((ndof, len(red_map)), dtype=float)
    for gd, rd in enumerate(red_ids):
        transform[gd, rd] = 1.0
    kr = transform.T @ stiffness @ transform
    mr = transform.T @ mass @ transform
    dyn = np.flatnonzero(np.diag(mr) > 1e-9)
    stat = np.array([i for i in range(kr.shape[0]) if i not in set(dyn)], dtype=int)
    kdd = kr[np.ix_(dyn, dyn)]
    kds = kr[np.ix_(dyn, stat)]
    ksd = kr[np.ix_(stat, dyn)]
    kss = kr[np.ix_(stat, stat)]
    kcond = kdd - kds @ np.linalg.solve(kss, ksd)
    mcond = mr[np.ix_(dyn, dyn)]
    eigvals = np.linalg.eigvals(np.linalg.solve(mcond, kcond)).real
    eigvals = np.sort(eigvals[eigvals > 1e-8])
    omegas = np.sqrt(eigvals)
    periods = (2.0 * math.pi / omegas).tolist()
    alpha = 2.0 * p["damping_ratio"] * omegas[0] * omegas[1] / (omegas[0] + omegas[1])
    beta = 2.0 * p["damping_ratio"] / (omegas[0] + omegas[1])
    return periods[0], periods, float(alpha), float(beta)


def read_and_normalize_motion() -> tuple[np.ndarray, dict[str, Any]]:
    raw = np.loadtxt(INPUT_PATH, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw)))
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


def sign_nonzero(value: float) -> float:
    return 1.0 if value > 0.0 else -1.0


def boucwen_update(u_new: float, u_old: float, z_old: float) -> tuple[float, float]:
    p = PARAMS
    du = u_new - u_old
    z = z_old
    for _ in range(40):
        abs_z = max(abs(z), 1.0e-14)
        psi = p["boucwen_gamma"] + p["boucwen_beta"] * sign_nonzero(du * z)
        phi = p["boucwen_Ao"] - abs_z ** p["boucwen_n"] * psi
        residual = z - z_old - du * phi
        dphi_dz = -p["boucwen_n"] * abs_z ** (p["boucwen_n"] - 1.0) * sign_nonzero(z) * psi
        dg_dz = 1.0 - du * dphi_dz
        dz = -residual / dg_dz
        z += dz
        if abs(dz) <= 1.0e-12 * max(1.0, abs(z)):
            break
    abs_z = max(abs(z), 1.0e-14)
    psi = p["boucwen_gamma"] + p["boucwen_beta"] * sign_nonzero(du * z)
    phi = p["boucwen_Ao"] - abs_z ** p["boucwen_n"] * psi
    dphi_dz = -p["boucwen_n"] * abs_z ** (p["boucwen_n"] - 1.0) * sign_nonzero(z) * psi
    dg_dz = 1.0 - du * dphi_dz
    dz_du = phi / dg_dz
    return z, dz_du


def boucwen_histories(disp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p = PARAMS
    force = np.zeros_like(disp, dtype=float)
    zhist = np.zeros_like(disp, dtype=float)
    for j in range(disp.shape[1]):
        u_old = 0.0
        z_old = 0.0
        for i in range(disp.shape[0]):
            if i == 0:
                u_old = disp[i, j]
                zhist[i, j] = z_old
            else:
                z_old, _ = boucwen_update(float(disp[i, j]), float(u_old), float(z_old))
                u_old = float(disp[i, j])
                zhist[i, j] = z_old
            force[i, j] = p["boucwen_alpha"] * p["isolator_k0_N_per_m"] * disp[i, j] + (
                1.0 - p["boucwen_alpha"]
            ) * p["isolator_k0_N_per_m"] * zhist[i, j]
    return force, zhist


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
    signed = []
    abs_area = []
    for i in range(3):
        s, a = loop_area(iso_force[:, i], iso_disp[:, i])
        signed.append(s)
        abs_area.append(a)
    return {
        "peak_abs_roof_displacement_m": float(np.max(np.abs(roof))),
        "max_abs_story_drift_ratio": float(np.max(np.abs(story))),
        "peak_abs_isolation_displacement_m": float(np.max(np.abs(iso))),
        "peak_abs_total_base_shear_N": float(np.max(np.abs(base_shear))),
        "isolator_loop_area_signed_Nm": signed,
        "isolator_loop_area_abs_Nm": abs_area,
        "total_isolator_loop_area_signed_Nm": float(np.sum(signed)),
        "total_isolator_loop_area_abs_Nm": float(np.sum(abs_area)),
    }


def export_history_from_odb() -> None:
    try:
        from odbAccess import openOdb  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "ODB extraction requires Abaqus Python. Run `abaqus python extract_results.py`."
        ) from exc

    odb = openOdb(str(ODB_PATH))
    try:
        step = odb.steps["DYNAMIC_GM"] if "DYNAMIC_GM" in odb.steps else next(iter(odb.steps.values()))
        labels = [101, 102, 103, 1, 2, 3, 11, 12, 13, 21, 22, 23, 31, 32, 33]
        rows = []
        for frame in step.frames:
            u = {}
            for value in frame.fieldOutputs["U"].values:
                if value.nodeLabel in labels:
                    u[int(value.nodeLabel)] = float(value.data[0])
            missing = [label for label in labels if label not in u]
            if missing:
                raise KeyError(f"Missing U output for node labels {missing} at time {frame.frameValue}")

            base = np.array([u[1], u[2], u[3]], dtype=float)
            floor1 = np.array([u[11], u[12], u[13]], dtype=float)
            floor2 = np.array([u[21], u[22], u[23]], dtype=float)
            floor3 = np.array([u[31], u[32], u[33]], dtype=float)
            iso_disp = np.array([u[1] - u[101], u[2] - u[102], u[3] - u[103]], dtype=float)
            base_mean = float(np.mean(base))
            f1 = float(np.mean(floor1))
            f2 = float(np.mean(floor2))
            roof = float(np.mean(floor3))
            rows.append(
                {
                    "time_s": float(frame.frameValue),
                    "roof_displacement_m": roof,
                    "story_1_drift_ratio": (f1 - base_mean) / PARAMS["story_height_m"],
                    "story_2_drift_ratio": (f2 - f1) / PARAMS["story_height_m"],
                    "story_3_drift_ratio": (roof - f2) / PARAMS["story_height_m"],
                    "isolation_displacement_m": float(np.mean(iso_disp)),
                    "iso1_disp_m": float(iso_disp[0]),
                    "iso2_disp_m": float(iso_disp[1]),
                    "iso3_disp_m": float(iso_disp[2]),
                }
            )
    finally:
        odb.close()

    iso_disp = np.column_stack(
        [
            np.array([row["iso1_disp_m"] for row in rows], dtype=float),
            np.array([row["iso2_disp_m"] for row in rows], dtype=float),
            np.array([row["iso3_disp_m"] for row in rows], dtype=float),
        ]
    )
    iso_force, iso_z = boucwen_histories(iso_disp)
    for i, row in enumerate(rows):
        row["iso1_force_N"] = float(iso_force[i, 0])
        row["iso2_force_N"] = float(iso_force[i, 1])
        row["iso3_force_N"] = float(iso_force[i, 2])
        row["iso1_z_m"] = float(iso_z[i, 0])
        row["iso2_z_m"] = float(iso_z[i, 1])
        row["iso3_z_m"] = float(iso_z[i, 2])
        row["total_base_shear_N"] = float(np.sum(iso_force[i, :]))

    fieldnames = [
        "time_s",
        "roof_displacement_m",
        "story_1_drift_ratio",
        "story_2_drift_ratio",
        "story_3_drift_ratio",
        "isolation_displacement_m",
        "total_base_shear_N",
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
    with HISTORY_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with HYSTERESIS_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(fieldnames[:1] + fieldnames[7:])
        for row in rows:
            writer.writerow([row[name] for name in fieldnames[:1] + fieldnames[7:]])


def read_csv_history() -> dict[str, np.ndarray]:
    with HISTORY_CSV.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows:
        raise ValueError(f"{HISTORY_CSV} has no data rows")
    required = [
        "time_s",
        "roof_displacement_m",
        "story_1_drift_ratio",
        "story_2_drift_ratio",
        "story_3_drift_ratio",
        "isolation_displacement_m",
        "total_base_shear_N",
        "iso1_disp_m",
        "iso1_force_N",
        "iso2_disp_m",
        "iso2_force_N",
        "iso3_disp_m",
        "iso3_force_N",
    ]
    missing = [field for field in required if field not in fields]
    if missing:
        raise ValueError(f"{HISTORY_CSV} is missing required columns: {missing}")
    return {field: np.array([float(row[field]) for row in rows], dtype=float) for field in fields}


def payload_from_csv() -> dict[str, Any]:
    ag, input_meta = read_and_normalize_motion()
    data = read_csv_history()
    time = data["time_s"]
    roof = data["roof_displacement_m"]
    story = np.column_stack([data["story_1_drift_ratio"], data["story_2_drift_ratio"], data["story_3_drift_ratio"]])
    iso = data["isolation_displacement_m"]
    base_shear = data["total_base_shear_N"]
    iso_disp = np.column_stack([data["iso1_disp_m"], data["iso2_disp_m"], data["iso3_disp_m"]])
    iso_force = np.column_stack([data["iso1_force_N"], data["iso2_force_N"], data["iso3_force_N"]])
    if all(name in data for name in ["iso1_z_m", "iso2_z_m", "iso3_z_m"]):
        iso_z = np.column_stack([data["iso1_z_m"], data["iso2_z_m"], data["iso3_z_m"]])
    else:
        _, iso_z = boucwen_histories(iso_disp)

    if time.size <= ag.size and np.allclose(time, np.arange(time.size) * PARAMS["dt"], rtol=0.0, atol=1e-7):
        accel = ag[: time.size]
    else:
        accel = np.interp(time, np.arange(ag.size) * PARAMS["dt"], ag)

    fundamental, modal_periods, alpha_m, beta_k = reference_modal_properties()
    peaks = peaks_from_histories(roof, story, iso, base_shear, iso_force, iso_disp)
    metadata = {}
    if METADATA_PATH.exists():
        metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))

    return {
        "schema_version": "example6-response-v1",
        "software": "ABAQUS",
        "status": "completed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "units": {"length": "m", "force": "N", "mass": "kg", "time": "s"},
        "input": input_meta,
        "model": PARAMS,
        "damping": {
            "type": "Abaqus-equivalent Rayleigh damping",
            "target_modes": [1, 2],
            "target_damping_ratio": PARAMS["damping_ratio"],
            "alpha_mass": alpha_m,
            "beta_stiffness_initial": beta_k,
            "implementation": [
                "alphaM*M represented by horizontal DASHPOT1 elements at floor mass nodes",
                "betaK*K represented by Abaqus material damping for B23 frame members",
                "betaK*Kinitial represented inside the Bouc-Wen UEL by betaK*k0 relative-velocity force",
            ],
        },
        "fundamental_period_s": fundamental,
        "modal_periods_s": modal_periods[:6],
        "time_s": time.tolist(),
        "ground_accel_m_per_s2": accel.tolist(),
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
            "analysis": "Abaqus/Standard direct transient dynamic analysis with compiled Bouc-Wen UEL",
            "attempted_steps": max(0, int(time.size) - 1),
            "completed_steps": max(0, int(time.size) - 1),
            "failed_step_indices": [],
            "failed_step_count": 0,
            "failed_step_rate": 0.0,
        },
        "abaqus_boucwen_implementation": {
            "direct_boucwen_state_evolution": True,
            "substituted_bilinear_model": False,
            "uel_source": "verification/abaqus/boucwen_uel.for",
            "state_update": "z_{n+1}-z_n-du*(Ao-abs(z_{n+1})^n*(gamma+beta*sign(du*z_{n+1})))=0",
            "force_convention": "F = alpha*k0*u + (1-alpha)*k0*z; post-processed base shear reports restoring force without Rayleigh damping force",
        },
        "assumptions": [
            "The Abaqus model uses B23 elastic beam elements, back-calculated rectangular sections matching the specified A and out-of-plane I, and horizontal-only anisotropic MASS elements.",
            "A tiny nonzero B23 material density is used only to satisfy Abaqus/Standard general dynamic input checks; frame self-mass is intentionally negligible to match the MATLAB/OpenSeesPy lumped horizontal mass model.",
            "Floor horizontal diaphragm constraints are implemented with *EQUATION equations tying U1 at each floor.",
            "Ground motion is applied as equivalent horizontal inertia force -m*ag(t), matching the MATLAB relative-coordinate formulation.",
            "Abaqus UEL provides the direct Bouc-Wen force during the solve; ODB post-processing recomputes the same restoring force from final displacement histories for auditable output.",
            "Gravity and P-Delta effects are intentionally excluded.",
        ],
        "output_files": {
            "json": "verification/abaqus/response_output.json",
            "time_histories_csv": "verification/abaqus/abaqus_history.csv",
            "isolator_hysteresis_csv": "verification/abaqus/isolator_hysteresis.csv",
            "input_deck": "verification/abaqus/example6_model.inp",
            "uel": "verification/abaqus/boucwen_uel.for",
        },
        "abaqus_model_metadata": metadata,
    }


def unavailable_payload(status: str, message: str) -> dict[str, Any]:
    try:
        _, input_meta = read_and_normalize_motion()
    except Exception as exc:
        input_meta = {"error": f"{type(exc).__name__}: {exc}"}
    return {
        "schema_version": "example6-response-v1",
        "software": "ABAQUS",
        "status": status,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "message": message,
        "input": input_meta,
        "model": PARAMS,
        "abaqus_boucwen_implementation": {
            "direct_boucwen_state_evolution": False,
            "substituted_bilinear_model": False,
        },
    }


def write_output(payload: dict[str, Any]) -> None:
    with OUTPUT_PATH.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def main() -> int:
    try:
        if ODB_PATH.exists():
            export_history_from_odb()
        if not HISTORY_CSV.exists():
            payload = unavailable_payload(
                "abaqus_results_missing",
                "example6_model.odb or abaqus_history.csv was not found; run the Abaqus job first.",
            )
            write_output(payload)
            return 2
        payload = payload_from_csv()
        write_output(payload)
        return 0
    except Exception as exc:
        write_output(unavailable_payload("failed", f"{type(exc).__name__}: {exc}"))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
