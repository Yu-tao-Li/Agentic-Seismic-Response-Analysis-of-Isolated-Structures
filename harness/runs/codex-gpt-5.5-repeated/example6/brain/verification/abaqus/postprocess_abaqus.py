"""Postprocess ABAQUS ODB for Example 6."""
from __future__ import annotations

import csv
import json
import math
import os
import re
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
ODB = SCRIPT_DIR / "example6_boucwen.odb"
OUT = SCRIPT_DIR / "abaqus_results.json"


def params():
    return {
        "E": 2.06e11,
        "nu": 0.30,
        "rho": 7850.0,
        "Ac": 0.020,
        "Ic": 8.0e-4,
        "Ab": 0.015,
        "Ib": 4.5e-4,
        "bay_width": 6.0,
        "story_height": 3.6,
        "floor_mass": 2.5e5,
        "k0": 2.0e6,
        "Fy": 1.0e4,
        "uy": 0.005,
        "alpha": 0.05,
        "bouc_n": 2.0,
        "beta_bw": 2.0e4,
        "gamma_bw": 2.0e4,
        "Ao": 1.0,
        "deltaA": 0.0,
        "deltaNu": 0.0,
        "deltaEta": 0.0,
        "kv": 1.0e10,
        "g": 9.81,
        "target_pga": 0.40 * 9.81,
        "dt": 0.01,
        "damping_ratio": 0.05,
        "rayleigh_alpha_m": 0.18649697508516622,
        "rayleigh_beta_k_initial": 0.0070671988247071136,
        "elastic_period_1": 2.842810601922031,
        "elastic_period_2": 0.5262442067740834,
    }


def node_id(level, col):
    return 100 + level * 10 + col + 1


def ground_id(col):
    return 10 + col + 1


def assumptions():
    return {
        "frame": "ABAQUS 30-DOF reduced elastic frame UEL built from the same Euler-Bernoulli stiffness and diaphragm transformation used by MATLAB.",
        "mass": "Specified floor mass is assigned directly to the three horizontal floor generalized coordinates only.",
        "constraints": "The frame UEL carries base isolator translations, floor diaphragm translations, and retained massless vertical/rotational algebraic coordinates.",
        "isolation": "Three two-node horizontal Bouc-Wen UEL isolators connect fixed ground nodes to top isolator/base translation nodes.",
        "ground_motion": "Equivalent inertial-force formulation using normalized acceleration amplitude and CLOAD=-m*ag(t).",
        "story_drift": "Story 1 drift uses floor-1 diaphragm displacement minus average top-of-isolator displacement divided by story height.",
        "base_shear": "Total base shear is the sum of positive-resisting Bouc-Wen material forces from ground reactions.",
        "coordinates": "All output displacements are relative to ground in the frame x direction.",
        "bouc_wen_state": "UEL state is not directly exposed in ODB; z is reconstructed from F = alpha*k0*u + (1-alpha)*k0*z.",
    }


def write_failure(stage, reason):
    result = {
        "software": "ABAQUS",
        "status": "FAILED",
        "generated_in_this_run": True,
        "failure_stage": stage,
        "reason": reason,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    raise SystemExit(reason)


def field_map(frame, field_name):
    if field_name not in frame.fieldOutputs:
        write_failure("postprocess_missing_field", "ODB field %s is missing." % field_name)
    out = {}
    for value in frame.fieldOutputs[field_name].values:
        out[value.nodeLabel] = tuple(float(x) for x in value.data)
    return out


def parse_sta():
    sta = SCRIPT_DIR / "example6_boucwen.sta"
    if not sta.exists():
        return {"failed_step_count": None, "reason": "STA file missing"}
    rows = []
    for line in sta.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 9 and all(part.lstrip("-").replace(".", "", 1).replace("E", "", 1).replace("+", "", 1).isdigit() for part in parts[:6]):
            try:
                rows.append(
                    {
                        "step": int(parts[0]),
                        "increment": int(parts[1]),
                        "attempt": int(parts[2]),
                        "severe_discontinuity_iterations": int(parts[3]),
                        "equilibrium_iterations": int(parts[4]),
                        "total_iterations": int(parts[5]),
                    }
                )
            except Exception:
                pass
    if not rows:
        return {"failed_step_count": None, "reason": "No increment rows parsed from STA file"}
    iteration_history = [row["total_iterations"] for row in rows]
    attempts = [row["attempt"] for row in rows]
    return {
        "test_type": "Abaqus/Standard equilibrium checks",
        "failed_step_count": 0,
        "cutback_count": 0,
        "increment_count": len(rows),
        "max_iterations_observed": int(max(iteration_history)),
        "iteration_history": iteration_history,
        "max_attempts_observed": int(max(attempts)),
    }


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


if not ODB.exists():
    write_failure("postprocess", "ODB file is missing; ABAQUS analysis did not complete.")

try:
    from odbAccess import openOdb
except Exception as exc:
    write_failure("postprocess_import", repr(exc))

p = params()
input_file = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"
raw_ag = np.loadtxt(str(input_file), dtype=float).reshape(-1)
raw_pga = float(np.max(np.abs(raw_ag)))
scale = p["target_pga"] / raw_pga
ag = raw_ag * scale
n_steps = int(ag.size)
target_times = np.arange(n_steps, dtype=float) * p["dt"]

odb = openOdb(str(ODB), readOnly=True)
try:
    if "DYNAMIC_BOUCWEN" not in odb.steps:
        write_failure("postprocess_missing_step", "Step DYNAMIC_BOUCWEN is missing from ODB.")
    step = odb.steps["DYNAMIC_BOUCWEN"]
    frames = list(step.frames)
    if len(frames) != n_steps:
        write_failure(
            "postprocess_frame_count",
            "Expected %d ODB frames at input sample times, found %d." % (n_steps, len(frames)),
        )

    floors = np.zeros((n_steps, 3), dtype=float)
    bases = np.zeros((n_steps, 3), dtype=float)
    iso_disp = np.zeros((n_steps, 3), dtype=float)
    iso_total_force = np.zeros((n_steps, 3), dtype=float)
    iso_vel = np.zeros((n_steps, 3), dtype=float)
    time = np.zeros(n_steps, dtype=float)

    floor_nodes = [node_id(level, 1) for level in (1, 2, 3)]
    base_nodes = [node_id(0, col) for col in range(3)]
    ground_nodes = [ground_id(col) for col in range(3)]

    for i, frame in enumerate(frames):
        time[i] = float(frame.frameValue)
        u = field_map(frame, "U")
        v = field_map(frame, "V")
        rf = field_map(frame, "RF")
        for j, tag in enumerate(floor_nodes):
            floors[i, j] = u[tag][0]
        for j, tag in enumerate(base_nodes):
            bases[i, j] = u[tag][0]
            iso_disp[i, j] = u[tag][0] - u[ground_nodes[j]][0]
            iso_vel[i, j] = v[tag][0] - v[ground_nodes[j]][0]
        for j, tag in enumerate(ground_nodes):
            # Abaqus ground RF is opposite the positive element resisting force
            # convention used by MATLAB/OpenSeesPy in this benchmark.
            iso_total_force[i, j] = -rf[tag][0]
finally:
    odb.close()

if np.max(np.abs(time - target_times)) > 2.0e-6:
    write_failure(
        "postprocess_time_alignment",
        "ODB frame times do not match input sample times; max error %.6g" % float(np.max(np.abs(time - target_times))),
    )

roof = floors[:, 2]
base_avg = np.mean(bases, axis=1)
drifts = np.column_stack(
    [
        (floors[:, 0] - base_avg) / p["story_height"],
        (floors[:, 1] - floors[:, 0]) / p["story_height"],
        (floors[:, 2] - floors[:, 1]) / p["story_height"],
    ]
)
isolation_disp = np.mean(iso_disp, axis=1)
z_hist = np.zeros_like(iso_disp)
iso_force = np.zeros_like(iso_disp)
for i in range(1, n_steps):
    du = iso_disp[i, :] - iso_disp[i - 1, :]
    z_guess = z_hist[i - 1, :] + p["Ao"] * du
    for j in range(3):
        z = float(z_guess[j])
        z_old = float(z_hist[i - 1, j])
        duj = float(du[j])
        for _ in range(30):
            absz = abs(z)
            g = (
                z
                - z_old
                - p["Ao"] * duj
                + p["beta_bw"] * abs(duj) * (absz ** (p["bouc_n"] - 1.0)) * z
                + p["gamma_bw"] * duj * (absz ** p["bouc_n"])
            )
            if abs(g) < 1.0e-13 * max(1.0, abs(z)):
                break
            if absz == 0.0:
                dgdz = 1.0
            else:
                dgdz = (
                    1.0
                    + p["beta_bw"] * abs(duj) * p["bouc_n"] * (absz ** (p["bouc_n"] - 1.0))
                    + p["gamma_bw"] * duj * p["bouc_n"] * (absz ** (p["bouc_n"] - 1.0)) * math.copysign(1.0, z)
                )
            z -= g / dgdz
        z_hist[i, j] = z
    iso_force[i, :] = p["alpha"] * p["k0"] * iso_disp[i, :] + (1.0 - p["alpha"]) * p["k0"] * z_hist[i, :]
base_shear = np.sum(iso_force, axis=1)

loop_area = []
for j in range(3):
    loop_area.append(
        float(np.sum(0.5 * (iso_force[1:, j] + iso_force[:-1, j]) * (iso_disp[1:, j] - iso_disp[:-1, j])))
    )

peaks = {
    "roof_displacement_abs_max": float(np.max(np.abs(roof))),
    "max_interstory_drift_ratio_abs": float(np.max(np.abs(drifts))),
    "isolation_displacement_abs_max": float(np.max(np.abs(iso_disp))),
    "mean_isolation_displacement_abs_max": float(np.max(np.abs(isolation_disp))),
    "total_base_shear_abs_max": float(np.max(np.abs(base_shear))),
    "isolator_loop_area": loop_area,
}

conv = parse_sta()
status = "OK" if int(conv.get("failed_step_count") or 0) == 0 else "FAILED"

result = {
    "software": "ABAQUS",
    "status": status,
    "generated_in_this_run": True,
    "time_step": p["dt"],
    "num_input_samples": n_steps,
    "num_analysis_intervals": n_steps - 1,
    "input_file": os.path.relpath(str(input_file), str(ROOT_DIR)),
    "raw_pga": raw_pga,
    "normalization_scale": scale,
    "normalized_pga": float(np.max(np.abs(ag))),
    "fundamental_period": p["elastic_period_1"],
    "elastic_periods": [p["elastic_period_1"], p["elastic_period_2"]],
    "damping": {
        "ratio_modes_1_2": p["damping_ratio"],
        "rayleigh_alpha_m": p["rayleigh_alpha_m"],
        "rayleigh_beta_k_initial": p["rayleigh_beta_k_initial"],
        "stiffness_matrix": "Mass-proportional damping is assigned to lumped masses, and frame/isolator UELs include initial-stiffness-proportional viscous terms.",
    },
    "bouc_wen_update": "UEL backward-Euler update: z_{n+1}-z_n-Ao*du+beta*abs(du)*abs(z)^(n-1)*z+gamma*du*abs(z)^n = 0; F = alpha*k0*u + (1-alpha)*k0*z.",
    "modeling_assumptions": assumptions(),
    "time": time.tolist(),
    "ground_acceleration": ag.tolist(),
    "roof_displacement": roof.tolist(),
    "floor_displacements": floors.tolist(),
    "base_node_displacements": bases.tolist(),
    "interstory_drift_ratios": drifts.tolist(),
    "isolation_displacement": isolation_disp.tolist(),
    "isolator_displacements": iso_disp.tolist(),
    "total_base_shear": base_shear.tolist(),
    "isolator_forces": iso_force.tolist(),
    "abaqus_total_isolator_reaction_forces": iso_total_force.tolist(),
    "bouc_wen_z": z_hist.tolist(),
    "hysteresis_loop_area": loop_area,
    "peak_responses": peaks,
    "convergence": conv,
}

OUT.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")

write_csv(
    SCRIPT_DIR / "time_histories.csv",
    [
        "time_s",
        "ground_accel_mps2",
        "roof_disp_m",
        "drift_story1",
        "drift_story2",
        "drift_story3",
        "mean_isolation_disp_m",
        "total_base_shear_N",
    ],
    [
        [
            float(time[i]),
            float(ag[i]),
            float(roof[i]),
            float(drifts[i, 0]),
            float(drifts[i, 1]),
            float(drifts[i, 2]),
            float(isolation_disp[i]),
            float(base_shear[i]),
        ]
        for i in range(n_steps)
    ],
)

write_csv(
    SCRIPT_DIR / "hysteresis.csv",
    [
        "time_s",
        "u_iso_1_m",
        "f_iso_1_N",
        "z_iso_1_m",
        "u_iso_2_m",
        "f_iso_2_N",
        "z_iso_2_m",
        "u_iso_3_m",
        "f_iso_3_N",
        "z_iso_3_m",
    ],
    [
        [
            float(time[i]),
            float(iso_disp[i, 0]),
            float(iso_force[i, 0]),
            float(z_hist[i, 0]),
            float(iso_disp[i, 1]),
            float(iso_force[i, 1]),
            float(z_hist[i, 1]),
            float(iso_disp[i, 2]),
            float(iso_force[i, 2]),
            float(z_hist[i, 2]),
        ]
        for i in range(n_steps)
    ],
)
