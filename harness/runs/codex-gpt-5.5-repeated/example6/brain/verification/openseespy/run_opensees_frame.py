"""OpenSeesPy 2D isolated frame benchmark for Example 6."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import openseespy.opensees as ops


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
LOG_FILE = SCRIPT_DIR / "execution_log.txt"


def log(message: str) -> None:
    print(message)
    with LOG_FILE.open("a", encoding="utf-8") as f:
        f.write(message + "\n")


def params() -> dict[str, float]:
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
    }


def node_tag(level: int, col: int) -> int:
    return 100 + level * 10 + col


def ground_node_tag(col: int) -> int:
    return 10 + col


def build_model(p: dict[str, float]) -> tuple[list[int], list[int], list[int]]:
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    for col in range(3):
        x = col * p["bay_width"]
        ops.node(ground_node_tag(col), x, 0.0)
        ops.fix(ground_node_tag(col), 1, 1, 1)

    for level in range(4):
        y = level * p["story_height"]
        for col in range(3):
            x = col * p["bay_width"]
            tag = node_tag(level, col)
            ops.node(tag, x, y)
            if level >= 1:
                ops.mass(tag, p["floor_mass"] / 3.0, 0.0, 0.0)

    # Rigid-diaphragm-type horizontal compatibility only.
    floor_masters = []
    for level in range(1, 4):
        master = node_tag(level, 1)
        floor_masters.append(master)
        ops.equalDOF(master, node_tag(level, 0), 1)
        ops.equalDOF(master, node_tag(level, 2), 1)

    ops.geomTransf("Linear", 1)
    ele = 1
    for col in range(3):
        for level in range(3):
            ops.element(
                "elasticBeamColumn",
                ele,
                node_tag(level, col),
                node_tag(level + 1, col),
                p["Ac"],
                p["E"],
                p["Ic"],
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
                p["Ab"],
                p["E"],
                p["Ib"],
                1,
            )
            ele += 1

    isolator_elements = []
    for col in range(3):
        bw_tag = 1000 + col
        vert_tag = 2000 + col
        ele_tag = 3000 + col
        ops.uniaxialMaterial(
            "BoucWen",
            bw_tag,
            p["alpha"],
            p["k0"],
            p["bouc_n"],
            p["gamma_bw"],
            p["beta_bw"],
            p["Ao"],
            p["deltaA"],
            p["deltaNu"],
            p["deltaEta"],
        )
        ops.uniaxialMaterial("Elastic", vert_tag, p["kv"])
        ops.element(
            "zeroLength",
            ele_tag,
            ground_node_tag(col),
            node_tag(0, col),
            "-mat",
            bw_tag,
            vert_tag,
            "-dir",
            1,
            2,
            "-doRayleigh",
            1,
        )
        isolator_elements.append(ele_tag)

    base_nodes = [node_tag(0, col) for col in range(3)]
    return floor_masters, base_nodes, isolator_elements


def get_periods(count: int = 2) -> tuple[np.ndarray, np.ndarray]:
    eigvals = np.array(ops.eigen("-fullGenLapack", count), dtype=float)
    eigvals = eigvals[eigvals > 0.0]
    omegas = np.sqrt(eigvals)
    periods = 2.0 * math.pi / omegas
    return periods, omegas


def rayleigh_coefficients(w1: float, w2: float, xi: float) -> tuple[float, float]:
    a = np.array([[1.0 / (2.0 * w1), w1 / 2.0], [1.0 / (2.0 * w2), w2 / 2.0]])
    b = np.array([xi, xi])
    alpha_m, beta_k = np.linalg.solve(a, b)
    return float(alpha_m), float(beta_k)


def material_stress_strain(ele_tag: int) -> tuple[float, float]:
    response = ops.eleResponse(ele_tag, "material", 1, "stressStrain")
    if response is None or len(response) < 2:
        return float("nan"), float("nan")
    return float(response[0]), float(response[1])


def collect_response(
    p: dict[str, float],
    floor_masters: list[int],
    base_nodes: list[int],
    isolator_elements: list[int],
) -> dict[str, Any]:
    floor = np.array([ops.nodeDisp(tag, 1) for tag in floor_masters], dtype=float)
    bases = np.array([ops.nodeDisp(tag, 1) for tag in base_nodes], dtype=float)
    forces = []
    defos = []
    for ele_tag in isolator_elements:
        force, deformation = material_stress_strain(ele_tag)
        forces.append(force)
        defos.append(deformation)
    forces_arr = np.array(forces, dtype=float)
    defos_arr = np.array(defos, dtype=float)
    base_avg = float(np.mean(bases))
    drifts = np.array(
        [
            (floor[0] - base_avg) / p["story_height"],
            (floor[1] - floor[0]) / p["story_height"],
            (floor[2] - floor[1]) / p["story_height"],
        ],
        dtype=float,
    )
    return {
        "roof": float(floor[2]),
        "floor": floor,
        "bases": bases,
        "drifts": drifts,
        "iso_displacements": defos_arr,
        "iso_forces": forces_arr,
        "total_base_shear": float(np.sum(forces_arr)),
    }


def assumptions() -> dict[str, str]:
    return {
        "frame": "OpenSeesPy 2D elasticBeamColumn elements, 3 DOF per frame node.",
        "mass": "Only specified lumped floor mass is assigned to horizontal floor DOFs; steel density is recorded but not added as self-mass.",
        "constraints": "At each floor, equalDOF constrains horizontal displacement of side nodes to the center node; vertical translations and rotations remain independent.",
        "isolation": "Three zeroLength elements connect fixed ground nodes to top isolator/base column nodes with BoucWen in local x and Elastic in local y.",
        "story_drift": "Story 1 drift uses floor-1 diaphragm displacement minus average top-of-isolator displacement divided by story height.",
        "base_shear": "Total base shear is the sum of horizontal BoucWen material resisting forces from the three isolator zeroLength elements.",
        "coordinates": "UniformExcitation response is relative to the moving support in the frame x direction.",
        "bouc_wen_state": "OpenSeesPy does not expose the internal z history through eleResponse; force-deformation histories are recorded.",
    }


def write_csv(path: Path, header: list[str], rows: list[list[float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    if LOG_FILE.exists():
        LOG_FILE.unlink()
    log("OpenSeesPy 2D frame Bouc-Wen benchmark started")

    p = params()
    input_file = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"
    raw_ag = np.loadtxt(input_file, dtype=float).reshape(-1)
    if raw_ag.size == 0:
        raise RuntimeError(f"Ground-motion file is empty: {input_file}")
    raw_pga = float(np.max(np.abs(raw_ag)))
    scale = p["target_pga"] / raw_pga
    ag = raw_ag * scale
    normalized_file = SCRIPT_DIR / "normalized_ground_accel_mps2.txt"
    np.savetxt(normalized_file, ag, fmt="%.16e")
    log(
        f"Input samples: {ag.size}, raw PGA: {raw_pga:.12g}, "
        f"scale: {scale:.12g}, normalized PGA: {np.max(np.abs(ag)):.12g} m/s^2"
    )

    floor_masters, base_nodes, isolator_elements = build_model(p)
    periods, omegas = get_periods(2)
    alpha_m, beta_k = rayleigh_coefficients(float(omegas[0]), float(omegas[1]), p["damping_ratio"])
    ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)
    log(f"Elastic modal periods used for damping: {periods[0]:.12g} s, {periods[1]:.12g} s")
    log(f"Rayleigh coefficients: alphaM = {alpha_m:.12g}, betaKinit = {beta_k:.12g}")

    # Passing values avoids OpenSees' C++ file-reader limitation with Unicode paths.
    ops.timeSeries("Path", 1, "-dt", p["dt"], "-values", *ag.tolist(), "-factor", 1.0)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-8, 35, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    n_steps = ag.size
    n_intervals = n_steps - 1
    times = np.arange(n_steps, dtype=float) * p["dt"]
    roof = np.zeros(n_steps)
    floors = np.zeros((n_steps, 3))
    bases = np.zeros((n_steps, 3))
    drifts = np.zeros((n_steps, 3))
    iso_disp = np.zeros((n_steps, 3))
    iso_force = np.zeros((n_steps, 3))
    base_shear = np.zeros(n_steps)
    iter_hist = np.zeros(n_intervals, dtype=int)
    failed_steps = 0
    fallback_steps = 0

    response = collect_response(p, floor_masters, base_nodes, isolator_elements)
    roof[0] = response["roof"]
    floors[0, :] = response["floor"]
    bases[0, :] = response["bases"]
    drifts[0, :] = response["drifts"]
    iso_disp[0, :] = response["iso_displacements"]
    iso_force[0, :] = response["iso_forces"]
    base_shear[0] = response["total_base_shear"]

    for step in range(n_intervals):
        ok = ops.analyze(1, p["dt"])
        iters = int(ops.testIter()) if hasattr(ops, "testIter") else -1
        if ok != 0:
            fallback_steps += 1
            log(f"Step {step + 1}: Newton failed with code {ok}; trying NewtonLineSearch")
            ops.algorithm("NewtonLineSearch", "-type", "Bisection")
            ok = ops.analyze(1, p["dt"])
            iters = int(ops.testIter()) if hasattr(ops, "testIter") else iters
            ops.algorithm("Newton")
        if ok != 0:
            failed_steps += 1
            log(f"Step {step + 1}: failed after fallback with code {ok}")
            break
        iter_hist[step] = iters
        response = collect_response(p, floor_masters, base_nodes, isolator_elements)
        idx = step + 1
        roof[idx] = response["roof"]
        floors[idx, :] = response["floor"]
        bases[idx, :] = response["bases"]
        drifts[idx, :] = response["drifts"]
        iso_disp[idx, :] = response["iso_displacements"]
        iso_force[idx, :] = response["iso_forces"]
        base_shear[idx] = response["total_base_shear"]

    if failed_steps:
        valid_len = int(np.count_nonzero(iter_hist) + 1)
        roof[valid_len:] = np.nan
        floors[valid_len:, :] = np.nan
        bases[valid_len:, :] = np.nan
        drifts[valid_len:, :] = np.nan
        iso_disp[valid_len:, :] = np.nan
        iso_force[valid_len:, :] = np.nan
        base_shear[valid_len:] = np.nan

    loop_area = np.array(
        [
            np.trapezoid(iso_force[:, j], iso_disp[:, j])
            for j in range(3)
        ],
        dtype=float,
    )

    peaks = {
        "roof_displacement_abs_max": float(np.nanmax(np.abs(roof))),
        "max_interstory_drift_ratio_abs": float(np.nanmax(np.abs(drifts))),
        "isolation_displacement_abs_max": float(np.nanmax(np.abs(iso_disp))),
        "mean_isolation_displacement_abs_max": float(np.nanmax(np.abs(np.nanmean(iso_disp, axis=1)))),
        "total_base_shear_abs_max": float(np.nanmax(np.abs(base_shear))),
        "isolator_loop_area": loop_area.tolist(),
    }

    result: dict[str, Any] = {
        "software": "OpenSeesPy",
        "status": "OK" if failed_steps == 0 else "FAILED",
        "generated_in_this_run": True,
        "time_step": p["dt"],
        "num_input_samples": int(n_steps),
        "num_analysis_intervals": int(n_intervals),
        "input_file": os.path.relpath(input_file, ROOT_DIR),
        "raw_pga": raw_pga,
        "normalization_scale": float(scale),
        "normalized_pga": float(np.max(np.abs(ag))),
        "fundamental_period": float(periods[0]),
        "elastic_periods": periods.tolist(),
        "damping": {
            "ratio_modes_1_2": p["damping_ratio"],
            "rayleigh_alpha_m": alpha_m,
            "rayleigh_beta_k_initial": beta_k,
            "stiffness_matrix": "OpenSees initial stiffness Rayleigh term; zeroLength isolators created with -doRayleigh 1.",
        },
        "bouc_wen_update": "OpenSees uniaxialMaterial BoucWen(alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu, deltaEta); internal z is not exposed by eleResponse.",
        "modeling_assumptions": assumptions(),
        "time": times.tolist(),
        "ground_acceleration": ag.tolist(),
        "roof_displacement": roof.tolist(),
        "floor_displacements": floors.tolist(),
        "base_node_displacements": bases.tolist(),
        "interstory_drift_ratios": drifts.tolist(),
        "isolation_displacement": np.nanmean(iso_disp, axis=1).tolist(),
        "isolator_displacements": iso_disp.tolist(),
        "total_base_shear": base_shear.tolist(),
        "isolator_forces": iso_force.tolist(),
        "bouc_wen_z": None,
        "hysteresis_loop_area": loop_area.tolist(),
        "peak_responses": peaks,
        "convergence": {
            "test_type": "NormUnbalance",
            "tolerance": 1.0e-8,
            "max_iterations_per_step": 35,
            "failed_step_count": int(failed_steps),
            "fallback_step_count": int(fallback_steps),
            "max_iterations_observed": int(np.max(iter_hist)) if iter_hist.size else 0,
            "iteration_history": iter_hist.tolist(),
        },
    }

    with (SCRIPT_DIR / "opensees_results.json").open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, allow_nan=True)

    time_rows = [
        [
            times[i],
            ag[i],
            roof[i],
            drifts[i, 0],
            drifts[i, 1],
            drifts[i, 2],
            float(np.nanmean(iso_disp[i, :])),
            base_shear[i],
        ]
        for i in range(n_steps)
    ]
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
        time_rows,
    )
    hyst_rows = [
        [
            times[i],
            iso_disp[i, 0],
            iso_force[i, 0],
            iso_disp[i, 1],
            iso_force[i, 1],
            iso_disp[i, 2],
            iso_force[i, 2],
        ]
        for i in range(n_steps)
    ]
    write_csv(
        SCRIPT_DIR / "hysteresis.csv",
        [
            "time_s",
            "u_iso_1_m",
            "f_iso_1_N",
            "u_iso_2_m",
            "f_iso_2_N",
            "u_iso_3_m",
            "f_iso_3_N",
        ],
        hyst_rows,
    )
    log(
        f"OpenSeesPy benchmark completed: failed_steps={failed_steps}, "
        f"fallback_steps={fallback_steps}, max iterations={int(np.max(iter_hist)) if iter_hist.size else 0}"
    )


if __name__ == "__main__":
    main()
