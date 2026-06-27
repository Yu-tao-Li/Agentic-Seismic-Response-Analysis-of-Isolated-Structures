import json
import math
import os
from pathlib import Path

import numpy as np
import openseespy.opensees as ops


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
GM_PATH = ROOT_DIR / "input_data" / "GM1.txt"

DT = 0.02
G = 9.81
TARGET_PGA = 0.40 * G
MASS_PER_STORY = 1000.0
STORY_STIFFNESS = 500.0e3
ZETA = 0.05


def rayleigh_coefficients(mass, stiffness, zeta):
    m = mass * np.eye(3)
    k = stiffness
    k_mat = np.array(
        [[2.0 * k, -k, 0.0], [-k, 2.0 * k, -k], [0.0, -k, k]],
        dtype=float,
    )
    eigvals = np.linalg.eigvals(np.linalg.solve(m, k_mat))
    omega = np.sort(np.sqrt(np.real(eigvals)))
    a = np.array(
        [[1.0 / (2.0 * omega[0]), omega[0] / 2.0],
         [1.0 / (2.0 * omega[1]), omega[1] / 2.0]],
        dtype=float,
    )
    alpha_m, beta_k = np.linalg.solve(a, np.array([zeta, zeta], dtype=float))
    return float(alpha_m), float(beta_k), omega, k_mat


def initial_relative_acceleration(ag0, alpha_m, beta_k, k_mat):
    m = MASS_PER_STORY * np.eye(3)
    c = alpha_m * m + beta_k * k_mat
    u0 = np.zeros(3)
    v0 = np.zeros(3)
    p0 = -m @ np.ones(3) * ag0
    return np.linalg.solve(m, p0 - c @ v0 - k_mat @ u0)


def main():
    os.chdir(SCRIPT_DIR)

    raw = np.loadtxt(GM_PATH, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw)))
    if raw_peak <= 0.0:
        raise ValueError("GM1.txt has zero peak acceleration.")

    ag = raw / raw_peak * TARGET_PGA
    n = int(ag.size)
    time = np.arange(n, dtype=float) * DT
    input_pga = float(np.max(np.abs(ag)))
    eq_force = -MASS_PER_STORY * ag

    alpha_m, beta_k, omega, k_mat = rayleigh_coefficients(
        MASS_PER_STORY, STORY_STIFFNESS, ZETA
    )

    np.savetxt(SCRIPT_DIR / "normalized_ground_acceleration.txt", np.column_stack([time, ag]), fmt="%.12e")
    np.savetxt(SCRIPT_DIR / "opensees_eq_force_N.txt", eq_force, fmt="%.12e")

    ops.wipe()
    ops.model("BasicBuilder", "-ndm", 1, "-ndf", 1)
    for node in range(4):
        ops.node(node, 0.0)
    ops.fix(0, 1)
    for node in (1, 2, 3):
        ops.mass(node, MASS_PER_STORY)

    ops.uniaxialMaterial("Elastic", 1, STORY_STIFFNESS)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1, "-doRayleigh", 1)
    ops.element("zeroLength", 2, 1, 2, "-mat", 1, "-dir", 1, "-doRayleigh", 1)
    ops.element("zeroLength", 3, 2, 3, "-mat", 1, "-dir", 1, "-doRayleigh", 1)
    ops.rayleigh(alpha_m, beta_k, 0.0, 0.0)

    ops.timeSeries("Path", 1, "-dt", DT, "-filePath", "opensees_eq_force_N.txt", "-factor", 1.0)
    ops.pattern("Plain", 1, 1)
    for node in (1, 2, 3):
        ops.load(node, 1.0)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    top_disp = np.zeros(n, dtype=float)
    top_acc_abs = np.zeros(n, dtype=float)
    rel_acc0 = initial_relative_acceleration(float(ag[0]), alpha_m, beta_k, k_mat)
    top_acc_abs[0] = rel_acc0[2] + ag[0]

    ok_steps = 0
    failures = []
    for i in range(1, n):
        ok = ops.analyze(1, DT)
        if ok != 0:
            failures.append({"step": i, "time": float(i * DT), "code": int(ok)})
            break
        ok_steps += 1
        top_disp[i] = ops.nodeDisp(3, 1)
        top_acc_abs[i] = ops.nodeAccel(3, 1) + ag[i]

    if failures:
        raise RuntimeError("OpenSees analysis failed: {}".format(failures))

    np.savetxt(SCRIPT_DIR / "topStoDis.txt", np.column_stack([time, top_disp]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "topStoAcc.txt", np.column_stack([time, top_acc_abs]), fmt="%.12e", delimiter="\t")

    response = {
        "software": "OpenSeesPy",
        "model": "three_story_linear_shear_building_relative_coordinates",
        "output_definition": {
            "top_story_displacement": "relative displacement of story 3 with respect to ground, m",
            "top_story_acceleration": "absolute acceleration of story 3, m/s^2",
        },
        "time": time.tolist(),
        "top_story_displacement": top_disp.tolist(),
        "top_story_acceleration": top_acc_abs.tolist(),
        "normalized_ground_acceleration": ag.tolist(),
        "input_pga_after_normalization": input_pga,
        "dt": DT,
        "num_samples": n,
        "mass_per_story_kg": MASS_PER_STORY,
        "story_stiffness_N_per_m": STORY_STIFFNESS,
        "damping_ratio": ZETA,
        "rayleigh_alpha_m": alpha_m,
        "rayleigh_beta_k": beta_k,
        "circular_frequencies_rad_s": omega.tolist(),
        "newmark_gamma": 0.5,
        "newmark_beta": 0.25,
        "opensees_completed_steps": ok_steps,
    }
    (SCRIPT_DIR / "response.json").write_text(json.dumps(response, indent=2), encoding="utf-8")

    print("OpenSeesPy analysis complete. Samples: {}, dt: {:.8f}, normalized PGA: {:.8f} m/s^2".format(n, DT, input_pga))
    print("Rayleigh alpha_m: {:.12g}, beta_k: {:.12g}".format(alpha_m, beta_k))
    print("Peak |top displacement|: {:.12g} m".format(float(np.max(np.abs(top_disp)))))
    print("Peak |top absolute acceleration|: {:.12g} m/s^2".format(float(np.max(np.abs(top_acc_abs)))))


if __name__ == "__main__":
    main()
