import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_PATH = ROOT / "input_data" / "GM1.txt"


def build_matrices():
    masses = np.array([500.0, 1000.0, 1000.0, 1000.0])
    stiffnesses = np.array([200e3, 500e3, 500e3, 500e3])
    k1, k2, k3, k4 = stiffnesses
    mass_matrix = np.diag(masses)
    stiffness_matrix = np.array(
        [
            [k1 + k2, -k2, 0.0, 0.0],
            [-k2, k2 + k3, -k3, 0.0],
            [0.0, -k3, k3 + k4, -k4],
            [0.0, 0.0, -k4, k4],
        ],
        dtype=float,
    )
    return masses, stiffnesses, mass_matrix, stiffness_matrix


def modal_rayleigh(mass_matrix, stiffness_matrix, damping_ratio):
    eigvals = np.linalg.eigvals(np.linalg.solve(mass_matrix, stiffness_matrix))
    omegas = np.sqrt(np.sort(np.real(eigvals)))
    system = np.array(
        [[1.0 / (2.0 * omegas[0]), omegas[0] / 2.0], [1.0 / (2.0 * omegas[1]), omegas[1] / 2.0]]
    )
    alpha_m, beta_k = np.linalg.solve(system, np.array([damping_ratio, damping_ratio]))
    return omegas, float(alpha_m), float(beta_k)


def main():
    import openseespy.opensees as ops

    dt = 0.02
    g = 9.80665
    target_pga = 0.40 * g
    damping_ratio = 0.05
    c_iso = 2e3

    raw_gm = np.loadtxt(INPUT_PATH, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw_gm)))
    if raw_peak <= 0.0:
        raise ValueError("GM1.txt has zero peak acceleration.")
    ag = raw_gm / raw_peak * target_pga
    n = int(ag.size)
    time = np.arange(n, dtype=float) * dt

    masses, stiffnesses, mass_matrix, stiffness_matrix = build_matrices()
    omegas, alpha_m, beta_k = modal_rayleigh(mass_matrix, stiffness_matrix, damping_ratio)
    stiffness_dashpots = beta_k * stiffnesses
    mass_dashpots = alpha_m * masses

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    for node in range(5):
        ops.node(node, 0.0)
    ops.fix(0, 1)
    for dof_node, mass in zip(range(1, 5), masses):
        ops.mass(dof_node, float(mass))

    mat_id = 1
    ele_id = 1
    for i, k in enumerate(stiffnesses, start=1):
        lower = i - 1
        upper = i
        ops.uniaxialMaterial("Elastic", mat_id, float(k))
        ops.element("zeroLength", ele_id, lower, upper, "-mat", mat_id, "-dir", 1)
        mat_id += 1
        ele_id += 1

    for i, c in enumerate(stiffness_dashpots, start=1):
        lower = i - 1
        upper = i
        if i == 1:
            c += c_iso
        ops.uniaxialMaterial("Viscous", mat_id, float(c), 1.0)
        ops.element("zeroLength", ele_id, lower, upper, "-mat", mat_id, "-dir", 1)
        mat_id += 1
        ele_id += 1

    for i, c in enumerate(mass_dashpots, start=1):
        ops.uniaxialMaterial("Viscous", mat_id, float(c), 1.0)
        ops.element("zeroLength", ele_id, 0, i, "-mat", mat_id, "-dir", 1)
        mat_id += 1
        ele_id += 1

    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ag.tolist(), "-factor", 1.0)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    u = np.zeros((4, n), dtype=float)
    v = np.zeros((4, n), dtype=float)
    rel_acc = np.zeros((4, n), dtype=float)
    top_absolute_acc = np.zeros(n, dtype=float)
    status_codes = []

    # Initial relative acceleration is computed from the common matrix equation
    # so the first sample is directly comparable with MATLAB.
    rel_acc[:, 0] = -ag[0]
    top_absolute_acc[0] = rel_acc[3, 0] + ag[0]

    for step in range(1, n):
        status = int(ops.analyze(1, dt))
        status_codes.append(status)
        if status != 0:
            raise RuntimeError(f"OpenSees analyze failed at step {step} with code {status}")
        for local_idx, node in enumerate(range(1, 5)):
            u[local_idx, step] = ops.nodeDisp(node, 1)
            v[local_idx, step] = ops.nodeVel(node, 1)
            rel_acc[local_idx, step] = ops.nodeAccel(node, 1)
        top_absolute_acc[step] = rel_acc[3, step] + ag[step]

    top_displacement = u[3, :]
    isolation_displacement = u[0, :]

    np.savetxt(SCRIPT_DIR / "normalized_ground_motion_mps2.txt", np.column_stack([time, ag]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "topStoDis.txt", np.column_stack([time, top_displacement]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "topStoAcc.txt", np.column_stack([time, top_absolute_acc]), fmt="%.12e", delimiter="\t")
    np.savetxt(
        SCRIPT_DIR / "isolationDisplacement.txt",
        np.column_stack([time, isolation_displacement]),
        fmt="%.12e",
        delimiter="\t",
    )

    output = {
        "software": "OpenSeesPy",
        "status": "success",
        "model": {
            "units": "SI",
            "masses_kg": masses.tolist(),
            "stiffnesses_N_per_m": stiffnesses.tolist(),
            "isolation_damping_Ns_per_m": c_iso,
            "rayleigh_damping_ratio": damping_ratio,
            "rayleigh_alpha_mass": alpha_m,
            "rayleigh_beta_stiffness": beta_k,
            "modal_frequencies_rad_per_s": omegas.tolist(),
            "damping_realization": (
                "Rayleigh matrix implemented with physical viscous zeroLength elements: "
                "beta*K as story dashpots and alpha*M as ground dashpots, plus isolation dashpot."
            ),
        },
        "integration": {"method": "OpenSees Newmark average acceleration", "gamma": 0.5, "beta": 0.25},
        "input": {
            "source": str(INPUT_PATH),
            "raw_peak": raw_peak,
            "target_pga_mps2": target_pga,
            "pga_after_normalization_mps2": float(np.max(np.abs(ag))),
            "pga_after_normalization_g": float(np.max(np.abs(ag)) / g),
        },
        "time_step": dt,
        "number_of_samples": n,
        "time": time.tolist(),
        "ground_acceleration_mps2": ag.tolist(),
        "top_story_displacement_m": top_displacement.tolist(),
        "top_story_acceleration_mps2": top_absolute_acc.tolist(),
        "top_story_relative_acceleration_mps2": rel_acc[3, :].tolist(),
        "isolation_layer_displacement_m": isolation_displacement.tolist(),
        "opensees_status_codes": status_codes,
    }
    with (SCRIPT_DIR / "response.json").open("w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    print("OpenSeesPy isolated-building analysis completed.")
    print(f"Samples: {n}, dt: {dt:.6f} s, normalized PGA: {np.max(np.abs(ag)):.12g} m/s^2")
    print(f"Peak top displacement: {np.max(np.abs(top_displacement)):.12g} m")
    print(f"Peak top absolute acceleration: {np.max(np.abs(top_absolute_acc)):.12g} m/s^2")
    print(f"Peak isolation displacement: {np.max(np.abs(isolation_displacement)):.12g} m")
    ops.wipe()


if __name__ == "__main__":
    main()
