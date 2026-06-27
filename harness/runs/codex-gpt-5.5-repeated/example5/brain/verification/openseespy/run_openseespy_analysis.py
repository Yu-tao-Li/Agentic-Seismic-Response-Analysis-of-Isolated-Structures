"""OpenSeesPy nonlinear time-history analysis for the isolated 4-DOF building."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import openseespy.opensees as ops


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
GM_FILE = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"


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
    if omega.size < 4:
        raise RuntimeError("Expected four positive elastic frequencies.")
    mat = np.array([[1.0 / (2.0 * omega[0]), omega[0] / 2.0],
                    [1.0 / (2.0 * omega[3]), omega[3] / 2.0]])
    alpha_m, beta_k = np.linalg.solve(mat, np.array([damping_ratio, damping_ratio]))
    return float(alpha_m), float(beta_k), omega


def material_response(ele_tag: int, lower_node: int, upper_node: int, stiffness: float) -> tuple[float, float, float]:
    strain = ops.nodeDisp(upper_node, 1) - (ops.nodeDisp(lower_node, 1) if lower_node != 0 else 0.0)
    stress_resp = ops.eleResponse(ele_tag, "material", 1, "stress")
    if stress_resp:
        force = float(stress_resp[0])
    else:
        forces = ops.eleForce(ele_tag)
        force = float(forces[-1])
    plastic = strain - force / stiffness
    return float(strain), force, float(plastic)


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

    alpha_m, beta_k, omega = modal_rayleigh(masses, stiffness, damping_ratio)

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    for node in range(5):
        ops.node(node, 0.0)
    ops.fix(0, 1)
    for i, mass in enumerate(masses, start=1):
        ops.mass(i, float(mass))
    for i in range(4):
        ops.uniaxialMaterial("Steel01", i + 1, float(fy[i]), float(stiffness[i]), float(b[i]))
        ops.element("zeroLength", i + 1, i, i + 1, "-mat", i + 1, "-dir", 1, "-doRayleigh", 1)

    ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ag.tolist())
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("NewtonLineSearch", "-type", "Bisection")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    disp = np.zeros((n, 4), dtype=float)
    rel_acc = np.zeros((n, 4), dtype=float)
    story_def = np.zeros((n, 4), dtype=float)
    story_force = np.zeros((n, 4), dtype=float)
    story_plastic = np.zeros((n, 4), dtype=float)
    failed_steps = 0
    ok_flags = np.ones(n, dtype=bool)

    # Store the initial kinematic state; OpenSees reports zero acceleration before the first step.
    for j in range(4):
        disp[0, j] = ops.nodeDisp(j + 1, 1)
        rel_acc[0, j] = ops.nodeAccel(j + 1, 1)
    for ele in range(1, 5):
        d, f, p = material_response(ele, ele - 1, ele, stiffness[ele - 1])
        story_def[0, ele - 1] = d
        story_force[0, ele - 1] = f
        story_plastic[0, ele - 1] = p

    for istep in range(1, n):
        ok = ops.analyze(1, dt)
        if ok != 0:
            # Retry with a more conservative algorithm before recording a true failed increment.
            ops.test("NormDispIncr", 1.0e-9, 100, 0)
            ops.algorithm("Newton")
            ok = ops.analyze(1, dt)
            ops.test("NormDispIncr", 1.0e-10, 60, 0)
            ops.algorithm("NewtonLineSearch", "-type", "Bisection")
        if ok != 0:
            failed_steps += 1
            ok_flags[istep] = False
            raise RuntimeError(f"OpenSeesPy failed to converge at step {istep}, time {time[istep]:.6f}, code {ok}")

        for j in range(4):
            disp[istep, j] = ops.nodeDisp(j + 1, 1)
            rel_acc[istep, j] = ops.nodeAccel(j + 1, 1)
        for ele in range(1, 5):
            d, f, p = material_response(ele, ele - 1, ele, stiffness[ele - 1])
            story_def[istep, ele - 1] = d
            story_force[istep, ele - 1] = f
            story_plastic[istep, ele - 1] = p

    top_disp = disp[:, 3]
    top_abs_acc = rel_acc[:, 3] + ag
    iso_disp = disp[:, 0]
    base_shear = story_force[:, 0]
    yielding = np.abs(story_plastic) > 1.0e-10
    loop_area = float(np.trapezoid(base_shear, iso_disp))

    np.savetxt(SCRIPT_DIR / "topStoDisIso2.txt", np.column_stack([time, top_disp]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "topStoAccIso2.txt", np.column_stack([time, top_abs_acc]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "isolation_hysteresis.txt", np.column_stack([time, iso_disp, base_shear]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "story_hysteresis_all.txt", np.column_stack([time, story_def, story_force]), fmt="%.12e", delimiter="\t")

    response = {
        "software": "OpenSeesPy",
        "generated_by": Path(__file__).name,
        "units": {"length": "m", "mass": "kg", "time": "s", "force": "N", "acceleration": "m/s^2"},
        "model": {
            "masses_kg": masses.tolist(),
            "stiffness_N_per_m": stiffness.tolist(),
            "yield_N": fy.tolist(),
            "post_yield_ratio": b.tolist(),
            "damping_ratio": damping_ratio,
            "rayleigh_modes": [1, 4],
            "rayleigh_alphaM": alpha_m,
            "rayleigh_betaKinit": beta_k,
            "elastic_omega_rad_s": omega.tolist(),
        },
        "input": {
            "ground_motion_file": str(GM_FILE),
            "dt": dt,
            "sample_count": int(n),
            "raw_peak": raw_peak,
            "target_pga_mps2": target_pga,
            "normalized_peak_mps2": float(np.max(np.abs(ag))),
        },
        "time": time.tolist(),
        "ground_acceleration": ag.tolist(),
        "top_story_displacement": top_disp.tolist(),
        "top_story_acceleration": top_abs_acc.tolist(),
        "isolation_displacement": iso_disp.tolist(),
        "base_shear": base_shear.tolist(),
        "isolation_hysteresis": {
            "displacement": iso_disp.tolist(),
            "force": base_shear.tolist(),
            "loop_area": loop_area,
        },
        "story_displacement": story_def.tolist(),
        "story_force": story_force.tolist(),
        "story_plastic_deformation": story_plastic.tolist(),
        "yielding": yielding.tolist(),
        "convergence": {
            "failed_step_count": int(failed_steps),
            "converged": ok_flags.tolist(),
            "test": "NormDispIncr 1e-10 60 with NewtonLineSearch; one Newton retry before failure",
        },
        "peaks": {
            "top_displacement_abs_max": float(np.max(np.abs(top_disp))),
            "top_acceleration_abs_max": float(np.max(np.abs(top_abs_acc))),
            "isolation_displacement_abs_max": float(np.max(np.abs(iso_disp))),
            "base_shear_abs_max": float(np.max(np.abs(base_shear))),
        },
    }
    (SCRIPT_DIR / "response.json").write_text(json.dumps(response, indent=2), encoding="utf-8")
    print(f"OpenSeesPy nonlinear analysis complete. Samples: {n}, failed steps: {failed_steps}")
    print(f"Rayleigh alphaM = {alpha_m:.12g}, betaKinit = {beta_k:.12g}")
    print(
        "Peaks: top displacement %.12g m, top acceleration %.12g m/s^2, "
        "isolation displacement %.12g m, base shear %.12g N"
        % (
            response["peaks"]["top_displacement_abs_max"],
            response["peaks"]["top_acceleration_abs_max"],
            response["peaks"]["isolation_displacement_abs_max"],
            response["peaks"]["base_shear_abs_max"],
        )
    )


if __name__ == "__main__":
    main()
