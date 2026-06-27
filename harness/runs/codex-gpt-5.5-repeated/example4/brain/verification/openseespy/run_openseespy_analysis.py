import json
import math
import sys
import time as walltime
from pathlib import Path

import numpy as np
import openseespy.opensees as ops


def model_parameters():
    masses = np.array([2.50e5, 2.70e5, 2.70e5, 1.80e5], dtype=float)
    stiffness = np.array([50e6, 245e6, 195e6, 98e6], dtype=float)
    yield_forces = np.array([math.inf, 1225e3, 975e3, 490e3], dtype=float)
    b = np.array([math.nan, 0.0, 0.0, 0.0], dtype=float)
    bmat = np.array(
        [[1.0, 0.0, 0.0, 0.0],
         [-1.0, 1.0, 0.0, 0.0],
         [0.0, -1.0, 1.0, 0.0],
         [0.0, 0.0, -1.0, 1.0]],
        dtype=float,
    )
    k0 = bmat.T @ np.diag(stiffness) @ bmat
    return {
        "masses": masses,
        "stiffness": stiffness,
        "yield_forces": yield_forces,
        "post_yield": b,
        "c_iso": 1000e3,
        "B": bmat,
        "K0": k0,
        "M": np.diag(masses),
    }


def rayleigh_from_initial_modes(mass, stiffness, zeta=0.05):
    vals, _ = np.linalg.eig(np.linalg.solve(mass, stiffness))
    omegas = np.sqrt(np.sort(np.real(vals)))
    w1 = omegas[0]
    w2 = omegas[-1]
    mat = np.array([[1.0 / (2.0 * w1), w1 / 2.0], [1.0 / (2.0 * w2), w2 / 2.0]])
    alpha_m, beta_k = np.linalg.solve(mat, np.array([zeta, zeta]))
    return float(alpha_m), float(beta_k), omegas


def write_json(path, payload):
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        return obj

    path.write_text(json.dumps(payload, indent=2, default=convert), encoding="utf-8")


def ele_scalar_force(ele_tag):
    values = ops.eleResponse(ele_tag, "force")
    if not values:
        return float("nan")
    if len(values) == 1:
        return float(values[0])
    # zeroLength returns resisting forces at both end nodes. The second end is
    # the positive story deformation side for the element orientation used here.
    return float(values[-1])


def main():
    script_dir = Path(__file__).resolve().parent
    root_dir = script_dir.parent.parent
    input_file = root_dir / "input_data" / "Northridge_01_NO_968.txt"
    internal_log = script_dir / "openseespy_internal_log.txt"

    log_lines = []
    log_lines.append(f"OpenSeesPy analysis started: {walltime.strftime('%Y-%m-%dT%H:%M:%S')}")
    try:
        raw = np.loadtxt(input_file, dtype=float)
        raw = np.asarray(raw, dtype=float).reshape(-1)
        if raw.size == 0:
            raise RuntimeError(f"Input ground motion file is empty: {input_file}")

        dt = 0.01
        g = 9.81
        target_pga = 0.40 * g
        scale = target_pga / np.max(np.abs(raw))
        ag = raw * scale
        n = raw.size
        time_vec = np.arange(n, dtype=float) * dt

        model = model_parameters()
        alpha_m, beta_k, omegas = rayleigh_from_initial_modes(model["M"], model["K0"], 0.05)
        stiffness = model["stiffness"]
        masses = model["masses"]
        yield_forces = model["yield_forces"]

        c_story = beta_k * stiffness
        c_mass = alpha_m * masses

        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        for nd in range(5):
            ops.node(nd, 0.0)
        ops.fix(0, 1)
        for i, mass in enumerate(masses, start=1):
            ops.mass(i, float(mass))

        # Restoring elements: one linear isolation element and three ideal EPP story elements.
        ops.uniaxialMaterial("Elastic", 1, float(stiffness[0]))
        ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
        spring_ele_tags = [1]
        for story in range(2, 5):
            mat_tag = story
            ele_tag = story
            dy = yield_forces[story - 1] / stiffness[story - 1]
            ops.uniaxialMaterial("ElasticPP", mat_tag, float(stiffness[story - 1]), float(dy), float(-dy), 0.0)
            ops.element("zeroLength", ele_tag, story - 1, story, "-mat", mat_tag, "-dir", 1)
            spring_ele_tags.append(ele_tag)

        # Explicit viscous network equivalent to C = alpha*M + beta*K_initial plus isolation dashpot.
        dashpot_ele_tags = []
        ele = 100
        for story in range(1, 5):
            mat = 100 + story
            coeff = c_story[story - 1]
            if story == 1:
                coeff += model["c_iso"]
            ops.uniaxialMaterial("Viscous", mat, float(coeff), 1.0)
            ops.element("zeroLength", ele, story - 1, story, "-mat", mat, "-dir", 1)
            dashpot_ele_tags.append(ele)
            ele += 1
        for i, coeff in enumerate(c_mass, start=1):
            mat = 200 + i
            ops.uniaxialMaterial("Viscous", mat, float(coeff), 1.0)
            ops.element("zeroLength", ele, 0, i, "-mat", mat, "-dir", 1)
            dashpot_ele_tags.append(ele)
            ele += 1

        ops.timeSeries("Path", 1, "-dt", dt, "-values", *[float(x) for x in ag])
        ops.pattern("Plain", 1, 1)
        for i, mass in enumerate(masses, start=1):
            ops.load(i, float(-mass))

        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-9, 30, 0)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")

        u = np.zeros((n, 4), dtype=float)
        rel_acc = np.zeros((n, 4), dtype=float)
        top_abs_acc = np.zeros(n, dtype=float)
        iso_disp = np.zeros(n, dtype=float)
        story_force = np.zeros((n, 4), dtype=float)
        base_shear = np.zeros(n, dtype=float)
        yielded_instant = np.zeros((n, 3), dtype=int)
        yielded_cumulative = np.zeros((n, 3), dtype=int)
        first_yield_time = [math.nan, math.nan, math.nan]
        iterations = np.zeros(n, dtype=int)
        converged = np.ones(n, dtype=int)
        failed_steps = 0
        fallback_steps = 0

        top_abs_acc[0] = 0.0
        story_force[0, 0] = stiffness[0] * u[0, 0]

        for step in range(1, n):
            ok = ops.analyze(1, dt)
            if ok != 0:
                fallback_steps += 1
                ops.test("NormDispIncr", 1.0e-8, 80, 0)
                ops.algorithm("ModifiedNewton")
                ok = ops.analyze(1, dt)
                ops.algorithm("Newton")
                ops.test("NormDispIncr", 1.0e-9, 30, 0)
            if ok != 0:
                # Last attempt with two half steps preserves the requested output time.
                fallback_steps += 1
                ops.test("NormDispIncr", 1.0e-8, 100, 0)
                ops.algorithm("NewtonLineSearch", "-type", "Bisection")
                ok1 = ops.analyze(1, dt / 2.0)
                ok2 = ops.analyze(1, dt / 2.0) if ok1 == 0 else -1
                ok = 0 if (ok1 == 0 and ok2 == 0) else ok
                ops.algorithm("Newton")
                ops.test("NormDispIncr", 1.0e-9, 30, 0)
            if ok != 0:
                failed_steps += 1
                converged[step] = 0
                log_lines.append(f"Failed convergence at target step {step}, time {time_vec[step]:.6f}, OpenSees code {ok}")
                break

            for nd in range(1, 5):
                u[step, nd - 1] = ops.nodeDisp(nd, 1)
                rel_acc[step, nd - 1] = ops.nodeAccel(nd, 1)
            iso_disp[step] = u[step, 0]
            top_abs_acc[step] = rel_acc[step, 3] + ag[step]
            for j, ele_tag in enumerate(spring_ele_tags):
                if j == 0:
                    story_force[step, j] = stiffness[0] * u[step, 0]
                else:
                    story_force[step, j] = ele_scalar_force(ele_tag)
            base_shear[step] = story_force[step, 0]

            for j in range(3):
                yielded = abs(story_force[step, j + 1]) >= 0.999 * yield_forces[j + 1]
                yielded_instant[step, j] = int(yielded)
                yielded_cumulative[step, j] = max(yielded_cumulative[step - 1, j], int(yielded))
                if yielded and math.isnan(first_yield_time[j]):
                    first_yield_time[j] = float(time_vec[step])

        if failed_steps:
            end = step + 1
            u[end:, :] = np.nan
            rel_acc[end:, :] = np.nan
            top_abs_acc[end:] = np.nan
            iso_disp[end:] = np.nan
            story_force[end:, :] = np.nan
            base_shear[end:] = np.nan
            converged[end:] = 0
        else:
            step = n - 1

        np.savetxt(script_dir / "topStoDisIso1.txt", np.column_stack([time_vec, u[:, 3]]), fmt="%.10e", delimiter="\t")
        np.savetxt(script_dir / "topStoAccIso1.txt", np.column_stack([time_vec, top_abs_acc]), fmt="%.10e", delimiter="\t")
        np.savetxt(script_dir / "isoDisIso1.txt", np.column_stack([time_vec, iso_disp]), fmt="%.10e", delimiter="\t")
        np.savetxt(script_dir / "baseShearIso1.txt", np.column_stack([time_vec, base_shear]), fmt="%.10e", delimiter="\t")

        def peak(values):
            finite = np.isfinite(values)
            if not np.any(finite):
                return float("nan"), float("nan")
            idx_local = np.nanargmax(np.abs(values))
            return float(abs(values[idx_local])), float(time_vec[idx_local])

        ptd, ptd_t = peak(u[:, 3])
        pta, pta_t = peak(top_abs_acc)
        pid, pid_t = peak(iso_disp)
        pbs, pbs_t = peak(base_shear)

        response = {
            "software": "OpenSeesPy",
            "generated_at": walltime.strftime("%Y-%m-%dT%H:%M:%S"),
            "input": {
                "file": str(input_file),
                "raw_count": int(n),
                "dt": dt,
                "scale_to_mps2": float(scale),
                "normalized_pga_mps2": float(np.max(np.abs(ag))),
                "target_pga_g": 0.40,
            },
            "model": {
                "masses_kg": masses,
                "story_stiffness_N_per_m": stiffness,
                "yield_forces_N": [None, 1225e3, 975e3, 490e3],
                "post_yield_ratios": [None, 0.0, 0.0, 0.0],
                "isolation_dashpot_Ns_per_m": model["c_iso"],
                "rayleigh": {
                    "zeta": 0.05,
                    "mode_pair": [1, 4],
                    "alpha_mass": alpha_m,
                    "beta_initial_stiffness": beta_k,
                    "c_story_dashpots_Ns_per_m": c_story,
                    "c_mass_dashpots_Ns_per_m": c_mass,
                    "initial_omega_rad_per_s": omegas,
                },
            },
            "method": {
                "integrator": "OpenSees Newmark average acceleration",
                "nonlinear_update": "ElasticPP uniaxial materials in zeroLength story springs",
                "damping": "Explicit Viscous zeroLength elements matching alpha*M + beta*K_initial and isolation dashpot",
            },
            "time": time_vec,
            "top_story_displacement": u[:, 3],
            "top_story_acceleration": top_abs_acc,
            "top_story_relative_acceleration": rel_acc[:, 3],
            "isolation_displacement": iso_disp,
            "restoring_force": story_force,
            "base_shear": base_shear,
            "yielding": {
                "instant": yielded_instant,
                "cumulative": yielded_cumulative,
                "first_yield_time": first_yield_time,
            },
            "convergence": {
                "failed_step_count": int(failed_steps),
                "converged": converged,
                "iterations": iterations,
                "fallback_step_count": int(fallback_steps),
            },
            "summary": {
                "peak_top_displacement_m": ptd,
                "peak_top_displacement_time_s": ptd_t,
                "peak_top_acceleration_mps2": pta,
                "peak_top_acceleration_g": pta / 9.81 if math.isfinite(pta) else math.nan,
                "peak_top_acceleration_time_s": pta_t,
                "peak_isolation_displacement_m": pid,
                "peak_isolation_displacement_time_s": pid_t,
                "peak_base_shear_N": pbs,
                "peak_base_shear_time_s": pbs_t,
                "first_yield_time_s": first_yield_time,
                "failed_step_count": int(failed_steps),
            },
        }
        write_json(script_dir / "response_openseespy.json", response)
        log_lines.append(f"OpenSeesPy analysis completed. Failed steps: {failed_steps}; fallback steps: {fallback_steps}")
    except Exception as exc:
        failure = {"software": "OpenSeesPy", "status": "failed", "message": str(exc)}
        write_json(script_dir / "response_openseespy.json", failure)
        log_lines.append(f"OpenSeesPy analysis failed: {exc}")
        internal_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        raise
    finally:
        try:
            ops.wipe()
        except Exception:
            pass
        internal_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
