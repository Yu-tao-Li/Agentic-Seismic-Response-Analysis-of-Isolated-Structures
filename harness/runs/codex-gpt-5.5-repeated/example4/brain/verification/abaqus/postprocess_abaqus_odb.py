import json
import math
import os
import sys
import time as walltime

try:
    import numpy as np
except Exception:
    np = None

from odbAccess import openOdb


def array(values):
    if np is not None:
        return np.array(values, dtype=float)
    return values


def dump_json(path, payload):
    def default(obj):
        if np is not None and isinstance(obj, np.ndarray):
            return obj.tolist()
        if hasattr(obj, "item"):
            return obj.item()
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        return obj

    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=default)
        f.write("\n")


def save_txt(path, columns):
    with open(path, "w", encoding="utf-8") as f:
        for row in zip(*columns):
            f.write("\t".join(f"{float(x):.10e}" for x in row) + "\n")


def component(value, idx):
    data = value.data
    try:
        return float(data[idx])
    except TypeError:
        return float(data)


def field_by_node(field, node_labels, component_index=0):
    out = {label: float("nan") for label in node_labels}
    for value in field.values:
        label = int(value.nodeLabel)
        if label in out:
            out[label] = component(value, component_index)
    return out


def field_by_element(field, element_labels):
    out = {label: float("nan") for label in element_labels}
    for value in field.values:
        label = int(value.elementLabel)
        if label in out:
            out[label] = component(value, 0)
    return out


def finite_peak(times, values):
    pairs = [(abs(v), t) for t, v in zip(times, values) if math.isfinite(v)]
    if not pairs:
        return float("nan"), float("nan")
    return max(pairs, key=lambda x: x[0])


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(os.path.dirname(script_dir))
    meta_path = os.path.join(script_dir, "abaqus_input_metadata.json")
    with open(meta_path, "r", encoding="utf-8") as f:
        meta = json.load(f)

    raw_path = os.path.join(root_dir, "input_data", "Northridge_01_NO_968.txt")
    if np is not None:
        raw = np.loadtxt(raw_path, dtype=float).reshape(-1)
        ag = raw * float(meta["scale_to_mps2"])
    else:
        with open(raw_path, "r", encoding="utf-8") as f:
            raw = [float(line.strip()) for line in f if line.strip()]
        ag = [x * float(meta["scale_to_mps2"]) for x in raw]
    dt = float(meta["dt"])

    odb_path = os.path.join(script_dir, "isolated_shear_building.odb")
    log = []
    log.append(f"ABAQUS ODB postprocessing started: {walltime.strftime('%Y-%m-%dT%H:%M:%S')}")
    odb = openOdb(odb_path, readOnly=True)
    try:
        step = odb.steps["Dynamic"]
        node_labels = [2, 3, 4, 5]
        element_labels = [1, 2, 3, 4]
        yield_forces = [float("inf"), 1225e3, 975e3, 490e3]

        times = []
        top_disp = []
        top_acc = []
        top_rel_acc = []
        iso_disp = []
        story_force = []
        base_shear = []
        yielded_instant = []
        yielded_cumulative = []
        cum = [0, 0, 0]
        first_yield = [float("nan"), float("nan"), float("nan")]

        for frame in step.frames:
            t = float(frame.frameValue)
            u_field = frame.fieldOutputs["U"]
            a_field = frame.fieldOutputs["A"]
            u_nodes = field_by_node(u_field, node_labels, 0)
            a_nodes = field_by_node(a_field, node_labels, 0)
            idx = int(round(t / dt))
            if idx < 0:
                idx = 0
            if idx >= len(ag):
                idx = len(ag) - 1
            ag_t = float(ag[idx])

            if "S" in frame.fieldOutputs:
                stresses = field_by_element(frame.fieldOutputs["S"], element_labels)
                forces = [stresses[e] for e in element_labels]
            else:
                forces = [
                    float(meta["story_stiffness_N_per_m"][0]) * (u_nodes[2] - 0.0),
                    float("nan"),
                    float("nan"),
                    float("nan"),
                ]

            instant = []
            for j in range(3):
                yielded = int(math.isfinite(forces[j + 1]) and abs(forces[j + 1]) >= 0.999 * yield_forces[j + 1])
                instant.append(yielded)
                if yielded:
                    cum[j] = 1
                    if math.isnan(first_yield[j]):
                        first_yield[j] = t

            times.append(t)
            top_disp.append(u_nodes[5])
            top_rel_acc.append(a_nodes[5])
            top_acc.append(a_nodes[5] + ag_t)
            iso_disp.append(u_nodes[2])
            story_force.append(forces)
            base_shear.append(forces[0])
            yielded_instant.append(instant)
            yielded_cumulative.append(list(cum))

        save_txt(os.path.join(script_dir, "topStoDisIso1.txt"), [times, top_disp])
        save_txt(os.path.join(script_dir, "topStoAccIso1.txt"), [times, top_acc])
        save_txt(os.path.join(script_dir, "isoDisIso1.txt"), [times, iso_disp])
        save_txt(os.path.join(script_dir, "baseShearIso1.txt"), [times, base_shear])

        ptd, ptd_t = finite_peak(times, top_disp)
        pta, pta_t = finite_peak(times, top_acc)
        pid, pid_t = finite_peak(times, iso_disp)
        pbs, pbs_t = finite_peak(times, base_shear)

        final_time_error = abs(times[-1] - float(meta["total_time"])) if times else float("inf")
        failed_steps = 0 if final_time_error <= 0.5 * dt else 1
        response = {
            "software": "ABAQUS",
            "generated_at": walltime.strftime("%Y-%m-%dT%H:%M:%S"),
            "input": {
                "file": meta["input_file"],
                "raw_count": int(meta["raw_count"]),
                "dt": dt,
                "scale_to_mps2": float(meta["scale_to_mps2"]),
                "normalized_pga_mps2": float(meta["normalized_pga_mps2"]),
                "target_pga_g": 0.40,
            },
            "model": {
                "masses_kg": meta["masses_kg"],
                "story_stiffness_N_per_m": meta["story_stiffness_N_per_m"],
                "yield_forces_N": meta["yield_forces_N"],
                "post_yield_ratios": meta["post_yield_ratios"],
                "isolation_dashpot_Ns_per_m": meta["isolation_dashpot_Ns_per_m"],
                "rayleigh": meta["rayleigh"],
            },
            "method": {
                "integrator": "ABAQUS/Standard implicit dynamic with HHT alpha = 0.0 (Newmark average-acceleration limit)",
                "nonlinear_update": "T2D2 truss elements with elastic-perfectly-plastic material",
                "damping": "DASHPOT2 network matching alpha*M + beta*K_initial and isolation dashpot",
            },
            "time": times,
            "top_story_displacement": top_disp,
            "top_story_acceleration": top_acc,
            "top_story_relative_acceleration": top_rel_acc,
            "isolation_displacement": iso_disp,
            "restoring_force": story_force,
            "base_shear": base_shear,
            "yielding": {
                "instant": yielded_instant,
                "cumulative": yielded_cumulative,
                "first_yield_time": first_yield,
            },
            "convergence": {
                "failed_step_count": failed_steps,
                "converged": [1] * len(times),
                "iterations": [],
                "odb_final_time_error_s": final_time_error,
            },
            "summary": {
                "peak_top_displacement_m": ptd,
                "peak_top_displacement_time_s": ptd_t,
                "peak_top_acceleration_mps2": pta,
                "peak_top_acceleration_g": pta / 9.81 if math.isfinite(pta) else float("nan"),
                "peak_top_acceleration_time_s": pta_t,
                "peak_isolation_displacement_m": pid,
                "peak_isolation_displacement_time_s": pid_t,
                "peak_base_shear_N": pbs,
                "peak_base_shear_time_s": pbs_t,
                "first_yield_time_s": first_yield,
                "failed_step_count": failed_steps,
            },
        }
        dump_json(os.path.join(script_dir, "response_abaqus.json"), response)
        log.append(f"Extracted {len(times)} frames; final time error {final_time_error:.6e} s.")
    finally:
        odb.close()

    with open(os.path.join(script_dir, "abaqus_postprocess.log"), "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
