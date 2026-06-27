"""Extract Abaqus ODB connector-model results into the harness response schema."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
from odbAccess import openOdb


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
GM_FILE = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"
ODB_FILE = SCRIPT_DIR / "isolated_shear_abaqus.odb"


def component_value(value, index: int = 0) -> float:
    data = value.data
    if isinstance(data, float):
        return float(data)
    return float(data[index])


def field_by_node(frame, variable: str, node_labels: list[int]) -> dict[int, float]:
    output = frame.fieldOutputs[variable]
    result: dict[int, float] = {}
    for val in output.values:
        label = int(val.nodeLabel)
        if label in node_labels:
            result[label] = component_value(val, 0)
    return result


def field_by_element(frame, preferred: list[str], element_labels: list[int]) -> tuple[str, dict[int, float]]:
    for variable in preferred:
        if variable not in frame.fieldOutputs:
            continue
        result: dict[int, float] = {}
        for val in frame.fieldOutputs[variable].values:
            label = int(val.elementLabel)
            if label in element_labels:
                result[label] = component_value(val, 0)
        if result:
            return variable, result
    raise KeyError("None of the requested connector field outputs are present: " + ", ".join(preferred))


def interpolate_ag(times: np.ndarray) -> tuple[np.ndarray, dict]:
    dt = 0.01
    g0 = 9.80665
    target_pga = 0.40 * g0
    raw = np.loadtxt(GM_FILE, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw)))
    ag = raw * (target_pga / raw_peak)
    source_time = np.arange(ag.size, dtype=float) * dt
    interp = np.interp(times, source_time, ag)
    return interp, {
        "ground_motion_file": str(GM_FILE),
        "dt": dt,
        "sample_count": int(ag.size),
        "raw_peak": raw_peak,
        "target_pga_mps2": target_pga,
        "normalized_peak_mps2": float(np.max(np.abs(ag))),
    }


def count_failed_steps() -> int:
    failed = 0
    for name in ["isolated_shear_abaqus.sta", "isolated_shear_abaqus.msg", "abaqus_run.log"]:
        path = SCRIPT_DIR / name
        if not path.exists():
            continue
        text = path.read_text(errors="ignore")
        failed += len(re.findall(r"\bfailed\b|too many attempts|analysis terminated", text, flags=re.IGNORECASE))
    return int(failed)


def integrate_loop(force: np.ndarray, displacement: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(force, displacement))
    return float(np.trapz(force, displacement))


def main() -> None:
    if not ODB_FILE.exists():
        raise RuntimeError(f"ODB file is missing: {ODB_FILE}")
    odb = openOdb(path=str(ODB_FILE), readOnly=True)
    try:
        step = odb.steps["EARTHQUAKE"]
        frames = list(step.frames)
        if not frames:
            raise RuntimeError("ODB contains no frames in EARTHQUAKE step.")

        node_labels = [1, 2, 3, 4]
        element_labels = [1, 2, 3, 4]
        times = []
        disp = []
        acc_rel = []
        story_def = []
        story_elastic_def = []
        story_force = []
        total_force_variable = None
        restoring_force_variable = None
        available_field_outputs = sorted(frames[-1].fieldOutputs.keys())

        for frame in frames:
            times.append(float(frame.frameValue))
            u_map = field_by_node(frame, "U", node_labels)
            a_map = field_by_node(frame, "A", node_labels)
            cu_var, cu_map = field_by_element(frame, ["CU"], element_labels)
            cue_var, cue_map = field_by_element(frame, ["CUE"], element_labels)
            force_var, force_map = field_by_element(frame, ["CEF", "CTF"], element_labels)
            total_var, _ = field_by_element(frame, ["CTF", "CEF"], element_labels)
            restoring_force_variable = force_var
            total_force_variable = total_var
            disp.append([u_map.get(i, math.nan) for i in node_labels])
            acc_rel.append([a_map.get(i, math.nan) for i in node_labels])
            story_def.append([cu_map.get(i, math.nan) for i in element_labels])
            story_elastic_def.append([cue_map.get(i, math.nan) for i in element_labels])
            story_force.append([force_map.get(i, math.nan) for i in element_labels])

        time = np.array(times, dtype=float)
        disp_arr = np.array(disp, dtype=float)
        acc_rel_arr = np.array(acc_rel, dtype=float)
        story_def_arr = np.array(story_def, dtype=float)
        story_elastic_def_arr = np.array(story_elastic_def, dtype=float)
        story_force_arr = np.array(story_force, dtype=float)
        ag, input_meta = interpolate_ag(time)

        stiffness = np.array([50e6, 245e6, 195e6, 98e6], dtype=float)
        masses = np.array([2.50e5, 2.70e5, 2.70e5, 1.80e5], dtype=float)
        fy = np.array([500e3, 1225e3, 975e3, 490e3], dtype=float)
        b = np.array([0.05, 0.0, 0.0, 0.0], dtype=float)
        story_plastic = story_def_arr - story_elastic_def_arr
        yielding = np.abs(story_plastic) > 1.0e-8
        top_disp = disp_arr[:, 3]
        top_acc = acc_rel_arr[:, 3] + ag
        iso_disp = disp_arr[:, 0]
        base_shear = story_force_arr[:, 0]
        loop_area = integrate_loop(base_shear, iso_disp)

        np.savetxt(SCRIPT_DIR / "topStoDisIso2.txt", np.column_stack([time, top_disp]), fmt="%.12e", delimiter="\t")
        np.savetxt(SCRIPT_DIR / "topStoAccIso2.txt", np.column_stack([time, top_acc]), fmt="%.12e", delimiter="\t")
        np.savetxt(SCRIPT_DIR / "isolation_hysteresis.txt", np.column_stack([time, iso_disp, base_shear]), fmt="%.12e", delimiter="\t")
        np.savetxt(SCRIPT_DIR / "story_hysteresis_all.txt", np.column_stack([time, story_def_arr, story_force_arr]), fmt="%.12e", delimiter="\t")

        metadata_path = SCRIPT_DIR / "abaqus_input_metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {}
        response = {
            "software": "ABAQUS",
            "generated_by": Path(__file__).name,
            "units": {"length": "m", "mass": "kg", "time": "s", "force": "N", "acceleration": "m/s^2"},
            "model": {
                "masses_kg": masses.tolist(),
                "stiffness_N_per_m": stiffness.tolist(),
                "yield_N": fy.tolist(),
                "post_yield_ratio": b.tolist(),
                "damping_ratio": 0.05,
                "rayleigh_modes": [1, 4],
                "rayleigh_alphaM": metadata.get("rayleigh_alphaM"),
                "rayleigh_betaKinit": metadata.get("rayleigh_betaKinit"),
                "damping_implementation": "physical connector dashpots: alphaM*m to fixed base node and betaKinit*k in story connectors",
            },
            "input": input_meta,
            "time": time.tolist(),
            "ground_acceleration": ag.tolist(),
            "top_story_displacement": top_disp.tolist(),
            "top_story_acceleration": top_acc.tolist(),
            "isolation_displacement": iso_disp.tolist(),
            "base_shear": base_shear.tolist(),
            "isolation_hysteresis": {
                "displacement": iso_disp.tolist(),
                "force": base_shear.tolist(),
                "loop_area": loop_area,
            },
            "story_displacement": story_def_arr.tolist(),
            "story_elastic_displacement": story_elastic_def_arr.tolist(),
            "story_force": story_force_arr.tolist(),
            "story_plastic_deformation": story_plastic.tolist(),
            "yielding": yielding.tolist(),
            "convergence": {
                "failed_step_count": count_failed_steps(),
                "frame_count": len(frames),
                "available_field_outputs": available_field_outputs,
                "restoring_force_variable": restoring_force_variable,
                "total_force_variable": total_force_variable,
            },
            "peaks": {
                "top_displacement_abs_max": float(np.nanmax(np.abs(top_disp))),
                "top_acceleration_abs_max": float(np.nanmax(np.abs(top_acc))),
                "isolation_displacement_abs_max": float(np.nanmax(np.abs(iso_disp))),
                "base_shear_abs_max": float(np.nanmax(np.abs(base_shear))),
            },
        }
        (SCRIPT_DIR / "response.json").write_text(json.dumps(response, indent=2), encoding="utf-8")
        print(f"ABAQUS extraction complete. Frames: {len(frames)}, failed step markers: {response['convergence']['failed_step_count']}")
        print(f"Restoring force output variable used: {restoring_force_variable}")
    finally:
        odb.close()


if __name__ == "__main__":
    main()
