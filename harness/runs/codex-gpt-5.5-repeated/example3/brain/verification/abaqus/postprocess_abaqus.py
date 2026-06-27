import json
from pathlib import Path

import numpy as np
from odbAccess import openOdb


SCRIPT_DIR = Path(__file__).resolve().parent
ODB_PATH = SCRIPT_DIR / "isolated_building.odb"
METADATA_PATH = SCRIPT_DIR / "model_metadata.json"
GM_PATH = SCRIPT_DIR / "normalized_ground_motion_mps2.txt"


def scalar_from_field(frame, field_name, node_label):
    field = frame.fieldOutputs[field_name]
    subset = field.getSubset(region=frame.frameValue and None) if False else field
    for value in subset.values:
        if value.nodeLabel == node_label:
            data = value.data
            if isinstance(data, float):
                return float(data)
            return float(data[0])
    raise KeyError(f"Node {node_label} not found in field {field_name}")


def main():
    with METADATA_PATH.open("r", encoding="utf-8") as f:
        metadata = json.load(f)

    gm = np.loadtxt(str(GM_PATH), dtype=float)
    time_input = gm[:, 0]
    ag = gm[:, 1]
    dt = float(metadata["time_step"])
    n = int(metadata["number_of_samples"])

    odb = openOdb(path=str(ODB_PATH), readOnly=True)
    try:
        step = odb.steps["GM_STEP"]
        frames = step.frames
        odb_times = np.array([float(frame.frameValue) for frame in frames], dtype=float)
        top_u = np.array([scalar_from_field(frame, "U", 5) for frame in frames], dtype=float)
        iso_u = np.array([scalar_from_field(frame, "U", 2) for frame in frames], dtype=float)
        top_rel_acc = np.array([scalar_from_field(frame, "A", 5) for frame in frames], dtype=float)
    finally:
        odb.close()

    # ABAQUS writes a frame at t=0 plus one frame at each direct-integration increment.
    # Interpolate only as a guard against harmless roundoff in frame times.
    if len(odb_times) != n or np.max(np.abs(odb_times - time_input)) > 1.0e-9:
        top_u = np.interp(time_input, odb_times, top_u)
        iso_u = np.interp(time_input, odb_times, iso_u)
        top_rel_acc = np.interp(time_input, odb_times, top_rel_acc)
        time = time_input
    else:
        time = odb_times

    top_abs_acc = top_rel_acc + ag

    np.savetxt(SCRIPT_DIR / "topStoDis.txt", np.column_stack([time, top_u]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "topStoAcc.txt", np.column_stack([time, top_abs_acc]), fmt="%.12e", delimiter="\t")
    np.savetxt(SCRIPT_DIR / "isolationDisplacement.txt", np.column_stack([time, iso_u]), fmt="%.12e", delimiter="\t")

    output = {
        "software": "ABAQUS",
        "status": "success",
        "model": {
            "units": "SI",
            "masses_kg": metadata["masses_kg"],
            "stiffnesses_N_per_m": metadata["stiffnesses_N_per_m"],
            "isolation_damping_Ns_per_m": metadata["isolation_damping_Ns_per_m"],
            "rayleigh_damping_ratio": metadata["rayleigh_damping_ratio"],
            "rayleigh_alpha_mass": metadata["rayleigh_alpha_mass"],
            "rayleigh_beta_stiffness": metadata["rayleigh_beta_stiffness"],
            "modal_frequencies_rad_per_s": metadata["modal_frequencies_rad_per_s"],
            "damping_realization": metadata["damping_realization"],
        },
        "integration": {"method": "Abaqus/Standard Dynamic Direct, HHT alpha=0", "gamma": 0.5, "beta": 0.25},
        "input": {
            "source": metadata["source_ground_motion"],
            "target_pga_mps2": metadata["target_pga_mps2"],
            "pga_after_normalization_mps2": metadata["pga_after_normalization_mps2"],
            "pga_after_normalization_g": metadata["pga_after_normalization_g"],
        },
        "time_step": dt,
        "number_of_samples": n,
        "time": time.tolist(),
        "ground_acceleration_mps2": ag.tolist(),
        "top_story_displacement_m": top_u.tolist(),
        "top_story_acceleration_mps2": top_abs_acc.tolist(),
        "top_story_relative_acceleration_mps2": top_rel_acc.tolist(),
        "isolation_layer_displacement_m": iso_u.tolist(),
        "odb_frame_count": int(len(odb_times)),
        "odb_time_start": float(odb_times[0]),
        "odb_time_end": float(odb_times[-1]),
    }
    with (SCRIPT_DIR / "response.json").open("w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    print("ABAQUS post-processing completed.")
    print(f"ODB frames: {len(odb_times)}, response samples: {len(time)}")
    print(f"Peak top displacement: {np.max(np.abs(top_u)):.12g} m")
    print(f"Peak top absolute acceleration: {np.max(np.abs(top_abs_acc)):.12g} m/s^2")
    print(f"Peak isolation displacement: {np.max(np.abs(iso_u)):.12g} m")


if __name__ == "__main__":
    main()
