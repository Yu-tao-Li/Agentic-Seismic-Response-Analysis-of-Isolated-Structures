from __future__ import print_function

import json
import math
from pathlib import Path

from odbAccess import openOdb


SCRIPT_DIR = Path(__file__).resolve().parent
ODB_PATH = SCRIPT_DIR / "shear_building_abaqus.odb"
META_PATH = SCRIPT_DIR / "abaqus_model_metadata.json"
GM_PATH = SCRIPT_DIR / "normalized_ground_acceleration.txt"


def read_two_column(path):
    data = []
    with open(str(path), "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append((float(parts[0]), float(parts[1])))
    return data


def interp_linear(x, y, target):
    if not x:
        raise RuntimeError("No data available for interpolation.")
    out = []
    j = 0
    last = len(x) - 1
    for t in target:
        if t <= x[0]:
            out.append(y[0])
            continue
        if t >= x[-1]:
            out.append(y[-1])
            continue
        while j < last - 1 and x[j + 1] < t:
            j += 1
        x0 = x[j]
        x1 = x[j + 1]
        y0 = y[j]
        y1 = y[j + 1]
        if abs(x1 - x0) < 1.0e-15:
            out.append(y1)
        else:
            r = (t - x0) / (x1 - x0)
            out.append(y0 + r * (y1 - y0))
    return out


def find_top_history(step):
    candidates = []
    for key, region in step.historyRegions.items():
        outputs = region.historyOutputs
        if "U1" in outputs and "A1" in outputs:
            score = 0
            key_upper = key.upper()
            if "TOP" in key_upper:
                score += 10
            if ".3" in key or " 3" in key or "N: 3" in key_upper:
                score += 5
            candidates.append((score, key, region))
    if not candidates:
        available = {}
        for key, region in step.historyRegions.items():
            available[key] = list(region.historyOutputs.keys())
        raise RuntimeError("Could not find TOP node U1/A1 history outputs. Available: {}".format(available))
    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1], candidates[0][2]


def write_two_column(path, time, values):
    with open(str(path), "w") as f:
        for t, v in zip(time, values):
            f.write("{:.12e}\t{:.12e}\n".format(t, v))


def main():
    with open(str(META_PATH), "r") as f:
        meta = json.load(f)
    gm = read_two_column(GM_PATH)
    target_time = [row[0] for row in gm]
    ag = [row[1] for row in gm]

    odb = openOdb(path=str(ODB_PATH), readOnly=True)
    try:
        step = odb.steps["EQ"]
        region_key, region = find_top_history(step)
        u_data = list(region.historyOutputs["U1"].data)
        a_data = list(region.historyOutputs["A1"].data)
    finally:
        odb.close()

    raw_u_time = [float(row[0]) for row in u_data]
    raw_u = [float(row[1]) for row in u_data]
    raw_a_time = [float(row[0]) for row in a_data]
    raw_a_rel = [float(row[1]) for row in a_data]

    top_disp = interp_linear(raw_u_time, raw_u, target_time)
    top_acc_rel = interp_linear(raw_a_time, raw_a_rel, target_time)

    if target_time and abs(target_time[0]) < 1.0e-15:
        top_disp[0] = 0.0
        top_acc_rel[0] = -ag[0]

    top_acc_abs = [a + g for a, g in zip(top_acc_rel, ag)]

    write_two_column(SCRIPT_DIR / "topStoDis.txt", target_time, top_disp)
    write_two_column(SCRIPT_DIR / "topStoAcc.txt", target_time, top_acc_abs)

    response = {
        "software": "ABAQUS",
        "model": "three_story_linear_shear_building_relative_coordinates",
        "output_definition": {
            "top_story_displacement": "relative displacement of story 3 with respect to ground, m",
            "top_story_acceleration": "absolute acceleration of story 3, m/s^2",
        },
        "time": target_time,
        "top_story_displacement": top_disp,
        "top_story_acceleration": top_acc_abs,
        "normalized_ground_acceleration": ag,
        "input_pga_after_normalization": meta["input_pga_after_normalization"],
        "dt": meta["dt"],
        "num_samples": meta["num_samples"],
        "mass_per_story_kg": meta["mass_per_story_kg"],
        "story_stiffness_N_per_m": meta["story_stiffness_N_per_m"],
        "damping_ratio": meta["damping_ratio"],
        "rayleigh_alpha_m": meta["rayleigh_alpha_m"],
        "rayleigh_beta_k": meta["rayleigh_beta_k"],
        "abaqus_history_region": region_key,
        "abaqus_raw_u_history_count": len(raw_u_time),
        "abaqus_raw_a_history_count": len(raw_a_time),
        "abaqus_raw_time_start": raw_u_time[0] if raw_u_time else None,
        "abaqus_raw_time_end": raw_u_time[-1] if raw_u_time else None,
        "postprocessing": "ODB top-node relative acceleration A1 was converted to absolute acceleration by adding normalized ground acceleration.",
    }
    with open(str(SCRIPT_DIR / "response.json"), "w") as f:
        json.dump(response, f, indent=2)

    peak_disp = max(abs(v) for v in top_disp)
    peak_acc = max(abs(v) for v in top_acc_abs)
    print("ABAQUS extraction complete from region: {}".format(region_key))
    print("Samples: {}, dt: {:.8f}, normalized PGA: {:.8f} m/s^2".format(meta["num_samples"], meta["dt"], meta["input_pga_after_normalization"]))
    print("Raw U history count: {}, raw A history count: {}".format(len(raw_u_time), len(raw_a_time)))
    print("Peak |top displacement|: {:.12g} m".format(peak_disp))
    print("Peak |top absolute acceleration|: {:.12g} m/s^2".format(peak_acc))


if __name__ == "__main__":
    main()
