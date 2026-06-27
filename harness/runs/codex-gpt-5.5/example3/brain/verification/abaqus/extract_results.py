from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from odbAccess import openOdb

out = Path(__file__).resolve().parent
odb = openOdb(path=str(out / "model.odb"), readOnly=True)
try:
    step = odb.steps["DYNAMIC"]
    top_label = 5
    time = []
    top_u = []
    top_a = []
    iso_u = []
    for frame in step.frames:
        time.append(float(frame.frameValue))
        u_field = frame.fieldOutputs["U"]
        a_field = frame.fieldOutputs["A"]
        u_by_label = {v.nodeLabel: v for v in u_field.values}
        a_by_label = {v.nodeLabel: v for v in a_field.values}
        top_u.append(float(u_by_label[top_label].data[1]))
        top_a.append(float(a_by_label[top_label].data[1]))
        iso_u.append(float(u_by_label[2].data[1]))
finally:
    odb.close()
time = np.asarray(time, dtype=float)
top_u = np.asarray(top_u, dtype=float)
top_a = np.asarray(top_a, dtype=float)
iso_u = np.asarray(iso_u, dtype=float)
np.savetxt(out / "topStoDis.txt", np.column_stack([time, top_u]), fmt="%.10e")
np.savetxt(out / "topStoAcc.txt", np.column_stack([time, top_a]), fmt="%.10e")
np.savetxt(out / "isolation_displacement.txt", np.column_stack([time, iso_u]), fmt="%.10e")
with open(out / "response_output.json", "w", encoding="utf-8") as f:
    json.dump({
        "software": "ABAQUS",
        "time": time.tolist(),
        "top_displacement": top_u.tolist(),
        "top_acceleration": top_a.tolist(),
        "isolation_displacement": iso_u.tolist(),
        "peak_top_displacement_m": float(np.max(np.abs(top_u))),
        "peak_top_acceleration_m_s2": float(np.max(np.abs(top_a))),
        "peak_isolation_displacement_m": float(np.max(np.abs(iso_u))),
        "n_samples": int(len(time)),
    }, f, indent=2)
