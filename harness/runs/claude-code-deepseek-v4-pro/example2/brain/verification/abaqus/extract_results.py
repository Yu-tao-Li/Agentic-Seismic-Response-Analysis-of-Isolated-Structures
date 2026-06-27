"""Extract roof (node 4) displacement U1 and acceleration A1 from ABAQUS ODB."""
from odbAccess import openOdb
import numpy as np
import json
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
odb_path = os.path.join(script_dir, 'three_story_shear.odb')

odb = openOdb(path=odb_path)
step = odb.steps['DYNAMIC']

hr = step.historyRegions['Node PART-1-1.4']
u1 = hr.historyOutputs['U1']
a1 = hr.historyOutputs['A1']

time   = np.array([t for (t, v) in u1.data])
topDis = np.array([v for (t, v) in u1.data])
topAcc = np.array([v for (t, v) in a1.data])

odb.close()

np.savetxt(os.path.join(script_dir, 'topStoDis.txt'), np.column_stack([time, topDis]), fmt='%.10e')
np.savetxt(os.path.join(script_dir, 'topStoAcc.txt'), np.column_stack([time, topAcc]), fmt='%.10e')

peak_dis = float(np.max(np.abs(topDis)))
peak_acc = float(np.max(np.abs(topAcc)))

n = len(time)
dt_val = float(time[1] - time[0]) if n > 1 else 0.02

data = {
    "software": "ABAQUS",
    "time_step": dt_val,
    "n_samples": n,
    "input_pga_normalized_m_s2": 3.924,
    "top_displacement_peak_m": peak_dis,
    "top_acceleration_peak_m_s2": peak_acc
}

with open(os.path.join(script_dir, 'response_output.json'), 'w') as f:
    json.dump(data, f, indent=2)

print(f"Extracted {n} points, dt={dt_val:.6f}s")
print(f"Peak displacement: {peak_dis:.6e} m")
print(f"Peak acceleration: {peak_acc:.6e} m/s^2")
