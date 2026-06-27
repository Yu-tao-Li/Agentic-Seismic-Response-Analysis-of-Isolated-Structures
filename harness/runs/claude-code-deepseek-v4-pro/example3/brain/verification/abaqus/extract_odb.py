"""Extract results from ABAQUS ODB history output."""
import sys
sys.path.insert(0, r'<abaqus_install>\EstProducts\2025\win_b64\code\python\lib')

from odbAccess import openOdb
import numpy as np
import json

odb = openOdb('isolated_building.odb')

step = odb.steps['DYNAMIC']

# History regions use auto-generated names: Node PART-1-1.<nodelabel>
# Node 5 = TOP (has U1, A1), Node 2 = ISO (has U1)
top_region = step.historyRegions['Node PART-1-1.5']
iso_region = step.historyRegions['Node PART-1-1.2']

# Extract history data: .data returns tuple of (time, value) pairs
u1_top_data = top_region.historyOutputs['U1'].data
a1_top_data = top_region.historyOutputs['A1'].data
u1_iso_data = iso_region.historyOutputs['U1'].data

n = len(u1_top_data)
t = np.array([p[0] for p in u1_top_data])
top_dis = np.array([p[1] for p in u1_top_data])
top_acc = np.array([p[1] for p in a1_top_data])
iso_dis = np.array([p[1] for p in u1_iso_data])

# Load ground motion for absolute acceleration
gm = np.loadtxt('../../input_data/GM1.txt')
g = 9.81
scale = 0.40 * g / np.max(np.abs(gm))
ag = gm * scale

# ABAQUS A is relative acceleration (base is fixed, forces applied as CLOAD)
# Verify: A1 should be relative acceleration
top_acc_abs = top_acc + ag[:n]

print(f"Data points: {n}")
print(f"Peak top displacement: {np.max(np.abs(top_dis)):.6f} m")
print(f"Peak top rel acceleration: {np.max(np.abs(top_acc)):.6f} m/s^2")
print(f"Peak top abs acceleration: {np.max(np.abs(top_acc_abs)):.6f} m/s^2")
print(f"Peak isolation displacement: {np.max(np.abs(iso_dis)):.6f} m")

# Write output
np.savetxt('topStoDis.txt', np.column_stack([t, top_dis]), fmt='%.10e')
np.savetxt('topStoAcc.txt', np.column_stack([t, top_acc]), fmt='%.10e')
np.savetxt('isoDis.txt', np.column_stack([t, iso_dis]), fmt='%.10e')

resp = {
    'time': t.tolist(),
    'top_displacement': top_dis.tolist(),
    'top_rel_acceleration': top_acc.tolist(),
    'top_abs_acceleration': top_acc_abs.tolist(),
    'isolation_displacement': iso_dis.tolist(),
    'input_pga': float(0.40 * g),
    'dt': 0.02,
    'n_samples': int(n),
    'peak_top_dis': float(np.max(np.abs(top_dis))),
    'peak_top_rel_acc': float(np.max(np.abs(top_acc))),
    'peak_top_abs_acc': float(np.max(np.abs(top_acc_abs))),
    'peak_iso_dis': float(np.max(np.abs(iso_dis))),
}

with open('response.json', 'w') as f:
    json.dump(resp, f, indent=2)

odb.close()
print("Extraction complete.")
