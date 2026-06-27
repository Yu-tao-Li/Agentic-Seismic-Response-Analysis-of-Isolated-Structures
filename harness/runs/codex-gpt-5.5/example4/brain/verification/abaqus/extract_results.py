"""Extract results from ABAQUS ODB file."""
from odbAccess import *
from abaqusConstants import *
import numpy as np
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

odb = openOdb('model.odb')
step = odb.steps['SEISMIC']

# Extract history output
region_top = step.historyRegions['Node PART-1-1.5']  # node 5
region_iso = step.historyRegions['Node PART-1-1.2']  # node 2

u1_top = np.array([(d[0], d[1]) for d in region_top.historyOutputs['U1'].data])
a1_top = np.array([(d[0], d[1]) for d in region_top.historyOutputs['A1'].data])
u1_iso = np.array([(d[0], d[1]) for d in region_iso.historyOutputs['U1'].data])

time_vals = u1_top[:, 0]
disp_top = u1_top[:, 1]
acc_top = a1_top[:, 1]
disp_iso = u1_iso[:, 1]

# Extract base reaction from field output
rf_data = []
for frame in step.frames:
    t = frame.frameValue
    rf = frame.fieldOutputs['RF'].values[0].data[0]  # node 1, DOF 1
    # Interpolate isolation displacement at field output times
    iso_disp_at_t = np.interp(t, time_vals, disp_iso)
    rf_data.append((t, rf, iso_disp_at_t))
rf_arr = np.array(rf_data)

# Save results
np.savetxt('topStoDisIso1.txt', np.column_stack([time_vals, disp_top]), fmt='%.10e')
np.savetxt('topStoAccIso1.txt', np.column_stack([time_vals, acc_top]), fmt='%.10e')
np.savetxt('iso_displacement.txt', np.column_stack([time_vals, disp_iso]), fmt='%.10e')
np.savetxt('hysteresis_layer1.txt', rf_arr, fmt='%.10e')

print(f"Top story: {len(disp_top)} history points")
print(f"Field output: {len(rf_arr)} frames")
print(f"Max top displacement: {np.max(np.abs(disp_top)):.6f} m")
print(f"Max top acceleration: {np.max(np.abs(acc_top)):.6f} m/s^2")
print(f"Max iso displacement: {np.max(np.abs(disp_iso)):.6f} m")
print(f"Max base reaction: {np.max(np.abs(rf_arr[:,1])):.2f} N")

odb.close()
print("Done.")
