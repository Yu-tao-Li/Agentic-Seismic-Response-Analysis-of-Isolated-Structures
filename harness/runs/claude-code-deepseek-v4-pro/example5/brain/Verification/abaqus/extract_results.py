"""
Extract results from ABAQUS ODB for Example 5.
Run with: abaqus python extract_results.py
"""
import numpy as np
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
os.chdir(os.path.dirname(__file__))

from odbAccess import openOdb
from abaqusConstants import *

odb_path = 'example5_model.odb'
print('Opening ODB:', odb_path)

odb = openOdb(odb_path)
step = odb.steps['DYNAMIC']

frame_times = [frame.frameValue for frame in step.frames]
n_frames = len(frame_times)
time = np.array(frame_times)
dt = 0.01
print('Number of frames:', n_frames)
print('Time range: %.4f to %.4f' % (time[0], time[-1]))

inst = odb.rootAssembly.instances['PART-1-1']
top_node = inst.nodes[4]   # node label 5
iso_node = inst.nodes[1]   # node label 2

# Check if E1 elset exists
try:
    ele1_set = inst.elementSets['E1']
    has_ele_set = True
except KeyError:
    ele1_set = None
    has_ele_set = False
    print('WARNING: E1 element set not found')

f0 = step.frames[0]
avail_fields = [f.name for f in f0.fieldOutputs.values()]
print('Available field outputs:', avail_fields)

top_disp = np.zeros(n_frames)
iso_disp = np.zeros(n_frames)
ele_force = np.zeros(n_frames)

has_ctf = 'CTF' in avail_fields
print('Has CTF field:', has_ctf)

for i, frame in enumerate(step.frames):
    u_field = frame.fieldOutputs['U']
    top_disp[i] = u_field.getSubset(region=top_node).values[0].data[0]
    iso_disp[i] = u_field.getSubset(region=iso_node).values[0].data[0]

    if has_ctf and has_ele_set:
        try:
            ctf_field = frame.fieldOutputs['CTF']
            vals = ctf_field.getSubset(region=ele1_set)
            if len(vals.values) > 0:
                ele_force[i] = vals.values[0].data[0]
        except Exception:
            pass

# Get acceleration: prefer direct ABAQUS output, fall back to central difference
has_acc = 'A' in avail_fields
print('Has A field:', has_acc)
top_accel = np.zeros(n_frames)
if has_acc:
    for i, frame in enumerate(step.frames):
        a_field = frame.fieldOutputs['A']
        top_accel[i] = a_field.getSubset(region=top_node).values[0].data[0]
else:
    if n_frames > 2:
        for i in range(1, n_frames-1):
            top_accel[i] = (top_disp[i+1] - 2*top_disp[i] + top_disp[i-1]) / (dt**2)
        top_accel[0] = top_accel[1]
        top_accel[-1] = top_accel[-2]

odb.close()

print('\n=== ABAQUS Results ===')
print('Peak top displacement: %.8f m' % np.max(np.abs(top_disp)))
print('Peak top acceleration: %.8f m/s^2' % np.max(np.abs(top_accel)))
print('Peak iso displacement: %.8f m' % np.max(np.abs(iso_disp)))
print('Peak base shear: %.4f N' % np.max(np.abs(ele_force)))

# Save results
np.savetxt('topStoDisIso2.txt', np.column_stack([time, top_disp]), fmt='%.10e')
np.savetxt('topStoAccIso2.txt', np.column_stack([time, top_accel]), fmt='%.10e')
np.savetxt('iso_displacement.txt', np.column_stack([time, iso_disp]), fmt='%.10e')
np.savetxt('base_shear.txt', np.column_stack([time, ele_force]), fmt='%.10e')
np.savetxt('hysteresis_layer1.txt',
           np.column_stack([time, iso_disp, ele_force]), fmt='%.10e')

response = {
    'software': 'ABAQUS 2025',
    'model_type': 'nonlinear_MDOF_CONN3D2_implicit',
    'n_stories': 4,
    'results': {
        'peak_top_displacement_m': float(np.max(np.abs(top_disp))),
        'peak_top_acceleration_ms2': float(np.max(np.abs(top_accel))),
        'peak_isolation_displacement_m': float(np.max(np.abs(iso_disp))),
        'peak_base_shear_N': float(np.max(np.abs(ele_force)))
    }
}
with open('abaqus_response.json', 'w') as f:
    json.dump(response, f, indent=2)

print('Output files saved.')
