"""Extract results from ABAQUS ODB for Example 5."""
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
base_node = inst.nodes[0]  # node label 1
try:
    spring1_set = inst.elementSets['SPRING1']
except KeyError:
    spring1_set = None
    print('WARNING: SPRING1 element set not found')

f0 = step.frames[0]
avail_fields = [f.name for f in f0.fieldOutputs.values()]
print('Available field outputs:', avail_fields)

top_disp = np.zeros(n_frames)
top_accel = np.zeros(n_frames)
iso_disp = np.zeros(n_frames)
restoring_shear = np.zeros(n_frames)
base_reaction = np.zeros(n_frames)

for i, frame in enumerate(step.frames):
    u_field = frame.fieldOutputs['U']
    top_disp[i] = u_field.getSubset(region=top_node).values[0].data[0]
    iso_disp[i] = u_field.getSubset(region=iso_node).values[0].data[0]

    if 'A' in frame.fieldOutputs:
        top_accel[i] = frame.fieldOutputs['A'].getSubset(region=top_node).values[0].data[0]
    elif n_frames > 2 and 0 < i < n_frames - 1:
        top_accel[i] = (top_disp[i+1] - 2*top_disp[i] + top_disp[i-1]) / (dt**2)

    if 'RF' in frame.fieldOutputs:
        base_reaction[i] = frame.fieldOutputs['RF'].getSubset(region=base_node).values[0].data[0]

    if spring1_set is not None and 'S' in frame.fieldOutputs:
        s_vals = frame.fieldOutputs['S'].getSubset(region=spring1_set).values
        if len(s_vals) > 0:
            data = s_vals[0].data
            restoring_shear[i] = data[0] if hasattr(data, '__len__') else data

if 'A' not in avail_fields and n_frames > 2:
    top_accel[0] = top_accel[1]
    top_accel[-1] = top_accel[-2]

odb.close()

print('\n=== ABAQUS Results ===')
print('Peak top displacement: %.8f m' % np.max(np.abs(top_disp)))
print('Peak top acceleration: %.8f m/s^2' % np.max(np.abs(top_accel)))
print('Peak iso displacement: %.8f m' % np.max(np.abs(iso_disp)))
print('Peak restoring shear: %.4f N' % np.max(np.abs(restoring_shear)))
print('Peak base reaction: %.4f N' % np.max(np.abs(base_reaction)))

# Save results
np.savetxt('topStoDisIso2.txt', np.column_stack([time, top_disp]), fmt='%.10e')
np.savetxt('topStoAccIso2.txt', np.column_stack([time, top_accel]), fmt='%.10e')
np.savetxt('iso_displacement.txt', np.column_stack([time, iso_disp]), fmt='%.10e')
np.savetxt('base_shear.txt', np.column_stack([time, restoring_shear]), fmt='%.10e')
np.savetxt('base_reaction.txt', np.column_stack([time, base_reaction]), fmt='%.10e')
np.savetxt('hysteresis_layer1.txt',
           np.column_stack([time, iso_disp, restoring_shear]), fmt='%.10e')

response = {
    'software': 'ABAQUS 2025',
    'model_type': 'nonlinear_MDOF_T2D2_kinematic_plasticity_with_explicit_alpha_beta_dashpots',
    'n_stories': 4,
    'damping_representation': {
        'alpha_M': 'explicit nodal DASHPOT2 elements to the fixed base',
        'beta_K0': 'explicit story-to-story DASHPOT2 elements'
    },
    'results': {
        'peak_top_displacement_m': float(np.max(np.abs(top_disp))),
        'peak_top_acceleration_ms2': float(np.max(np.abs(top_accel))),
        'peak_isolation_displacement_m': float(np.max(np.abs(iso_disp))),
        'peak_base_shear_N': float(np.max(np.abs(restoring_shear))),
        'peak_base_reaction_N': float(np.max(np.abs(base_reaction)))
    }
}
with open('abaqus_response.json', 'w') as f:
    json.dump(response, f, indent=2)

print('Output files saved.')
