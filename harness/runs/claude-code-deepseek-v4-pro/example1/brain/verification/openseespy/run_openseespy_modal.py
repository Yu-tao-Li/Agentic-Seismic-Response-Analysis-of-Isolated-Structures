"""
Three-story shear building modal analysis using OpenSeesPy.
Model: lumped mass 1000 kg per floor, story stiffness 500 kN/m.
SI units, 3D model with 6 DOFs per node, only X translation free at floors.
"""
import json
import os
import sys

import numpy as np

import openseespy.opensees as ops

m_floor = 1000.0     # kg
k_story = 500000.0   # N/m (500 kN/m)

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# --- Define nodes along X axis ---
ops.node(1, 0.0, 0.0, 0.0)   # base
ops.node(2, 1.0, 0.0, 0.0)   # 1st floor
ops.node(3, 2.0, 0.0, 0.0)   # 2nd floor
ops.node(4, 3.0, 0.0, 0.0)   # 3rd floor (roof)

# --- Boundary conditions ---
ops.fix(1, 1, 1, 1, 1, 1, 1)           # base: fully fixed
for nd in [2, 3, 4]:
    ops.fix(nd, 0, 1, 1, 1, 1, 1)      # floors: only X free

# --- Lumped masses (X direction only, other DOFs zero) ---
for nd in [2, 3, 4]:
    ops.mass(nd, m_floor, 0.0, 0.0, 0.0, 0.0, 0.0)

# --- Elastic spring material ---
ops.uniaxialMaterial('Elastic', 1, k_story)

# --- twoNodeLink elements: spring in DOF 1 (X) between floors ---
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.element('twoNodeLink', 2, 2, 3, '-mat', 1, '-dir', 1)
ops.element('twoNodeLink', 3, 3, 4, '-mat', 1, '-dir', 1)

# --- Eigenvalue analysis (fullGenLapack avoids ARPACK NCV constraint) ---
eigenvalues = ops.eigen('-fullGenLapack', 3)

omega = np.sqrt(np.array(eigenvalues))          # rad/s
periods = 2.0 * np.pi / omega                    # s

# --- Extract mode shapes (X DOF at floor nodes) ---
mode_shapes = np.zeros((3, 3))
for mode_num in range(1, 4):
    for i, nd in enumerate([2, 3, 4]):
        ev = ops.nodeEigenvector(nd, mode_num)
        mode_shapes[i, mode_num - 1] = ev[0]    # DOF 1 = X translation

# --- Normalize: max absolute = 1 ---
for j in range(3):
    mode_shapes[:, j] /= np.max(np.abs(mode_shapes[:, j]))

# --- Display ---
print('=== OpenSeesPy Modal Analysis Results ===')
print(f'Circular frequencies (rad/s): {omega}')
print(f'Periods (s):                 {periods}')
print(f'Natural frequencies (Hz):    {omega / (2*np.pi)}')
print(f'Normalized mode shapes:\n{mode_shapes}')

# --- Write JSON ---
results = {
    'circular_frequencies_rad_per_s': omega.tolist(),
    'periods_s': periods.tolist(),
    'normalized_mode_shapes': mode_shapes.tolist(),
}

script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, 'modal_results.json')
with open(out_path, 'w') as f:
    json.dump(results, f, indent=2)

print(f'\nResults saved to: {out_path}')
