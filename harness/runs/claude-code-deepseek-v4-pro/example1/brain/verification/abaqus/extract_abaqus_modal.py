"""
Extract modal analysis results from ABAQUS .dat file.
Uses regex parsing of the ABAQUS text output file.
Can be run with any Python (not restricted to ABAQUS Python).
"""
import json
import os
import re
import sys

import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
dat_path = os.path.join(script_dir, 'three_story_modal.dat')
out_path = os.path.join(script_dir, 'modal_results.json')

if not os.path.exists(dat_path):
    print(f'ERROR: .dat file not found at {dat_path}')
    sys.exit(1)

with open(dat_path, 'r', encoding='utf-8', errors='replace') as f:
    text = f.read()

# --- Extract eigenvalues (frequencies in CYCLES/TIME = Hz) ---
# Table format after "MODE NO  EIGENVALUE  FREQUENCY ...":
# "  1  99.031  9.9514  1.5838  ..."
# "  2  777.48  27.883  4.4378  ..."
# "  3  1623.5  40.293  6.4128  ..."
eigen_pat = (
    r'MODE NO\s+EIGENVALUE\s+FREQUENCY[^\n]*\n'
    r'[^\n]*\n\s*\n\s*\n'
    r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)[^\n]*\n'
    r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)[^\n]*\n'
    r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)'
)
m = re.search(eigen_pat, text)
if not m:
    print('ERROR: Could not parse eigenvalue table from .dat file')
    sys.exit(1)

freqs_cycles = [float(m.group(4)), float(m.group(8)), float(m.group(12))]

# --- Extract eigenvectors from per-mode NODE OUTPUT blocks ---
sections = re.split(r'E I G E N V A L U E\s+N U M B E R\s+\d+', text)[1:4]

modes_raw = []
for sec in sections:
    # Data lines after "NODE FOOT-  U1 ..." header
    vec_pat = (
        r'NODE FOOT-.*?U1[^\n]*\n'
        r'[^\n]*\n\s*\n'
        r'\s*(\d+)\s+([\d.E+\-]+)[^\n]*\n'
        r'\s*(\d+)\s+([\d.E+\-]+)[^\n]*\n'
        r'\s*(\d+)\s+([\d.E+\-]+)'
    )
    m2 = re.search(vec_pat, sec)
    if m2:
        modes_raw.append([float(m2.group(2)), float(m2.group(4)), float(m2.group(6))])

if len(modes_raw) < 3:
    print(f'ERROR: Only found {len(modes_raw)} mode shape sets')
    sys.exit(1)

# --- Process ---
omega = [2.0 * np.pi * f for f in freqs_cycles]
periods = [1.0 / f for f in freqs_cycles]

mode_shapes = np.array(modes_raw).T  # 3x3: rows=floors, cols=modes
for j in range(3):
    mode_shapes[:, j] /= np.max(np.abs(mode_shapes[:, j]))

# --- Display ---
print('=== ABAQUS Modal Analysis Results ===')
print(f'Circular frequencies (rad/s): {omega}')
print(f'Periods (s):                 {periods}')
print(f'Natural frequencies (Hz):    {freqs_cycles}')
print(f'Normalized mode shapes:\n{mode_shapes}')

# --- Write JSON ---
results = {
    'circular_frequencies_rad_per_s': [float(w) for w in omega],
    'periods_s': [float(p) for p in periods],
    'normalized_mode_shapes': mode_shapes.tolist(),
}

with open(out_path, 'w') as f:
    json.dump(results, f, indent=2)

print(f'\nResults saved to: {out_path}')
