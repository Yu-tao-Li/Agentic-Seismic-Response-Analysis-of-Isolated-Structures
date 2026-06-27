"""
Extract modal analysis results from ABAQUS ODB file.
Run with: abaqus python extract_abaqus_modal.py
"""

from odbAccess import openOdb
import json
import math
import sys

odb_path = 'three_story_modal.odb'

print("=== ABAQUS Modal Results Extraction ===\n")

odb = openOdb(path=odb_path, readOnly=True)
step = odb.steps['FREQ']

# Get eigenvalues from history output EIGFREQ (returns frequencies in Hz)
historyRegion = step.historyRegions['Assembly Assembly-1']
freq_hz = [row[1] for row in historyRegion.historyOutputs['EIGFREQ'].data]

# Compute circular frequencies and periods
omega = [2 * math.pi * f for f in freq_hz]
periods = [1.0 / f for f in freq_hz]

print("Frequencies (Hz): {}".format(["{:.4f}".format(f) for f in freq_hz]))
print("Circular frequencies (rad/s): {}".format(["{:.4f}".format(w) for w in omega]))
print("Periods (s): {}".format(["{:.6f}".format(T) for T in periods]))

# Extract mode shapes from displacement output
# Floor nodes are 2, 3, 4 (DOF 1 = X-direction)
nModes = len(freq_hz)
mode_shapes = []

for modeIdx in range(nModes):
    frame = step.frames[modeIdx + 1]  # frames 1,2,3 (frame 0 is base state)
    displacement = frame.fieldOutputs['U']

    shape = []
    for nodeLabel in [2, 3, 4]:
        found = False
        for value in displacement.values:
            if value.nodeLabel == nodeLabel:
                shape.append(value.data[0])  # X displacement (DOF 1)
                found = True
                break
        if not found:
            shape.append(0.0)

    mode_shapes.append([float(x) for x in shape])

# Normalize mode shapes (top-floor amplitude = 1)
for mode in range(nModes):
    top_val = mode_shapes[mode][2]  # floor 3 is index 2
    if abs(top_val) < 1e-15:
        print("Warning: Mode {} has near-zero top displacement".format(mode + 1))
        continue
    for j in range(3):
        mode_shapes[mode][j] /= top_val

print("\nNormalized mode shapes (floor 1, 2, 3):")
for mode in range(nModes):
    print("  Mode {}: [{:.6f}, {:.6f}, {:.6f}]".format(
        mode + 1, mode_shapes[mode][0], mode_shapes[mode][1], mode_shapes[mode][2]))

# Export to JSON
results = {
    "circular_frequencies_rad_per_s": omega,
    "periods_s": periods,
    "normalized_mode_shapes": mode_shapes
}

with open("modal_results.json", "w") as f:
    json.dump(results, f, indent=2)

print("\nResults written to modal_results.json")

odb.close()
