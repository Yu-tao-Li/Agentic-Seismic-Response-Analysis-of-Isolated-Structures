"""
Post-process ABAQUS ODB to extract modal frequencies and mode shapes.
"""
from abaqus import session
import json
import math

odb_path = 'three_story_modal.odb'
odb = session.openOdb(name=odb_path)
step = odb.steps['FREQ']

frames = step.frames
print("Number of frames:", len(frames))

# Extract frequencies from frames (frame 0 is base state, frames 1-3 are modes)
freqs_hz = []
for i in range(1, len(frames)):
    frame = frames[i]
    freq_hz = frame.frequency  # cycles/time
    freqs_hz.append(freq_hz)
    print(f"Mode {i}: freq = {freq_hz:.4f} Hz")

num_modes = len(freqs_hz)
circular_freqs = [2.0 * math.pi * f for f in freqs_hz]
periods = [1.0 / f if f > 0 else 0.0 for f in freqs_hz]

# Extract displacement output for mode shapes
floor_nodes = [2, 3, 4]

mode_shapes = []
for mode_idx in range(num_modes):
    frame = frames[mode_idx + 1]
    disp_field = frame.fieldOutputs['U']

    mode_shape = []
    for node_label in floor_nodes:
        for val in disp_field.values:
            if val.nodeLabel == node_label:
                mode_shape.append(float(val.data[0]))
                break
        else:
            mode_shape.append(0.0)

    mode_shapes.append(mode_shape)

# Normalize each mode shape (max absolute value = 1)
for j in range(num_modes):
    max_val = max(abs(v) for v in mode_shapes[j])
    if max_val > 0:
        mode_shapes[j] = [v / max_val for v in mode_shapes[j]]

print("\n=== 3-Story Shear Building Modal Analysis (ABAQUS) ===\n")
for j in range(num_modes):
    print(f"Mode {j + 1}:")
    print(f"  omega = {circular_freqs[j]:.6f} rad/s")
    print(f"  T     = {periods[j]:.6f} s")
    print(f"  phi   = [{mode_shapes[j][0]:.6f}, {mode_shapes[j][1]:.6f}, {mode_shapes[j][2]:.6f}]\n")

results = {
    "circular_frequencies_rad_per_s": [round(x, 10) for x in circular_freqs],
    "periods_s": [round(x, 10) for x in periods],
    "normalized_mode_shapes": [[round(x, 10) for x in row] for row in mode_shapes]
}

with open("modal_results.json", "w") as f:
    json.dump(results, f, indent=2)

print("Results saved to modal_results.json")
