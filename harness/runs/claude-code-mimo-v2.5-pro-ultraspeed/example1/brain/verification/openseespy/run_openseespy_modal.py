"""
Modal analysis of a 3-story shear building using OpenSeesPy.
Each story: lumped mass m = 1000 kg, story stiffness k = 500 kN/m = 500000 N/m
"""
import json
import numpy as np
import openseespy.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Parameters
m = 1000.0      # kg per floor
k = 500000.0    # N/m per story
num_floors = 3

# Nodes: node 0 = ground (fixed), nodes 1-3 = floors
ops.node(0, 0.0)
ops.fix(0, 1)

for i in range(1, num_floors + 1):
    ops.node(i, 0.0)
    ops.mass(i, m)

# Create elastic springs between consecutive nodes
# Material must be defined BEFORE the element that references it
for i in range(1, num_floors + 1):
    mat_tag = i
    ops.uniaxialMaterial('Elastic', mat_tag, k)
    ops.element('zeroLength', i, i - 1, i, '-mat', mat_tag, '-dir', 1)

# Modal analysis
num_modes = 3
eigenvalues = ops.eigen('-fullGenLapack', num_modes)

omega = np.sqrt(np.array(eigenvalues, dtype=float))
T = 2.0 * np.pi / omega

# Extract mode shapes (nodal displacements)
# nodeEigenvector returns a list; for 1-DOF model, take element [0]
mode_shapes = np.zeros((num_floors, num_modes))
for mode_idx in range(num_modes):
    for floor_idx in range(1, num_floors + 1):
        eigvec = ops.nodeEigenvector(floor_idx, mode_idx + 1)
        mode_shapes[floor_idx - 1, mode_idx] = float(eigvec[0]) if isinstance(eigvec, (list, tuple)) else float(eigvec)

# Normalize each mode shape to unit max absolute value
for j in range(num_modes):
    max_val = np.max(np.abs(mode_shapes[:, j]))
    if max_val > 0:
        mode_shapes[:, j] /= max_val

print("=== 3-Story Shear Building Modal Analysis (OpenSeesPy) ===\n")
for j in range(num_modes):
    print(f"Mode {j + 1}:")
    print(f"  omega = {omega[j]:.6f} rad/s")
    print(f"  T     = {T[j]:.6f} s")
    print(f"  phi   = [{mode_shapes[0, j]:.6f}, {mode_shapes[1, j]:.6f}, {mode_shapes[2, j]:.6f}]\n")

# Save JSON
results = {
    "circular_frequencies_rad_per_s": omega.tolist(),
    "periods_s": T.tolist(),
    "normalized_mode_shapes": mode_shapes.T.tolist()
}

with open("modal_results.json", "w") as f:
    json.dump(results, f, indent=2)

print("Results saved to modal_results.json")
