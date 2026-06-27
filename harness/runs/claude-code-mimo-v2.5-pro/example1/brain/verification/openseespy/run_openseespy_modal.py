"""
Modal analysis of a 3-story shear building (OpenSeesPy).
Each story: lumped mass m = 1000 kg, story stiffness k = 500 kN/m.
SI units (N, m, s, kg).

Model: elasticBeamColumn elements (vertical columns) with rotational and
vertical DOFs fixed at all floor nodes, leaving only horizontal translation
as free DOF. Story lateral stiffness k = 12*EI/L^3.
"""

import openseespy.opensees as ops
import json
import math

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

m = 1000.0       # kg per floor
k = 500e3        # N/m per story
L = 1.0          # story height (m)
nStories = 3

# --- Nodes ---
ops.node(0, 0.0, 0.0)   # base
ops.node(1, 0.0, 1.0)   # floor 1
ops.node(2, 0.0, 2.0)   # floor 2
ops.node(3, 0.0, 3.0)   # floor 3

# --- Boundary conditions ---
# Fix base: all DOFs
ops.fix(0, 1, 1, 1)
# Fix vertical (DOF 2) and rotation (DOF 3) at floor nodes
# Only horizontal translation (DOF 1) is free
ops.fix(1, 0, 1, 1)
ops.fix(2, 0, 1, 1)
ops.fix(3, 0, 1, 1)

# --- Lumped masses at floor nodes (horizontal DOF = DOF 1) ---
ops.mass(1, m, 0.0, 0.0)
ops.mass(2, m, 0.0, 0.0)
ops.mass(3, m, 0.0, 0.0)

# --- Column properties ---
# Story lateral stiffness: k = 12*E*I/L^3  =>  I = k*L^3/(12*E)
E = 1.0e8                 # Young's modulus (Pa)
A = 1.0                   # cross-sectional area (m²)
I = k * L**3 / (12.0 * E) # moment of inertia (m⁴) — gives k_lateral = 12*E*I/L³ = k

ops.geomTransf('Linear', 1)

# --- Elastic beam-column elements (story springs) ---
for i in range(nStories):
    ops.element('elasticBeamColumn', i + 1, i, i + 1, A, E, I, 1)

# --- Eigenvalue analysis ---
nEigen = 3
eigenValues = ops.eigen('-fullGenLapack', nEigen)

omega = [math.sqrt(ev) for ev in eigenValues]
periods = [2 * math.pi / w for w in omega]

# --- Extract mode shapes ---
mode_shapes = []
for mode in range(nEigen):
    shape = []
    for node in range(1, nStories + 1):
        shape.append(ops.nodeEigenvector(node, mode + 1, 1))  # DOF 1 (horizontal)
    mode_shapes.append(shape)

# Normalize: top-floor amplitude = 1
for mode in range(nEigen):
    top_val = mode_shapes[mode][nStories - 1]
    if abs(top_val) < 1e-15:
        continue
    for j in range(nStories):
        mode_shapes[mode][j] /= top_val

# --- Display results ---
print("=== OpenSeesPy Modal Analysis Results ===\n")
for mode in range(nEigen):
    print(f"Mode {mode + 1}:")
    print(f"  omega  = {omega[mode]:.4f} rad/s")
    print(f"  f      = {omega[mode] / (2 * math.pi):.4f} Hz")
    print(f"  T      = {periods[mode]:.6f} s")
    shape_str = ", ".join(f"{mode_shapes[mode][j]:.6f}" for j in range(nStories))
    print(f"  phi    = [{shape_str}]\n")

# --- Export to JSON ---
results = {
    "circular_frequencies_rad_per_s": omega,
    "periods_s": periods,
    "normalized_mode_shapes": mode_shapes
}

with open("modal_results.json", "w") as f:
    json.dump(results, f, indent=2)

print("Results written to modal_results.json")
