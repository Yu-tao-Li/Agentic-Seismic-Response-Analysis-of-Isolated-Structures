# Example 1 — 3-Story Shear Building Modal Analysis

## Conclusion

All three software implementations (MATLAB, OpenSeesPy, ABAQUS) produce consistent results. Circular frequencies agree to within 1.3×10⁻⁵ relative error, and mode shapes achieve MAC = 1.0 across all pairs.

## Model

- 3-story shear building, 1 DOF per floor (lateral translation)
- Lumped mass per floor: 1000 kg
- Story stiffness: 500 kN/m = 500,000 N/m
- SI units throughout

| Mode | ω (rad/s) | T (s) |
|------|-----------|-------|
| 1    | 9.9514    | 0.6314 |
| 2    | 27.8833   | 0.2253 |
| 3    | 40.2926   | 0.1559 |

## Software files

- `verification/matlab/main.m` — Analytical eigenvalue solution via `eig(K,M)`
- `verification/openseespy/run_openseespy_modal.py` — ZeroLength elastic springs + lumped mass; `eigen('-fullGenLapack')`
- `verification/abaqus/three_story_modal.inp` — SPRING2 elements + MASS elements; Lanczos frequency extraction

## Problems encountered and fixes

1. **OpenSeesPy ARPACK failure**: Requesting 3 eigenvalues from a 3-DOF system failed (NCV must exceed NEV). Fix: switched to `-fullGenLapack`.
2. **OpenSeesPy `nodeEigenvector` return type**: Returns a list in newer versions, not a scalar. Fix: index `[0]`.
3. **OpenSeesPy JSON transpose bug**: `mode_shapes` stored modes as columns; `tolist()` transposed the matrix. Fix: use `mode_shapes.T.tolist()`.
4. **ABAQUS `*MASS` syntax**: `*Mass, elset=` requires an element set of MASS elements, not a node set. Fix: created `*Element, type=MASS` elements and referenced them via `*Elset`.
5. **ABAQUS `*SPRING` data line**: SPRING2 requires a DOF specification line before the stiffness value. Fix: added `1, 1` line.
6. **ABAQUS `*INERTIA`**: Misinterpreted as `*InertiaRelief`. Replaced with MASS element approach.
7. **MATLAB JSON output**: Wrote literal `\n` instead of newlines. Fix: rewrote JSON file directly.

## Consistency criteria

- Frequency: relative error < 0.1% — covers floating-point and solver-level differences
- Mode shape: MAC ≥ 0.99 — ensures shapes are essentially identical (signs arbitrary)

## Final comparison

| Pair | Max freq rel. error | Min MAC |
|------|-------------------|---------|
| MATLAB vs OpenSeesPy | 0 | 1.0 |
| MATLAB vs ABAQUS | 1.3×10⁻⁵ | 1.0 |
| OpenSeesPy vs ABAQUS | 1.3×10⁻⁵ | 1.0 |

**Overall: PASS**
