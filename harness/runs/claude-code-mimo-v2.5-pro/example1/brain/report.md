# Modal Analysis Cross-Validation Report

## Conclusion

The three-story shear building model was analyzed independently in MATLAB, OpenSeesPy, and ABAQUS. All three software packages produce identical natural frequencies and mode shapes within numerical precision. The cross-validation **PASSES** all consistency thresholds.

| Mode | Frequency (Hz) | Period (s) | Mode Shape (floor 1, 2, 3) |
|------|----------------|------------|---------------------------|
| 1 | 1.5838 | 0.631385 | [0.445, 0.802, 1.000] |
| 2 | 4.4378 | 0.225339 | [-1.247, -0.555, 1.000] |
| 3 | 6.4128 | 0.155939 | [1.802, -2.247, 1.000] |

## Model Assumptions

- 3-DOF shear building, each floor: m = 1000 kg, k = 500 kN/m
- Lumped mass model (mass concentrated at floor levels)
- Linear elastic behavior, no damping
- Fixed base, horizontal motion only
- SI units (N, m, s, kg) throughout

## Software Files

| Software | Script | Key Approach |
|----------|--------|--------------|
| MATLAB | `verification/matlab/main.m` | Direct eigenvalue solution of 3x3 matrices |
| OpenSeesPy | `verification/openseespy/run_openseespy_modal.py` | `elasticBeamColumn` with `fullGenLapack` solver |
| ABAQUS | `verification/abaqus/three_story_modal.inp` | `SPRING2` + `MASS` elements, Lanczos frequency extraction |

## Problems Encountered and Fixes

### OpenSeesPy

1. **Initial attempt with truss elements failed.** Vertical `truss` elements provide axial stiffness in the Y-direction (DOF 2), but mass was applied to DOF 1 (X). The model had no stiffness coupling the mass DOF, causing solver failure.

2. **Fix:** Switched to `elasticBeamColumn` elements with rotational and vertical DOFs fixed at floor nodes, leaving only horizontal translation free. Story lateral stiffness k = 12EI/L³.

3. **ARPACK solver failure.** Requesting 3 eigenvalues from a 3-DOF system fails because ARPACK requires NCV > NEV, impossible when N = NEV = 3.

4. **Fix:** Used `-fullGenLapack` direct solver instead.

5. **Incorrect I value.** Initially passed EI (41666.67) as the moment-of-inertia parameter; the element multiplied by E again, giving 10000x too-high frequencies.

6. **Fix:** Computed I = kL³/(12E) = 4.167×10⁻⁴ m⁴ correctly.

### ABAQUS

1. **Multiple keyword syntax errors.** `*CMASS1` not recognized; `*NODEOUTPUT` should be `*NODE OUTPUT`; `*FREQUENCY` parameter syntax wrong; `*SPRING` DOF specification format incorrect.

2. **SPRING2 direction issue.** SPRING2 acts along the element axis. With vertically-aligned nodes, the spring acts in Y, not X. Solution: offset nodes horizontally (X = 0, 1, 2, 3) so SPRING2 acts in the X-direction.

3. ***SPRING DOF format.** ABAQUS requires `*SPRING` data as two lines: first line = DOF pair (1, 1), second line = stiffness value.

4. **B21 beam attempt failed.** B21 is a Timoshenko beam (includes shear deformation), giving incorrect stiffness. Also, with 3D nodes, out-of-plane DOFs (3, 4, 5) had no stiffness, creating rigid body modes. Switched back to SPRING2 with horizontal node offsets.

## Consistency Thresholds

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| Frequency relative difference | < 1% | Linear eigenvalue problem should give near-identical results |
| Period relative difference | < 1% | Same as frequency |
| MAC (Modal Assurance Criterion) | > 0.99 | Confirms mode shape agreement |
| Max element-wise shape difference | < 0.05 | After sign correction |

## Final Comparison Metrics

| Pair | Freq Rel Diff | MAC | Max Shape Diff | Result |
|------|---------------|-----|----------------|--------|
| MATLAB vs OpenSeesPy | < 10⁻¹⁵ | 1.000000 | < 10⁻¹⁵ | PASS |
| MATLAB vs ABAQUS | < 10⁻⁷ | 1.000000 | < 10⁻⁷ | PASS |
| OpenSeesPy vs ABAQUS | < 10⁻⁷ | 1.000000 | < 10⁻⁷ | PASS |

MATLAB and OpenSeesPy agree to machine precision (both use direct eigenvalue solvers). ABAQUS agrees within 10⁻⁷ relative difference (Lanczos iterative solver with finite-element discretization of beam columns vs. direct spring formulation).
