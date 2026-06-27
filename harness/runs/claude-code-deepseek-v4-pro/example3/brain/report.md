# Example 3 Cross-Software Verification Report

## Conclusion

**PASS** — MATLAB, OpenSeesPy, and ABAQUS all independently produce numerically consistent results for the same isolated four-DOF building model. All comparison metrics fall within adopted engineering thresholds. This run satisfies the strict three-software cross-validation criterion.

## Model And Input

- Model: four-DOF linear isolated shear building.
- Isolation layer: `m = 500 kg`, `k = 200 kN/m`, added viscous damping `c = 2000 N*s/m`.
- Upper stories 2 to 4: `m = 1000 kg` each, `k = 500 kN/m` each.
- Damping: Rayleigh damping with `xi = 0.05`, calibrated from modes 1 and 2.
- Computed Rayleigh coefficients: `alpha = 0.481987`, `beta = 0.00356281`.
- Ground motion: `input_data/GM1.txt`, 2686 rows, normalized to `0.40 g = 3.924 m/s^2`.
- Time step: `0.02 s`.
- Time integration: Newmark average acceleration (all three software packages).

## Software Implementations

| Software | Status | Approach | Files |
| --- | --- | --- | --- |
| MATLAB | completed | Direct Newmark integration with Rayleigh + isolation damping matrix | `verification/matlab/isolated_building.m`, `response.json`, `topStoDis.txt`, `topStoAcc.txt`, `isoDis.txt` |
| OpenSeesPy | completed | Equivalent forces via zeroLength elements, Viscous material for isolation damper | `verification/openseespy/isolated_building.py`, `response.json`, `topStoDis.txt`, `topStoAcc.txt`, `isoDis.txt` |
| ABAQUS | completed | Direct implicit dynamics (*DYNAMIC ALPHA=0), explicit DASHPOT2 elements for Rayleigh + isolation damping, SPRING2 for stiffness, CLOAD equivalent forces | `verification/abaqus/isolated_building.inp`, `response.json`, `topStoDis.txt`, `topStoAcc.txt`, `isoDis.txt` |

### Execution Environments

- **MATLAB**: Run natively. `verification/matlab/isolated_building.m`
- **OpenSeesPy**: `conda run -n ops312 python verification\openseespy\isolated_building.py`
- **ABAQUS**: Abaqus 2025 (Standard license). `abaqus job=isolated_building input=isolated_building.inp interactive`. Results extracted via `abaqus python extract_odb.py`

## Problems Encountered And Fixes

1. **ABAQUS SPRING2/DASHPOT2 data format**: Correct format requires `DOF1, DOF2` on first data line, stiffness/coefficient on second line. Previously attempted `TYPE=SPRING2` (invalid parameter), single-line `1, 200000.0` (misparsed as integer DOF2), and `1, 1, 200000.0` on one line (stiffness not found). Confirmed working format via minimal two-node frequency test.

2. **ABAQUS step increment limit**: Default max is 100 increments. For 53.70s at 0.02s (2685 increments), added `INC=3000` to `*STEP` keyword.

3. **ABAQUS output extraction**: Field output written every 100 increments (28 frames only). History output at every increment used instead, accessed via ODB history regions (`Node PART-1-1.5` for TOP, `Node PART-1-1.2` for ISO) yielding 2686 data points each.

4. **OpenSeesPy output time labeling**: The first OpenSeesPy implementation advanced the transient analysis with `analyze(1, dt)` and then stored the response against the current array index starting at `t=0`. That labeled the response at `t=dt` as if it were at `t=0`, creating a one-step phase shift. Peak values were nearly unchanged, but curve NRMSE was artificially inflated. The script now records the initial state at `t=0` and stores each post-step response at the next time sample, matching MATLAB and ABAQUS.

## Thresholds

| Metric | Threshold | Rationale |
| --- | ---: | --- |
| Displacement peak error | 2.0% | Engineering tolerance for peak response |
| Acceleration peak error | 2.0% | Engineering tolerance for peak response |
| Isolation displacement peak error | 2.0% | Engineering tolerance for peak response |
| Displacement NRMSE | 0.02 | Conservative for identical-formulation models |
| Acceleration NRMSE | 0.03 | Wider tolerance for acceleration (more sensitive to numerical differences) |

## Final Comparison

### Peak Values

| Quantity | MATLAB | OpenSeesPy | ABAQUS |
| --- | ---: | ---: | ---: |
| Peak top displacement (m) | 0.142062 | 0.142076 | 0.142062 |
| Peak top relative acceleration (m/s^2) | 7.916518 | 7.916839 | 7.916518 |
| Peak isolation displacement (m) | 0.079860 | 0.079869 | 0.079860 |

### MATLAB vs OpenSeesPy

| Metric | Value | Threshold | Pass |
| --- | ---: | ---: | --- |
| Displacement peak error | 0.0098% | 2.0% | YES |
| Acceleration peak error | 0.0041% | 2.0% | YES |
| Isolation disp. peak error | 0.0107% | 2.0% | YES |
| Displacement NRMSE | 0.0000765 | 0.02 | YES |
| Acceleration NRMSE | 0.000104 | 0.03 | YES |

### MATLAB vs ABAQUS

| Metric | Value | Threshold | Pass |
| --- | ---: | ---: | --- |
| Displacement peak error | 0.0000048% | 2.0% | YES |
| Acceleration peak error | 0.0000014% | 2.0% | YES |
| Isolation disp. peak error | 0.0000032% | 2.0% | YES |
| Displacement NRMSE | 2.27e-09 | 0.02 | YES |
| Acceleration NRMSE | 1.85e-09 | 0.03 | YES |

> MATLAB and ABAQUS agree to machine precision: both use mathematically identical direct Newmark integration, equivalent force loading, and explicit Rayleigh damping calculated from the same alpha/beta coefficients.

### OpenSeesPy vs ABAQUS

| Metric | Value | Threshold | Pass |
| --- | ---: | ---: | --- |
| Displacement peak error | 0.0098% | 2.0% | YES |
| Acceleration peak error | 0.0041% | 2.0% | YES |
| Isolation disp. peak error | 0.0107% | 2.0% | YES |
| Displacement NRMSE | 0.0000765 | 0.02 | YES |
| Acceleration NRMSE | 0.000104 | 0.03 | YES |

### Reference Comparison

Run MATLAB results match origin-extracted reference: displacement difference ~0.009%, acceleration difference ~0.003%.

## Verdict

**`overall_pass = true`** — **PASS**

All three independent structural analysis packages (MATLAB, OpenSeesPy, ABAQUS) produce numerically consistent results for the Example 3 isolated four-DOF building model. MATLAB and ABAQUS agree to machine precision via identical mathematical formulations. After correcting OpenSeesPy's output time convention, the maximum pairwise peak error is about 0.0107% and the maximum pairwise curve NRMSE is about 0.0104%. The earlier large curve error was a one-step output time-labeling artifact, not a physical model discrepancy. All metrics fall within conservatively adopted engineering thresholds. This run demonstrates auditable three-software cross-validation.
