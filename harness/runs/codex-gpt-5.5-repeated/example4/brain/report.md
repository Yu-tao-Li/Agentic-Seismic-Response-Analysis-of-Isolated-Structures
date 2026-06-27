# Cross-Validation Report

## Scope
Generated and executed three independent analysis paths for the four-story isolated nonlinear shear-building model using the supplied `Northridge_01_NO_968.txt` record, `dt = 0.01 s`, and PGA normalized to `0.40 g`.

## Model and Numerical Assumptions
- Masses, stiffnesses, yield forces, and post-yield ratios match the fixed prompt values in all generated inputs.
- Damping uses `zeta = 0.05` Rayleigh coefficients from initial modes 1 and 4, implemented as `alpha*M + beta*K_initial`; the isolation dashpot `c = 1000e3 N*s/m` is added explicitly.
- MATLAB uses average-acceleration Newmark with Newton iteration and return mapping for ideal elastic-plastic stories.
- OpenSeesPy uses zeroLength springs with Elastic/ElasticPP materials and explicit viscous elements for the same damping matrix.
- ABAQUS uses lumped mass elements, T2D2 truss spring elements with elastic-perfectly-plastic material behavior, DASHPOT2 elements for damping, and `*DYNAMIC, ALPHA=0.0` to match the Newmark average-acceleration limit. A negligible density `1.0e-12 kg/m^3` was added to truss materials because ABAQUS/Standard requires at least one active material density in general dynamic steps; inertia remains governed by the specified lumped masses.

## Execution Notes and Corrections
- MATLAB completed with zero failed steps.
- The default Python 3.13 environment contained an OpenSeesPy package linked to `python312.dll`, causing an import failure. The analysis was rerun successfully under `<python_env>/python.exe` with OpenSeesPy 3.8.0.
- ABAQUS first failed preprocessing due to zero active material density. After adding negligible truss density and setting `ALPHA=0.0` for Newmark-equivalent integration, ABAQUS completed successfully with 1899 increments and zero cutbacks. The wrapper was updated to verify success from the `.sta` file before postprocessing.

## Consistency Criteria
- Peak response relative error limit: 1.0%.
- Normalized RMS curve error limit: 2.0%; normalized max curve error limit: 6.0%.
- Isolation peak relative error limit: 1.0%.
- First-yield time tolerance: 0.02 s; failed step count must be zero.
- Reasoning: MATLAB and OpenSeesPy use average-acceleration Newmark and should match nearly exactly. ABAQUS/Standard is configured with HHT alpha = 0.0, the Newmark average-acceleration limit. The retained 1% peak and 2% normalized RMS curve limits are tight enough to catch modeling errors while allowing small solver implementation and ODB interpolation differences.

## Generated Outputs Used
- matlab: `verification/matlab/response_matlab.json`, `verification/matlab/topStoDisIso1.txt`, `verification/matlab/topStoAccIso1.txt`; time points = 1900, failed steps = 0.
- openseespy: `verification/openseespy/response_openseespy.json`, `verification/openseespy/topStoDisIso1.txt`, `verification/openseespy/topStoAccIso1.txt`; time points = 1900, failed steps = 0.
- abaqus: `verification/abaqus/response_abaqus.json`, `verification/abaqus/topStoDisIso1.txt`, `verification/abaqus/topStoAccIso1.txt`; time points = 1900, failed steps = 0.

## Peak Response Summary
| Software | Peak top disp (m) | Peak top acc (m/s^2) | Peak iso disp (m) | Peak base shear (N) | First yield times s2/s3/s4 |
|---|---:|---:|---:|---:|---|
| matlab | 0.03589847 | 2.78979538 | 0.02504130 | 1252065.167 | none/none/6.580 |
| openseespy | 0.03589846 | 2.78979545 | 0.02504130 | 1252065.185 | none/none/6.580 |
| abaqus | 0.03589847 | 2.78979532 | 0.02504130 | 1252065.125 | none/none/6.580 |

## Pairwise Metrics
| Pair | Max peak rel error | Max RMS curve norm | Max curve norm | Yield max time err (s) | Pair pass |
|---|---:|---:|---:|---:|---|
| matlab_vs_openseespy | 1.287285e-07 | 4.141132e-05 | 1.353491e-03 | 0.000000e+00 | True |
| matlab_vs_abaqus | 3.353020e-08 | 2.153414e-06 | 6.766599e-05 | 7.629395e-08 | True |
| openseespy_vs_abaqus | 1.319112e-07 | 4.142517e-05 | 1.353756e-03 | 7.629395e-08 | True |

## Final Verdict
Strict three-software pass: **True**.
All required generated outputs are present and pairwise metrics meet adopted thresholds.

Detailed machine-readable metrics are in `cross_validation_report.json`.
