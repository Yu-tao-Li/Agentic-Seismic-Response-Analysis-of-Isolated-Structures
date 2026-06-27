# Cross-Software Verification Report

## Scope

This run generated and executed independent verification paths under:

- `verification/matlab/`
- `verification/openseespy/`
- `verification/abaqus/`

The final comparison used only generated outputs from this run. No reference result files from other runs were used as generated software results.

## Model

All solvers used SI units and a 4-DOF shear-building model in relative ground coordinates. The four DOFs are the isolation mass and stories 2 to 4. The earthquake input was read from `input_data/Northridge_01_NO_968.txt`, sampled at `dt = 0.01 s`, and scaled to a peak ground acceleration of `0.40 g = 3.92266 m/s^2`.

Fixed properties:

- Masses: `[2.50e5, 2.70e5, 2.70e5, 1.80e5] kg`
- Story stiffnesses: `[50e6, 245e6, 195e6, 98e6] N/m`
- Yield forces: `[500e3, 1225e3, 975e3, 490e3] N`
- Post-yield ratios: `[0.05, 0, 0, 0]`

Rayleigh damping used `xi = 0.05` in elastic modes 1 and 4. The initial-stiffness coefficients were:

- `alphaM = 0.580729296596`
- `betaKinit = 0.00173379820464`

Top-story acceleration is reported as absolute acceleration: relative structural acceleration plus ground acceleration.

## Software Paths

MATLAB:

- Source: `verification/matlab/run_matlab_analysis.m`
- Method: Newmark average acceleration with Newton iteration and return-mapping material update.
- Output: `response.json`, `topStoDisIso2.txt`, `topStoAccIso2.txt`, `isolation_hysteresis.txt`, `story_hysteresis_all.txt`
- Log: `verification/matlab/matlab_run.log`

OpenSeesPy:

- Source: `verification/openseespy/run_openseespy_analysis.py`
- Method: 1D zeroLength shear chain with `Steel01` materials and uniform excitation.
- Python used: `<python_env>\python.exe`, because the active Python 3.13 environment had an incompatible OpenSeesPy native extension.
- Output: `response.json`, `topStoDisIso2.txt`, `topStoAccIso2.txt`, `isolation_hysteresis.txt`, `story_hysteresis_all.txt`
- Log: `verification/openseespy/openseespy_run.log`

ABAQUS:

- Generator: `verification/abaqus/write_abaqus_input.py`
- Input: `verification/abaqus/isolated_shear_abaqus.inp`
- Extractor: `verification/abaqus/extract_abaqus_results.py`
- Method: point masses, `CONN3D2` connector springs with connector plasticity/kinematic hardening, and physical connector dashpots matching `alphaM*M + betaK*K0`.
- Output: `response.json`, `topStoDisIso2.txt`, `topStoAccIso2.txt`, `isolation_hysteresis.txt`, `story_hysteresis_all.txt`
- Logs: `abaqus_input_generation.log`, `abaqus_run.log`, `abaqus_extract.log`, plus `.dat`, `.msg`, `.sta`

## Errors and Corrections

- OpenSeesPy initially failed under the active Python 3.13 environment because the native extension could not load. The OpenSees-specific conda environment `<python_env_root>\envs\opensees` imported OpenSeesPy successfully and was used for execution.
- OpenSeesPy postprocessing initially used `np.trapz`, which is absent in the installed NumPy 2.4 environment. It was changed to `np.trapezoid`.
- MATLAB initially wrote the assembled first-DOF restoring force as base shear. The output was corrected to the local isolation restoring force so it is physically comparable to OpenSeesPy and ABAQUS isolation hysteresis.
- ABAQUS first input attempt used node label `0` and requested invalid connector output `CVF`; Abaqus rejected preprocessing. The base node was changed to positive label `1000`, and invalid output requests were removed.
- ABAQUS extraction initially used `np.trapezoid`, unavailable in Abaqus Python's older NumPy. A compatibility helper now uses `np.trapz` when needed.
- ABAQUS yield indicators were changed to connector plastic motion `CU - CUE` instead of a force-derived estimate, avoiding false superstructure-yield flags from small connector output interpolation/roundoff.

## Criteria

The adopted thresholds were recorded before final pass/fail in `cross_validation_report.json`:

- Displacement and isolation/base-shear peak relative error: `5%`
- Top acceleration peak relative error: `10%`
- Displacement curve RMS relative error: `5%`
- Acceleration curve RMS relative error: `15%`
- Hysteresis loop area relative error: `10%`
- First-yield time difference: `0.02 s`

These tolerances reflect independent nonlinear integration and ABAQUS connector/output interpolation while still requiring close physical agreement.

## Results

All three paths produced non-empty machine-readable outputs and manuscript-style files. All convergence/failed-step counts are zero.

Peak responses:

| Software | Top disp. max (m) | Top abs. accel. max (m/s^2) | Isolation disp. max (m) | Base shear max (N) |
|---|---:|---:|---:|---:|
| MATLAB | 0.0274721664 | 2.5604880635 | 0.0218942374 | 529735.5935 |
| OpenSeesPy | 0.0274721407 | 2.5604995651 | 0.0218942140 | 529735.5351 |
| ABAQUS | 0.0277290754 | 2.5656342581 | 0.0220415965 | 530104.0000 |

Pairwise metric highlights:

| Pair | Top disp peak err | Top accel peak err | Top disp RMS err | Top accel RMS err | Hysteresis area err | Pass |
|---|---:|---:|---:|---:|---:|---|
| MATLAB vs OpenSeesPy | 9.35e-7 | 4.49e-6 | 1.96e-4 | 4.65e-4 | 2.49e-5 | true |
| MATLAB vs ABAQUS | 9.24e-3 | 1.98e-3 | 1.73e-2 | 5.77e-2 | 1.15e-2 | true |
| OpenSeesPy vs ABAQUS | 9.24e-3 | 1.98e-3 | 1.73e-2 | 5.77e-2 | 1.15e-2 | true |

Yielding sequence:

- Isolation layer first yielded at about `3.68 s` in all three solvers.
- Stories 2 to 4 did not yield in any solver.

## Final Status

`cross_validation_report.json` reports:

- `status`: `PASS`
- `strict_three_software_pass`: `true`

The required MATLAB, OpenSeesPy, and ABAQUS outputs are present, non-empty, generated in this run, physically comparable, and within the adopted thresholds.
