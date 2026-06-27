# Three-Software Seismic Response Verification

## Scope

This run generated and executed independent analysis paths for a three-story linear shear building using the files in this run directory only. The required input record was read from `input_data/GM1.txt`, using `dt = 0.02 s` and normalized to `0.40 g = 3.924 m/s^2`.

## Model

- Stories: 3
- Mass per story: `1000 kg`
- Story stiffness: `500 kN/m`
- Equation form: relative floor displacement under base excitation, `M u_ddot + C u_dot + K u = -M r ag(t)`
- Reported displacement: top-story relative displacement, m
- Reported acceleration: top-story absolute acceleration, `u_ddot_top + ag`, m/s^2
- Rayleigh damping target: `zeta = 0.05` in modes 1 and 2
- Rayleigh coefficients: `alpha_m = 0.733397376345`, `beta_k = 0.00264307281555`

## Generated Files Used For Validation

- MATLAB: `verification/matlab/run_matlab_analysis.m`, `execution.log`, `response.json`, `topStoDis.txt`, `topStoAcc.txt`
- OpenSeesPy: `verification/openseespy/run_openseespy_analysis.py`, `execution.log`, `response.json`, `topStoDis.txt`, `topStoAcc.txt`
- ABAQUS: `verification/abaqus/create_abaqus_input.py`, `shear_building_abaqus.inp`, `extract_abaqus_results.py`, `abaqus_execution.log`, `abaqus_extract.log`, `response.json`, `topStoDis.txt`, `topStoAcc.txt`
- Cross-validation: `cross_validation_report.json`

## Execution Notes And Corrections

MATLAB executed successfully on the first run using Newmark average acceleration.

OpenSeesPy initially failed under the active Python 3.13 environment because the native `openseespywin` extension DLL could not load. The same script was rerun with the existing `<python_env>\python.exe` Python 3.12 environment, where OpenSeesPy imported correctly. A later OpenSees run failed because the Path time series could not open a file through the non-ASCII absolute path; the script was corrected to change into its local verification directory and use a relative file path. The first successful OpenSees result was too lightly damped because zeroLength elements did not contribute stiffness-proportional Rayleigh damping by default; `-doRayleigh 1` was added to each zeroLength spring element.

ABAQUS initially terminated at 2.00 s because the default 100 increment limit was too low for the 53.70 s record. The step was corrected to `inc=3000`. The dynamic procedure was also set to `alpha=0` so ABAQUS used Newmark `beta = 0.25`, `gamma = 0.5`, matching MATLAB and OpenSeesPy. The corrected ABAQUS run completed successfully, and ODB post-processing extracted 2686 top-node history samples.

## Consistency Criteria

The adopted pass thresholds are recorded in `cross_validation_report.json`:

- Time-vector max difference: `<= 1.0e-9 s`
- Peak displacement relative error: `<= 2%`
- Peak acceleration relative error: `<= 5%`
- Displacement normalized L2 error: `<= 3%`
- Acceleration normalized L2 error: `<= 5%`
- Trend correlation: `>= 0.995`

These thresholds are intentionally tight because the problem is linear and all three paths use the same mass, stiffness, damping, normalized input, time step, and response definitions. Slight tolerance is retained for ABAQUS output and solver interpolation details.

## Final Results

All outputs contain `2686` samples, `dt = 0.02 s`, and normalized input PGA `3.924 m/s^2`.

| Software | Peak top displacement (m) | Peak top absolute acceleration (m/s^2) |
| --- | ---: | ---: |
| MATLAB | 0.109023891639 | 11.684750128 |
| OpenSeesPy | 0.109021808661 | 11.684487267 |
| ABAQUS | 0.109023890389 | 11.684750642 |

Pairwise validation:

| Pair | Disp. peak err | Acc. peak err | Disp. L2 err | Acc. L2 err | Disp. corr. | Acc. corr. | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| MATLAB vs OpenSeesPy | 1.91e-05 | 2.25e-05 | 3.86e-04 | 5.66e-04 | 0.999999926 | 0.999999841 | PASS |
| MATLAB vs ABAQUS | 1.15e-08 | 4.40e-08 | 2.73e-06 | 4.09e-06 | 1.000000000 | 1.000000000 | PASS |
| OpenSeesPy vs ABAQUS | 1.91e-05 | 2.25e-05 | 3.86e-04 | 5.66e-04 | 0.999999926 | 0.999999841 | PASS |

## Final Agreement

Strict three-software pass: `PASS`.

No required generated output is missing or empty. The final pass/fail decision is based only on outputs generated in this run under `verification/matlab`, `verification/openseespy`, and `verification/abaqus`.
