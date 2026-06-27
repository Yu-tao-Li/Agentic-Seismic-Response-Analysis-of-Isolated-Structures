# Three-Story Shear Building Modal Verification

## Model

The verified model is a three-story linear shear building in SI units.
Each floor has a 1000 kg lumped mass and each story has a 500000 N/m lateral stiffness.
The floor DOF order is floor 1, floor 2, floor 3.

The MATLAB mass and stiffness matrices are:

```text
M = 1000 * eye(3)
K = 500000 * [ 2 -1  0; -1  2 -1; 0 -1 1 ]
```

Mode shapes in all JSON files are 3x3 matrices with rows as floors and columns as modes. Each mode is normalized to max(abs(component)) = 1 and signed so the roof component is positive.

## Generated Outputs Used

- matlab: `harness/runs/codex-gpt-5.5-repeated/example1/brain/verification/matlab/modal_results.json`
- openseespy: `harness/runs/codex-gpt-5.5-repeated/example1/brain/verification/openseespy/modal_results.json`
- abaqus: `harness/runs/codex-gpt-5.5-repeated/example1/brain/verification/abaqus/modal_results.json`
- cross validation: `cross_validation_report.json`

## Software Implementations

- MATLAB builds the 3x3 mass and stiffness matrices directly and solves `K phi = lambda M phi`.
- OpenSeesPy uses a 1D chain of floor nodes with lumped masses and zeroLength elastic story springs.
- ABAQUS uses MASS elements at the floor nodes and SPRING2 elements in global U1, followed by Lanczos frequency extraction and ODB post-processing.

## Consistency Criteria

All three implementations represent the same three-DOF linear undamped eigenproblem in SI units. The exact model is small and well conditioned, so MATLAB and OpenSeesPy should agree to near roundoff. In this ABAQUS run the ODB API exposes modal frequency as cycles/time rounded to four decimal places, giving a first-mode relative quantization scale of about 3.2e-5. The adopted relative frequency and period tolerance of 5e-5 and mode-shape component tolerance of 1e-4 account for that documented ABAQUS output precision while remaining strict enough to catch model or unit inconsistencies.

Adopted thresholds:
- circular_frequency_relative_tolerance: `5e-05`
- period_relative_tolerance: `5e-05`
- mode_shape_component_abs_tolerance: `0.0001`
- mode_shape_mac_minimum: `0.999999`

## Errors Encountered and Corrections

- The default Python environment was Python 3.13 while the installed `openseespywin` binary imports `python312.dll`. This caused an OpenSeesPy DLL import failure during environment probing. The actual OpenSeesPy analysis was run with the matching Python 3.12 environment at `<python_env>\python.exe`; the probe is saved in `verification/openseespy/environment_probe.log`.
- ABAQUS completed successfully and produced an ODB. The `.dat` file reports warnings that requested displacement normalization was replaced by mass normalization under the default SIM architecture and that some constrained DOFs were inactive; the post-processor renormalizes extracted U1 mode shapes, and the participating floor DOF remains the intended U1 shear-building DOF.
- ABAQUS ODB frame frequencies are exposed as four-decimal cycles/time values in this run. The final tolerance was selected to account for that generated-output precision before computing the pass/fail result.

## Modal Results

### matlab

| Mode | Circular frequency (rad/s) | Period (s) | Floor 1 | Floor 2 | Floor 3 |
|---:|---:|---:|---:|---:|---:|
| 1 | 9.95143869486 | 0.631384616822 | 0.445041867913 | 0.801937735805 | 1 |
| 2 | 27.8833116047 | 0.225338560794 | -1 | -0.445041867913 | 0.801937735805 |
| 3 | 40.2925526848 | 0.155939122456 | 0.801937735805 | -1 | 0.445041867913 |

### openseespy

| Mode | Circular frequency (rad/s) | Period (s) | Floor 1 | Floor 2 | Floor 3 |
|---:|---:|---:|---:|---:|---:|
| 1 | 9.95143869486 | 0.631384616822 | 0.445041867913 | 0.801937735805 | 1 |
| 2 | 27.8833116047 | 0.225338560794 | -1 | -0.445041867913 | 0.801937735805 |
| 3 | 40.2925526848 | 0.155939122456 | 0.801937735805 | -1 | 0.445041867913 |

### abaqus

| Mode | Circular frequency (rad/s) | Period (s) | Floor 1 | Floor 2 | Floor 3 |
|---:|---:|---:|---:|---:|---:|
| 1 | 9.95130888951 | 0.631392852633 | 0.445041853346 | 0.801937722823 | 1 |
| 2 | 27.8835197562 | 0.225336878634 | -1 | -0.445041853346 | 0.801937722823 |
| 3 | 40.2928107379 | 0.155938123752 | 0.801937722823 | -1 | 0.445041853346 |

## Pairwise Metrics

| Pair | Max freq rel diff | Max period rel diff | Max mode component diff | Min MAC | Pass |
|---|---:|---:|---:|---:|---:|
| matlab vs openseespy | 5.35508e-16 | 5.27518e-16 | 3.33067e-16 | 1 | True |
| matlab vs abaqus | 1.30439e-05 | 1.30439e-05 | 1.45671e-08 | 1 | True |
| openseespy vs abaqus | 1.30439e-05 | 1.30439e-05 | 1.45671e-08 | 1 | True |

## Final Agreement

Final status: strict three-software PASS. MATLAB, OpenSeesPy, and ABAQUS generated outputs are present and agree within the adopted thresholds.
