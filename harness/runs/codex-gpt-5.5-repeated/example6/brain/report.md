# Example 6 Cross-Software Verification Report

## Outcome

Strict three-software result: **PASS**.

MATLAB, OpenSeesPy, and ABAQUS all ran in this run directory. The earlier ABAQUS compiler-blocked conclusion was rechecked and overturned: Visual Studio 2022 BuildTools plus Intel oneAPI successfully compile the Bouc-Wen user subroutine, and ABAQUS/Standard completes with no analysis warnings or errors.

## Outputs Used

- MATLAB: `verification/matlab/matlab_results.json`
- OpenSeesPy: `verification/openseespy/opensees_results.json`
- ABAQUS: `verification/abaqus/abaqus_results.json`
- Cross-validation: `cross_validation_report.json`

No outputs were copied from other runs.

## ABAQUS Compiler Evidence

The ABAQUS wrapper command loads the local build environment:

```text
set "VS2022INSTALLDIR=<visual_studio_buildtools>" &&
call "<visual_studio_buildtools>\VC\Auxiliary\Build\vcvars64.bat" &&
call "<intel_oneapi>\setvars.bat" intel64 --force &&
"abaqus" ...
```

`verification/abaqus/abaqus_make.log` shows Visual Studio and oneAPI initialized, Intel Fortran Classic `ifort` found, and `Abaqus JOB boucwen_uel.for COMPLETED`.

`verification/abaqus/abaqus_job.log` shows ABAQUS 2025 compiling/linking the user subroutine, running Standard, and `Abaqus JOB example6_boucwen COMPLETED`.

`verification/abaqus/example6_boucwen.msg` reports:

```text
THE ANALYSIS HAS BEEN COMPLETED
0 WARNING MESSAGES DURING ANALYSIS
0 ERROR MESSAGES
```

## ABAQUS Model

The final ABAQUS model uses compiled user elements:

- `U1`: two-node horizontal Bouc-Wen isolator UEL.
- `U2`: 30-DOF reduced elastic-frame UEL generated from the same transformed Euler-Bernoulli stiffness used by MATLAB.
- Parallel `DASHPOT2` elements provide the isolator initial-stiffness damping term.
- Floor masses are assigned to the three horizontal floor generalized coordinates.
- Ground motion is applied as equivalent inertial loads with endpoint timing aligned to the MATLAB/OpenSeesPy Newmark integration.
- ODB postprocessing computes roof displacement, drifts, isolation displacement, Bouc-Wen material force, reconstructed `z`, base shear, peaks, and CSV histories.

## Key Metrics

MATLAB vs OpenSeesPy: **PASS**.

- Peak roof displacement relative difference: `1.077477191544167e-08`
- Peak isolation displacement relative difference: `1.4293817641999786e-08`
- Peak total base shear relative difference: `2.0891628003337892e-10`
- Roof displacement NRMSE: `7.255510740701096e-08`

MATLAB vs ABAQUS: **PASS**.

- Peak roof displacement relative difference: `1.8320735134103193e-05`
- Maximum inter-story drift ratio relative difference: `2.059469853266933e-05`
- Peak isolation displacement relative difference: `7.894733270682764e-06`
- Peak total base shear relative difference: `2.0078727280391436e-06`
- Total hysteresis loop area relative difference: `6.434739527771795e-06`
- Roof displacement NRMSE: `2.5677850048343697e-05`
- Isolation displacement NRMSE: `2.3947855829177854e-05`
- Base shear NRMSE: `3.880533629355637e-05`

OpenSeesPy vs ABAQUS: **PASS**.

- Peak roof displacement relative difference: `1.833150970861689e-05`
- Peak isolation displacement relative difference: `7.909026975478886e-06`
- Peak total base shear relative difference: `2.007663812178544e-06`
- Roof displacement NRMSE: `2.567751635223817e-05`

## Final Statement

The final `cross_validation_report.json` has `final_result: "PASS"` and `strict_three_software_pass: true`. MATLAB, OpenSeesPy, and ABAQUS agree within the strict thresholds in this Example 6 run.
