# Seismic Isolation Cross-Validation Report
## Model and Input
- Four dynamic DOFs: isolation mass plus stories 2, 3, and 4.
- Masses: `[500, 1000, 1000, 1000]` kg.
- Stiffnesses: `[200e3, 500e3, 500e3, 500e3]` N/m.
- Isolation damping: `2000` N*s/m added to the isolation damping term.
- Rayleigh damping ratio: `0.05`, fitted to modes 1 and 2.
- `GM1.txt` was normalized to `0.40 g = 3.92266 m/s^2`; `dt = 0.02 s`; samples = `2686`.
## Solver Implementations
- MATLAB: direct MDOF equation of motion with average-acceleration Newmark integration; isolation dashpot added after Rayleigh damping.
- OpenSeesPy: 1D shear-building zeroLength elastic springs and viscous dashpots; Rayleigh matrix realized as `beta*K` story dashpots plus `alpha*M` ground dashpots and the isolation dashpot.
- ABAQUS: equivalent `MASS`, `SPRING2`, and `DASHPOT2` model; direct dynamic step with HHT `alpha=0`; ground motion applied as equivalent inertial nodal loads and ODB post-processed to absolute top acceleration.
## Generated Outputs Used
- matlab: verification/matlab/run_matlab.m, verification/matlab/execution.log, verification/matlab/response.json, verification/matlab/topStoDis.txt, verification/matlab/topStoAcc.txt, verification/matlab/isolationDisplacement.txt
- openseespy: verification/openseespy/run_openseespy.py, verification/openseespy/execution.log, verification/openseespy/response.json, verification/openseespy/topStoDis.txt, verification/openseespy/topStoAcc.txt, verification/openseespy/isolationDisplacement.txt
- abaqus: verification/abaqus/generate_abaqus_input.py, verification/abaqus/isolated_building.inp, verification/abaqus/execution.log, verification/abaqus/postprocess_abaqus.py, verification/abaqus/postprocess.log, verification/abaqus/response.json, verification/abaqus/topStoDis.txt, verification/abaqus/topStoAcc.txt, verification/abaqus/isolationDisplacement.txt
## Thresholds
The models are linear and use the same SI matrices, normalized input, and average-acceleration Newmark family integration. Exact input consistency is therefore required. A 1% response tolerance is strict for independent solver implementations while allowing small differences in OpenSeesPy and ABAQUS solver internals, time-history bookkeeping, and output extraction. Adopted limits: peak relative error <= 1%, normalized L2 <= 1%, NRMSE <= 1%, trend correlation >= 0.999, and exact input consistency within numerical roundoff.
## Final Result
Strict three-software pass: **PASS**.
MATLAB, OpenSeesPy, and ABAQUS all generated nonempty response files in this run. The final comparison used only the three generated `response.json` files under `verification/`.
## Peak Responses
| Software | Top disp. (m) | Top abs. acc. (m/s^2) | Isolation disp. (m) |
| --- | ---: | ---: | ---: |
| matlab | 0.142013250702 | 5.71542917762 | 0.0798328024929 |
| openseespy | 0.142027167479 | 5.71586557443 | 0.0798413448295 |
| abaqus | 0.142013249688 | 5.71542666543 | 0.0798327958859 |
## Pairwise Metrics
| Pair | Response | Peak err. | Norm. L2 | NRMSE | Corr. |
| --- | --- | ---: | ---: | ---: | ---: |
| matlab_vs_openseespy | top_story_displacement | 9.79867e-05 | 0.000453179 | 4.1768e-05 | 0.999999905535 |
| matlab_vs_openseespy | top_story_acceleration | 7.63483e-05 | 0.000863376 | 8.02806e-05 | 0.999999634437 |
| matlab_vs_openseespy | isolation_layer_displacement | 0.000106991 | 0.000455471 | 4.14467e-05 | 0.999999904775 |
| matlab_vs_abaqus | top_story_displacement | 7.14285e-09 | 1.80524e-06 | 1.66383e-07 | 0.999999999998 |
| matlab_vs_abaqus | top_story_acceleration | 4.39545e-07 | 4.54626e-06 | 4.22732e-07 | 0.99999999999 |
| matlab_vs_abaqus | isolation_layer_displacement | 8.27605e-08 | 1.80457e-06 | 1.64212e-07 | 0.999999999998 |
| openseespy_vs_abaqus | top_story_displacement | 9.79939e-05 | 0.000453133 | 4.17627e-05 | 0.999999905537 |
| openseespy_vs_abaqus | top_story_acceleration | 7.67878e-05 | 0.000863311 | 8.02753e-05 | 0.999999634422 |
| openseespy_vs_abaqus | isolation_layer_displacement | 0.000107074 | 0.000455422 | 4.14411e-05 | 0.999999904778 |
## Errors Encountered and Corrections
- MATLAB: plain `matlab -batch` printed the test message but crashed during shutdown in this environment; rerun with `matlab -nojvm -batch`, which completed and wrote `execution.log`.
- OpenSeesPy: the active Python 3.13 environment had an OpenSeesPy Windows DLL import failure; rerun in the installed `opensees` conda environment with Python 3.12, which imported and executed successfully.
- ABAQUS: a syntax smoke test confirmed `SPRING2`/`DASHPOT2`. The first full run stopped at the default 100-increment limit; the input generator was corrected to set `inc=2696`, and the rerun completed 2685 increments with zero analysis errors.
## Final Assessment
All required generated outputs are present and nonempty. The input records are identical across solvers, the response trends are essentially perfectly correlated, and every pairwise response metric satisfies the adopted thresholds.
