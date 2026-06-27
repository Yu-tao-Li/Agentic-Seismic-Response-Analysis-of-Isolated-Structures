# example4 codex-gpt-5.5 cross-software report

## Conclusion

- Overall verdict: `PASS`
- Overall pass: `true`
- Summary: MATLAB, OpenSeesPy, and ABAQUS now agree for Example 4 after fixing OpenSeesPy Rayleigh participation and output time alignment.

## Software files

| Software | Status | Main artifact |
| --- | --- | --- |
| matlab | completed | `verification/matlab/nonlinear_mdof_analysis.m` |
| openseespy | completed | `verification/openseespy/nonlinear_mdof_analysis.py` |
| abaqus | completed | `verification/abaqus/model.inp` |

## Final responses

| Quantity | Values |
| --- | --- |
| top_displacement | matlab: 0.0360226<br>openseespy: 0.0360226<br>abaqus: 0.0360223 |
| top_acceleration | matlab: 5.39531<br>openseespy: 5.39539<br>abaqus: 5.39534 |
| isolation_displacement | matlab: 0.0250453<br>openseespy: 0.0250453<br>abaqus: 0.0250451 |

## Thresholds

```json
{
  "peak_response_error_percent": 1.0,
  "curve_nrmse": 0.01,
  "rationale": "No yielding occurred in this 0.40 g run, so the response is effectively linear; after matching Rayleigh participation and time-vector convention, a 1% peak threshold and 1% curve NRMSE are conservative."
}
```

## Comparison Summary

### matlab_vs_openseespy
- pass: `true`
#### top_displacement
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 0.0360226 |
| comparison_peak | 0.0360226 |
| peak_error_percent | 5.509845e-06 |
| nrmse | 6.538451e-05 |
| relative_l2_error | 0.000366331 |
#### top_acceleration
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 5.39531 |
| comparison_peak | 5.39539 |
| peak_error_percent | 0.00148421 |
| nrmse | 3.647028e-05 |
| relative_l2_error | 0.000331499 |
#### isolation_displacement
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 0.0250453 |
| comparison_peak | 0.0250453 |
| peak_error_percent | 3.926699e-06 |
| nrmse | 6.651338e-05 |
| relative_l2_error | 0.000369919 |

### matlab_vs_abaqus
- pass: `true`
#### top_displacement
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 0.0360226 |
| comparison_peak | 0.0360223 |
| peak_error_percent | 0.000907141 |
| nrmse | 4.363379e-05 |
| relative_l2_error | 0.000244468 |
#### top_acceleration
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 5.39531 |
| comparison_peak | 5.39534 |
| peak_error_percent | 0.000481265 |
| nrmse | 1.437762e-05 |
| relative_l2_error | 0.000130686 |
#### isolation_displacement
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 0.0250453 |
| comparison_peak | 0.0250451 |
| peak_error_percent | 0.00092297 |
| nrmse | 4.438644e-05 |
| relative_l2_error | 0.000246859 |

### openseespy_vs_abaqus
- pass: `true`
#### top_displacement
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 0.0360226 |
| comparison_peak | 0.0360223 |
| peak_error_percent | 0.000912651 |
| nrmse | 2.192784e-05 |
| relative_l2_error | 0.000122855 |
#### top_acceleration
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 5.39539 |
| comparison_peak | 5.39534 |
| peak_error_percent | 0.00100293 |
| nrmse | 3.065295e-05 |
| relative_l2_error | 0.000278624 |
#### isolation_displacement
| Metric | Value |
| --- | ---: |
| samples | 1900 |
| reference_peak | 0.0250453 |
| comparison_peak | 0.0250451 |
| peak_error_percent | 0.000926897 |
| nrmse | 2.230653e-05 |
| relative_l2_error | 0.000124059 |

## Problems and fixes

1. Initial discrepancy was caused by OpenSeesPy elements not participating in Rayleigh damping; MATLAB and ABAQUS both included alpha*M + beta*K0.
2. After enabling -doRayleigh, peak responses matched but curve comparison still showed an artificial error because OpenSeesPy output started at 0.01 s.
3. Recording the initial state at t=0 removed the time shift and produced auditable three-software agreement.

## Figure note

- The saved `hysteresis_layer1.txt` files use mixed force conventions: MATLAB stores only the linear isolation spring force `keq*u_iso`, while OpenSeesPy and ABAQUS store reaction-type quantities that include damping contributions. Therefore, panel (d) in `example4_isolated_response_comparison.pdf` recomputes the same isolation-layer total force for all three software paths as `F_iso = keq*u_iso + c_iso*v_iso`, using `keq = 50e6 N/m` and `c_iso = 1000e3 N/(m/s)`.
- This plotting correction affects only the force-displacement visualization. The pass/fail validation in this report is still based on the directly compared top displacement, top acceleration, and isolation displacement time histories.

## Notes

- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
