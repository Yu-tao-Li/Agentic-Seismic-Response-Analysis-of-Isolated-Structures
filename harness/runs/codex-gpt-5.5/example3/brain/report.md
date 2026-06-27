# example3 codex-gpt-5.5 cross-software report

## Conclusion

- Overall verdict: `PASS`
- Overall pass: `true`
- Summary: The three software paths reached the adopted agreement criteria.

## Software files

| Software | Status | Main artifact |
| --- | --- | --- |
| matlab | completed | `verification/matlab/main.m` |
| openseespy | completed | `verification/openseespy/main.py` |
| abaqus | completed | `verification/abaqus/model.inp` |

## Thresholds

```json
{
  "peak_displacement_relative_error": 0.02,
  "peak_acceleration_relative_error": 0.03,
  "curve_nrmse": 0.03,
  "rationale": "Linear systems with matched matrices and input should agree closely. Example 3 allows slightly looser acceleration/curve error because isolation damping is represented by discrete dashpot elements in OpenSeesPy/ABAQUS and matrix damping in MATLAB."
}
```

## Comparison Summary

### matlab_vs_openseespy
- pass: `true`
#### displacement
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 0.142062 |
| comparison_peak | 0.142076 |
| peak_relative_error | 9.799632e-05 |
| nrmse | 4.176801e-05 |
| relative_l2_error | 0.000453179 |
#### acceleration
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 7.91652 |
| comparison_peak | 7.91684 |
| peak_relative_error | 4.061461e-05 |
| nrmse | 8.675205e-05 |
| relative_l2_error | 0.00118139 |

### matlab_vs_abaqus
- pass: `true`
#### displacement
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 0.142062 |
| comparison_peak | 0.142014 |
| peak_relative_error | 0.000335082 |
| nrmse | 0.000168774 |
| relative_l2_error | 0.00183118 |
#### acceleration
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 7.91652 |
| comparison_peak | 7.89001 |
| peak_relative_error | 0.00334895 |
| nrmse | 0.00116758 |
| relative_l2_error | 0.0159001 |

### openseespy_vs_abaqus
- pass: `true`
#### displacement
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 0.142076 |
| comparison_peak | 0.142014 |
| peak_relative_error | 0.000433036 |
| nrmse | 0.000180246 |
| relative_l2_error | 0.00195571 |
#### acceleration
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 7.91684 |
| comparison_peak | 7.89001 |
| peak_relative_error | 0.00338943 |
| nrmse | 0.00117313 |
| relative_l2_error | 0.0159753 |

## Notes

- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
