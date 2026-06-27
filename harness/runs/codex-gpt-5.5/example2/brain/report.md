# example2 codex-gpt-5.5 cross-software report

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
  "peak_displacement_relative_error": 0.01,
  "peak_acceleration_relative_error": 0.02,
  "curve_nrmse": 0.02,
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
| reference_peak | 0.109024 |
| comparison_peak | 0.109022 |
| peak_relative_error | 1.910570e-05 |
| nrmse | 3.358954e-05 |
| relative_l2_error | 0.000385656 |
#### acceleration
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 13.7288 |
| comparison_peak | 13.7289 |
| peak_relative_error | 7.380131e-06 |
| nrmse | 5.408914e-05 |
| relative_l2_error | 0.00081322 |

### matlab_vs_abaqus
- pass: `true`
#### displacement
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 0.109024 |
| comparison_peak | 0.108829 |
| peak_relative_error | 0.00178918 |
| nrmse | 0.000543005 |
| relative_l2_error | 0.00623447 |
#### acceleration
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 13.7288 |
| comparison_peak | 13.6889 |
| peak_relative_error | 0.00290606 |
| nrmse | 0.00145695 |
| relative_l2_error | 0.021905 |

### openseespy_vs_abaqus
- pass: `true`
#### displacement
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 0.109022 |
| comparison_peak | 0.108829 |
| peak_relative_error | 0.00177011 |
| nrmse | 0.000549598 |
| relative_l2_error | 0.00631022 |
#### acceleration
| Metric | Value |
| --- | ---: |
| samples | 2686 |
| reference_peak | 13.7289 |
| comparison_peak | 13.6889 |
| peak_relative_error | 0.00291342 |
| nrmse | 0.00146119 |
| relative_l2_error | 0.0219699 |

## Notes

- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
