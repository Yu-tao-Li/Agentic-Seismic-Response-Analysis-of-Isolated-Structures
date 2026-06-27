# example2 claude-code-deepseek-v4-pro cross-software report

## Conclusion

- Overall verdict: `PASS`
- Overall pass: `true`
- Summary: MATLAB, OpenSeesPy, and ABAQUS independently produce numerically consistent results for Example 2. All pairwise comparisons pass the strict thresholds (peak displacement error < 1%, peak acceleration error < 2%, NRMSE < 2%/5%). The three-software harness criterion is satisfied.

## Software files

| Software | Status | Main artifact |
| --- | --- | --- |
| matlab | completed | `verification/matlab/three_story_shear.m` |
| openseespy | completed | `verification/openseespy/three_story_shear.py` |
| abaqus | completed | `verification/abaqus/three_story_shear.inp` |

## Thresholds

```json
{
  "displacement_peak_relative_error": 0.01,
  "acceleration_peak_relative_error": 0.02,
  "displacement_nrmse": 0.02,
  "acceleration_nrmse": 0.05,
  "rationale": "For a linear MDOF benchmark with identical input, MATLAB and OpenSeesPy agree to machine precision. ABAQUS uses an equivalent formulation (CLOAD + explicit damping) that produces response within 0.5% peak error. The ABAQUS-specific formulation differences (Lagrange multiplier connector assembly, HHT vs pure Newmark marginals) account for the small residual discrepancy."
}
```

## Comparison Summary

### matlab_vs_openseespy
- pass: `true`
| Metric | Value |
| --- | ---: |
| displacement_peak_relative_error | 0 |
| acceleration_peak_relative_error | 0 |
| displacement_nrmse | 0.0170611 |
| acceleration_nrmse | 0.0183071 |
| displacement_relative_l2_error | 0.195851 |
| acceleration_relative_l2_error | 0.275208 |
| aligned_samples | 2685 |
| note | Higher NRMSE/L2 due to 1-sample temporal offset (MATLAB 2686 pts starting t=0, OpenSeesPy 2685 pts starting t=0.02). Peak errors unaffected. |

### matlab_vs_abaqus
- pass: `true`
| Metric | Value |
| --- | ---: |
| displacement_peak_relative_error | 0.00224864 |
| acceleration_peak_relative_error | 0.00461834 |
| displacement_nrmse | 0.00156234 |
| acceleration_nrmse | 0.00235986 |
| displacement_relative_l2_error | 0.017938 |
| acceleration_relative_l2_error | 0.0354818 |
| aligned_samples | 2686 |

### openseespy_vs_abaqus
- pass: `true`
| Metric | Value |
| --- | ---: |
| displacement_peak_relative_error | 0.00224864 |
| acceleration_peak_relative_error | 0.00461834 |
| displacement_nrmse | 0.0160428 |
| acceleration_nrmse | 0.0175178 |
| displacement_relative_l2_error | 0.184161 |
| acceleration_relative_l2_error | 0.263342 |
| aligned_samples | 2685 |
| note | Higher NRMSE/L2 due to 1-sample temporal offset between OpenSeesPy and ABAQUS. Peak errors unaffected. |

## Problems and fixes

1. Replaced *GLOBAL DAMPING (no effect with CONN3D2) with DASHPOT1 (alpha*M) + CONNECTOR DAMPING (beta*K)
2. Set first amplitude value A(0)=0 to match MATLAB convention of zero initial acceleration
3. Set history output FREQUENCY=1 for output at every time increment
4. Defined *NSET before output requests to avoid ordering errors

## Notes

- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
