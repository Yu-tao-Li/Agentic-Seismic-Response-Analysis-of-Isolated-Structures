# example1 codex-gpt-5.5 cross-software report

## Conclusion

- Overall verdict: `PASS`
- Overall pass: `true`
- Summary: MATLAB, OpenSeesPy, and ABAQUS reached auditable agreement for Example 1 modal frequencies, periods, and mode shapes.

## Software files

| Software | Status | Main artifact |
| --- | --- | --- |
| matlab | completed | `verification/matlab/modal_results.json` |
| openseespy | completed | `verification/openseespy/modal_results.json` |
| abaqus | completed | `verification/abaqus/modal_results.json` |

## Metrics

| Metric | Value |
| --- | ---: |
| max_frequency_relative_error_against_reference | 1.304339e-05 |
| max_period_relative_error_against_reference | 1.303901e-05 |
| max_mode_shape_absolute_error_against_reference | 1.892258e-08 |
| max_pairwise_frequency_relative_difference | 1.304405e-05 |
| max_pairwise_period_relative_difference | 1.304388e-05 |
| max_pairwise_mode_shape_absolute_difference | 2.311774e-08 |

## Pairwise Summary

### MATLAB vs OpenSeesPy
| Metric | Value |
| --- | ---: |
| max_frequency_relative_difference | 5.355075e-16 |
| max_period_relative_difference | 5.275182e-16 |
| max_mode_shape_absolute_difference | 3.330669e-16 |

### MATLAB vs ABAQUS
| Metric | Value |
| --- | ---: |
| max_frequency_relative_difference | 1.304405e-05 |
| max_period_relative_difference | 1.304388e-05 |
| max_mode_shape_absolute_difference | 2.311774e-08 |

### OpenSeesPy vs ABAQUS
| Metric | Value |
| --- | ---: |
| max_frequency_relative_difference | 1.304405e-05 |
| max_period_relative_difference | 1.304388e-05 |
| max_mode_shape_absolute_difference | 2.311774e-08 |

## Notes

- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
