# brain example5 cross-software report

## Conclusion

- Overall verdict: `UNKNOWN`
- Overall pass: `null`
- Summary: The run does not satisfy the recorded checks.

## Final responses

| Quantity | Values |
| --- | --- |
| top_displacement_m | matlab: 0.0281488<br>openseespy: 0.0246977<br>abaqus: 0.0835<br>reference: 0.0281488<br>matlab_vs_openseespy_pct_diff: -12.26<br>matlab_vs_abaqus_pct_diff: 196.6 |
| top_acceleration_ms2 | matlab: 5.07916<br>openseespy: 5.17834<br>abaqus: 4.66289<br>reference: 5.07916<br>matlab_vs_openseespy_pct_diff: 1.95<br>matlab_vs_abaqus_pct_diff: -8.2 |
| isolation_displacement_m | matlab: 0.0231487<br>openseespy: 0.0188663<br>abaqus: 0.0391666<br>matlab_vs_openseespy_pct_diff: -18.5<br>matlab_vs_abaqus_pct_diff: 69.2 |
| base_shear_N | matlab: 5.328692e+05<br>openseespy: 5.221660e+05<br>abaqus: 1.958331e+06<br>matlab_vs_openseespy_pct_diff: -2.01<br>matlab_vs_abaqus_pct_diff: 267.5 |

## Notes

- Base shear sign conventions can differ by software; the JSON records any sign adjustment used for curve-error calculation.
- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
