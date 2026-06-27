# brain example4 cross-software report

## Conclusion

- Overall verdict: `PASS`
- Overall pass: `true`
- Summary: MATLAB, OpenSeesPy, and ABAQUS agree within 0.005% on peak response metrics. All three software paths produce auditable, numerically consistent results.

## Root Cause of Original OpenSeesPy Discrepancy

This OpenSeesPy version's `rayleigh()` command does NOT apply stiffness-proportional Rayleigh damping (`beta_k * K`) to `twoNodeLink` elements with uniaxial materials. Only mass-proportional Rayleigh (`alpha_m * M`) was active, causing a 5-8% under-damping and correspondingly larger peak responses.

**Fix:** Viscous uniaxialMaterial elements (alpha=1.0) placed in parallel with each spring element: `C_visc_i = beta_k * k_i`. Mass-proportional Rayleigh retained via `ops.rayleigh(alpha_m, 0.0, 0.0, 0.0)`.

## Comparison Summary

### matlab_vs_abaqus
| Metric | Value |
| --- | ---: |
| disp_error_pct | 0.003 |
| acc_error_pct | 0.0005 |
| iso_disp_error_pct | 0 |
| verdict | PASS |

### matlab_vs_openseespy
| Metric | Value |
| --- | ---: |
| disp_error_pct | 0.000 |
| acc_error_pct | 0.0015 |
| iso_disp_error_pct | 0 |
| verdict | PASS - after Rayleigh damping fix |

### abaqus_vs_openseespy
| Metric | Value |
| --- | ---: |
| disp_error_pct | 0.003 |
| acc_error_pct | 0.0010 |
| iso_disp_error_pct | 0 |
| verdict | PASS |

### matlab_vs_reference
| Metric | Value |
| --- | ---: |
| disp_error_pct | 0.001 |
| acc_error_pct | 0.0005 |
| verdict | PASS |

### abaqus_vs_reference
| Metric | Value |
| --- | ---: |
| disp_error_pct | 0.002 |
| acc_error_pct | 0 |
| verdict | PASS |

## Diagnostics

1. `undamped_compare.py` - structural model (M, K) confirmed identical (< 0.02%)
2. `damped_compare.py` - isolated difference to Rayleigh damping (16.6% offset in Rayleigh-only test)
3. `diagnose_component.py` - bug isolated: stiffness-proportional Rayleigh broken for twoNodeLink/zeroLength; mass-proportional works
4. `diagnose_viscous.py` - confirmed Viscous elements produce exact damping match with MATLAB
5. `verify_4dof_fix.py` - fix validated (< 0.01% error on full 4-DOF linear system)

## Notes

- This report is generated from the corresponding `cross_validation_report.json`; the JSON remains the auditable source for exact values.
