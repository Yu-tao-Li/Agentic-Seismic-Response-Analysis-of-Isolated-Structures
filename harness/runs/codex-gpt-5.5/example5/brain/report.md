# example5 codex-gpt-5.5 cross-software report

## Conclusion

- Overall verdict: `PASS`
- Overall pass: `true`
- Summary: MATLAB, OpenSeesPy, and ABAQUS reach tight nonlinear-response agreement for Example 5 after matching the Steel01 bilinear kinematic rule and representing alpha*M plus beta*K0 damping explicitly in ABAQUS.

## Steel01 implementation check

- OpenSees Steel01 uses a bilinear kinematic-hardening rule with `b` as the post-yield tangent ratio.
- The MATLAB `computeRsStiffness` clamp rule matches the Steel01 trial-stress update, which is why MATLAB and OpenSeesPy agree nearly exactly.
- The ABAQUS mismatch came from the damping representation, not from the Steel01 rule itself: `alpha*M` was not active in the point-mass/truss model when represented only by global damping.

## Software files

| Software | Status | Main artifact |
| --- | --- | --- |
| matlab | completed | `verification/matlab/nonlinear_mdof_analysis.m` |
| openseespy | completed | `verification/openseespy/shear_building_analysis.py` |
| abaqus | completed | `verification/abaqus/example5_model.inp` |

## Final responses

| Quantity | MATLAB | OpenSeesPy | ABAQUS |
| --- | ---: | ---: | ---: |
| top_displacement (m) | 0.0281488 | 0.0281488 | 0.0281502 |
| top_acceleration (m/s^2) | 5.07916 | 5.07924 | 5.07914 |
| isolation_displacement (m) | 0.0231477 | 0.0231476 | 0.0231494 |
| base_shear (N) | 532869 | 532869 | 532874 |

## Thresholds

```json
{
  "peak_response_error_percent": 5.0,
  "curve_nrmse": 0.05,
  "convergence_rate": ">=95%",
  "rationale": "Steel01 is a bilinear kinematic-hardening rule. After representing both alpha*M and beta*K0 damping explicitly in ABAQUS, the nonlinear response should satisfy a stricter 5% peak-response and 5% curve-NRMSE agreement criterion."
}
```

## Comparison Summary

### matlab_vs_openseespy
- pass: `true`
#### top_displacement
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 0.028148812 |
| comparison_peak | 0.028148759 |
| peak_error_percent | 0.00018984093 |
| nrmse | 4.2243868e-05 |
| relative_l2_error | 0.00019383506 |
| sign_applied_to_comparison | 1 |
#### top_acceleration
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 5.0791572 |
| comparison_peak | 5.0792449 |
| peak_error_percent | 0.0017262927 |
| nrmse | 3.2787939e-05 |
| relative_l2_error | 0.00028569075 |
| sign_applied_to_comparison | 1 |
#### isolation_displacement
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 0.023147669 |
| comparison_peak | 0.023147619 |
| peak_error_percent | 0.00021293555 |
| nrmse | 3.768239e-05 |
| relative_l2_error | 0.0001820314 |
| sign_applied_to_comparison | 1 |
#### base_shear
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 532869.17 |
| comparison_peak | 532869.05 |
| peak_error_percent | 2.3125572e-05 |
| nrmse | 6.8401097e-05 |
| relative_l2_error | 0.00022414163 |
| sign_applied_to_comparison | -1 |

### matlab_vs_abaqus
- pass: `true`
#### top_displacement
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 0.028148812 |
| comparison_peak | 0.028150175 |
| peak_error_percent | 0.0048396157 |
| nrmse | 2.6535754e-05 |
| relative_l2_error | 0.00012175872 |
| sign_applied_to_comparison | 1 |
#### top_acceleration
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 5.0791572 |
| comparison_peak | 5.0791136 |
| peak_error_percent | 0.00085947731 |
| nrmse | 5.496209e-06 |
| relative_l2_error | 4.7890051e-05 |
| sign_applied_to_comparison | 1 |
#### isolation_displacement
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 0.023147669 |
| comparison_peak | 0.023149391 |
| peak_error_percent | 0.0074394426 |
| nrmse | 3.2513653e-05 |
| relative_l2_error | 0.00015706291 |
| sign_applied_to_comparison | 1 |
#### base_shear
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 532869.17 |
| comparison_peak | 532873.46 |
| peak_error_percent | 0.00080425735 |
| nrmse | 7.4740087e-06 |
| relative_l2_error | 2.4491369e-05 |
| sign_applied_to_comparison | 1 |

### openseespy_vs_abaqus
- pass: `true`
#### top_displacement
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 0.028148759 |
| comparison_peak | 0.028150175 |
| peak_error_percent | 0.0050294662 |
| nrmse | 5.0142568e-05 |
| relative_l2_error | 0.00023007735 |
| sign_applied_to_comparison | 1 |
#### top_acceleration
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 5.0792449 |
| comparison_peak | 5.0791136 |
| peak_error_percent | 0.0025857254 |
| nrmse | 3.3266172e-05 |
| relative_l2_error | 0.00028985911 |
| sign_applied_to_comparison | 1 |
#### isolation_displacement
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 0.023147619 |
| comparison_peak | 0.023149391 |
| peak_error_percent | 0.0076523945 |
| nrmse | 5.0144166e-05 |
| relative_l2_error | 0.00024223012 |
| sign_applied_to_comparison | 1 |
#### base_shear
| Metric | Value |
| --- | ---: |
| samples | 1899 |
| reference_peak | 532869.05 |
| comparison_peak | 532873.46 |
| peak_error_percent | 0.00082738311 |
| nrmse | 6.8844928e-05 |
| relative_l2_error | 0.00022559426 |
| sign_applied_to_comparison | -1 |

## Problems and fixes

1. MATLAB and OpenSeesPy already agreed because the MATLAB clamp update reproduces the OpenSees Steel01 bilinear kinematic-hardening rule.
2. The previous ABAQUS truss-plasticity model used explicit beta*K0 dashpots but relied on *GLOBAL DAMPING for alpha*M; changing the alpha value did not affect the results, showing that alpha*M was not represented consistently for this point-mass/truss model.
3. ABAQUS was corrected by adding explicit DASHPOT2 elements from each floor node to the fixed base with c_i = alpha*m_i, while retaining story-to-story DASHPOT2 elements with c_i = beta*k_i for beta*K0.
4. After the explicit alpha*M correction, MATLAB, OpenSeesPy, and ABAQUS satisfy the stricter 5% peak-response and 5% curve-NRMSE criteria for all tracked response quantities.

## Notes

- `base_shear` is the restoring shear extracted from the first nonlinear spring; `base_reaction` is kept separately because it includes damping forces.
- Base-shear sign conventions can differ by software; the JSON records sign alignment used for pairwise curve-error calculation.
- Trial ABAQUS variants are preserved with suffixes such as `_beta150`, `_beta165`, `_beta200`, and `_explicit_alpha_beta100`.
