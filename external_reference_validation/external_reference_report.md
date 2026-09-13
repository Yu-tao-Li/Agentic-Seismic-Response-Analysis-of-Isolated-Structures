# External reference report

This report consolidates existing author-confirmed results. No Chen program or implementation associated with Cui et al. was rerun or modified during this packaging step.

| Example | External reference | Responses checked | Peak/modal metric | Curve/shape metric | Result |
| --- | --- | --- | --- | --- | --- |
| 1 | Chen compiled program | Frequencies; periods; mode shapes | max frequency error=4.885e-10; max period error=1.575e-08 | min MAC=1.000000000 | Consistent |
| 2 | Chen compiled program | Top displacement | peak error=3.018e-04 | NRMSE=4.957e-05 | Consistent |
| 3 | Not available | - | - | - | Internal only |
| 4 | Cui et al. published implementation | Top displacement; top relative acceleration; isolation displacement | max peak error=1.498e-02 | max NRMSE=7.760e-03 | Consistent |
| 5 | Cui et al. published implementation | Top response; isolation displacement and force; hysteresis | max peak error=4.302e-12; loop-area error=2.096e-12 | max NRMSE=3.824e-12 | Consistent |
| 6 | Not available | - | - | - | Internal only |

## Interpretation notes

- Example 1 uses whole-mode sign alignment only because an eigenvector's global sign is arbitrary.
- Example 2 uses the common 0.00-39.98 s interval because the archived Chen output contains 2000 samples.
- The archived Cui et al. reference acceleration quantity for Example 4 is relative acceleration, so it is compared with the matching corrected harness MATLAB quantity.
- The archived Example 4 force series uses a different force-output convention and is retained but not summarized as a force discrepancy.
- Example 5 compares top response, isolation displacement and restoring force, and hysteresis-loop area over 0.00-18.99 s.
- No response amplitude scaling, fitted time shift, or response sign adjustment was applied.
