# External Independent Reference Comparison

The corrected Example 5 MATLAB implementation was compared against the archived output from Cui, Shen, and Yang (2022). The Cui output predates and was generated outside the coding-agent harness. SHA-256 hashes identify the exact deposited reference files used in this comparison.

The comparison uses the common 0.00-18.99 s window (1900 samples at 0.01 s). No sign alignment, amplitude rescaling, or fitted time shift was applied.

| Metric | Value |
| --- | ---: |
| top displacement peak relative error | 4.30217e-12 |
| top displacement NRMSE | 3.82276e-12 |
| top displacement relative L2 error | 1.75385e-11 |
| top acceleration peak relative error | 3.93766e-12 |
| top acceleration NRMSE | 1.69595e-12 |
| top acceleration relative L2 error | 1.47803e-11 |
| hysteresis displacement NRMSE | 3.82384e-12 |
| hysteresis force NRMSE | 2.37002e-12 |
| hysteresis loop-area relative error | 2.09593e-12 |

Reference provenance and file hashes are recorded in `cui_external_reference_comparison.json`.
