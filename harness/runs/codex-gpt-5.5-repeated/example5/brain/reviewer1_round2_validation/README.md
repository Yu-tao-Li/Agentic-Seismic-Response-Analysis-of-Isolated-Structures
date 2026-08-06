# Reviewer 1 Round-2 Validation

This directory is the canonical corrected Example 5 verification path prepared for the second revision.

## Numerical criteria

- Newton residual: `norm(residual,1) / max(norm(deltaPeff,1), norm(restoringForce,1), norm(fy,1)) <= 1e-7`.
- Displacement correction: `norm(deltaU,inf) / max(norm(uTrial,inf), max(fy./k)) <= 1e-10`.
- Elastic/post-yield state transition: `abs(shear-elasticTrial) / max(abs(shear), abs(elasticTrial), fy) <= 1e-10`.
- Seismic inertia: `-M*r*ag`, with the full mass matrix and influence vector.

The three criteria are dimensionless. `M = diag(m)` constructs the diagonal mass matrix; it is not used as a replacement for the matrix-vector inertial product.

## Reproduction

From this directory in MATLAB R2025b or a compatible MATLAB release:

```matlab
run_tolerance_sensitivity
```

The command regenerates the corrected baseline outputs and a 27-run full-factorial tolerance sweep. The generated reports are:

- `response_summary.json`
- `tolerance_sensitivity_results.csv`
- `tolerance_sensitivity_summary.json`
- `tolerance_sensitivity_report.md`

## External reference

`compare_cui_reference.py` compares the corrected MATLAB outputs with the deposited Cui et al. reference outputs outside the coding-agent run archive. The comparison uses the common 0.00--18.99 s window without sign alignment, amplitude rescaling, or fitted time shifting.

Reference source:

J. Cui, X. Shen, and M. Yang, *Structural Seismic Response Analysis Programming and Application*, China Architecture & Building Press, Beijing, 2022 (in Chinese).

The exact reference file paths and SHA-256 hashes are recorded in `cui_external_reference_comparison.json`.

The historical generated run directories remain unchanged for auditability. This directory is the single canonical path for the corrected criteria cited in the second-revision response.
