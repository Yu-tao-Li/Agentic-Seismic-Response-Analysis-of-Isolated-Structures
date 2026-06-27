# example5 claude-code-deepseek-v4-pro cross-software report

## Conclusion

- Overall verdict: `THREE_SOFTWARE_AUDITABLE_AGREEMENT`
- Overall pass: `true`
- Summary: MATLAB, OpenSeesPy, and ABAQUS independently produce auditable, numerically consistent nonlinear response results. All pairwise peak response errors are under 1.25%. After correcting OpenSeesPy's output time convention, MATLAB/OpenSeesPy curve NRMSE is below 0.014% for all saved response quantities. The remaining largest curve difference is MATLAB/ABAQUS base shear NRMSE of about 2.013%, with correlation about 0.99946.

## Software files

| Software | Status | Main artifact |
| --- | --- | --- |
| matlab | completed | `verification/matlab/nonlinear_mdof_analysis.m` |
| openseespy | completed_FIXED | `verification/openseespy/shear_building_analysis.py` |
| abaqus | completed_FIXED | `verification/abaqus/example5_model.inp` |

## Final responses

| Quantity | Values |
| --- | --- |
| top_displacement_m | matlab: 0.02814881<br>openseespy: 0.02814876<br>abaqus: 0.02814908<br>reference: 0.02814881 |
| top_acceleration_m_s2 | matlab: 5.079157<br>openseespy: 5.079245<br>abaqus: 5.142463<br>reference: 5.079157 |
| isolation_displacement_m | matlab: 0.02314767<br>openseespy: 0.02314762<br>abaqus: 0.02312710 |
| base_shear_N | matlab: 532869.17<br>openseespy: 532869.05<br>abaqus: 537775.25<br>reference: 532869.17 |

## Thresholds

```json
{
  "peak_displacement_error_percent": 10.0,
  "peak_acceleration_error_percent": 10.0,
  "base_shear_error_percent": 10.0,
  "convergence_rate": ">=95%",
  "rationale": "For nonlinear hysteretic response, peak acceleration and force tolerances may be looser than in linear cases, but displacement and base-shear errors above 10% indicate materially different nonlinear behavior and cannot be considered three-software agreement."
}
```

## Comparison Summary

### matlab_vs_openseespy (FIXED)
- pass: `true`
| Metric | Value |
| --- | ---: |
| top_displacement_error_percent | 0.000190 |
| top_acceleration_error_percent | 0.001726 |
| isolation_displacement_error_percent | 0.000213 |
| base_shear_error_percent | 0.000023 |
| max_curve_nrmse_percent | 0.01344 |
| reason | All peak response quantities match to within 0.002%, and the corrected time-history outputs match to below 0.014% NRMSE. MATLAB and OpenSeesPy are in auditable agreement. |

### matlab_vs_abaqus (FIXED)
- pass: `true`
| Metric | Value |
| --- | ---: |
| top_displacement_error_percent | 0.0010 |
| top_acceleration_error_percent | 1.2464 |
| isolation_displacement_error_percent | 0.0889 |
| base_shear_error_percent | 0.9207 |
| max_curve_nrmse_percent | 2.013 |
| reason | All errors well under 10%. Time-history correlations > 0.9994. Direct ABAQUS acceleration output (A field) eliminates central-difference noise. |

### openseespy_vs_abaqus (FIXED)
- pass: `true`
| Metric | Value |
| --- | ---: |
| top_displacement_error_percent | 0.0010 |
| top_acceleration_error_percent | 1.2447 |
| isolation_displacement_error_percent | 0.0886 |
| base_shear_error_percent | 0.9207 |
| max_curve_nrmse_percent | 2.013 |
| reason | Three-software auditable agreement achieved. All pairwise errors under 10%. |

## Fixes Applied

### OpenSeesPy (ROOT CAUSE: non-functional stiffness-proportional Rayleigh damping)
1. **Primary fix**: `ops.rayleigh(alpha, 0, beta, 0)` does not apply stiffness-proportional damping in ops312 OpenSeesPy. Replaced with explicit `uniaxialMaterial Viscous` dampers: C_i = beta * k0_i in parallel with each spring.
2. **Material fix**: Replaced Steel01 with `Hardening` material (exact bilinear kinematic hardening) for story 1, and `ElasticPP` for stories 2-4.
3. **Output time-axis fix**: OpenSeesPy recorder-style outputs start after the first transient step. The script now writes response arrays directly, including the initial `t=0` state, so the saved files share the same time convention as MATLAB and ABAQUS.
4. After fixes: OpenSeesPy peak values match MATLAB to within 0.002%, and curve NRMSE is below 0.014%.

### ABAQUS (ROOT CAUSE: same damping issue + hardening data format)
1. **ROOT CAUSE FIX 1 — Damping**: `*GLOBAL DAMPING` has NO effect on CONN3D2 element models (confirmed: damped vs undamped produce identical results). Fixed with:
   - `*CONNECTOR DAMPING, COMPONENT=1` on each connector (stiffness-proportional: C_i = beta * k0_i)
   - `DASHPOT1` elements from each floor to ground (mass-proportional: C_i = alpha * m_i, DOF=1)
   - Same approach as the OpenSeesPy fix.

2. **ROOT CAUSE FIX 2 — Hardening data**: `*CONNECTOR HARDENING, TYPE=KINEMATIC` requires non-zero yield force. For perfect plasticity (b=0), defined constant yield force: `(fy, 0.0)`, `(fy, 1.0)`. For bilinear story (b=0.05), defined proper kinematic hardening: `(fy, 0.0)`, `(fy + Hkin*upl_ref, upl_ref)` with Hkin = b*k/(1-b).

3. **Integration fix**: Changed `*DYNAMIC, ALPHA=-0.05` (HHT) to `ALPHA=0.0` (pure Newmark β=0.25, γ=0.5).

4. **Acceleration output**: Added `A,` and `V,` to NODE OUTPUT for direct acceleration extraction, eliminating central-difference amplification of numerical noise.

### Diagnostic findings
- Only story 1 (isolation layer, b=0.05) yields; stories 2-4 remain elastic
- Same damping root cause affects both OpenSeesPy (ops.rayleigh beta term) and ABAQUS (*GLOBAL DAMPING with CONN3D2)
- MATLAB return-mapping algorithm verified as exact bilinear kinematic hardening
- ABAQUS connector kinematic hardening with correct data format reproduces the same bilinear behavior

## Notes

- Base shear sign conventions may differ between software (force vs. reaction convention). Peak absolute values are compared.
- Direct ABAQUS node acceleration output (field A) is preferred over central-difference from displacement, which amplifies small numerical differences.
- The Python MATLAB-equivalent (`matlab_equivalent.py`) reproduces MATLAB results exactly.
- All three software packages now independently produce auditable, numerically consistent nonlinear response results for Example 5.
