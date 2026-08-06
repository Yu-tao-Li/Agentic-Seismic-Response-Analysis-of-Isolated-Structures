# Dimensionless Tolerance Sensitivity

The adopted force-residual, displacement-correction, and state-transition tolerances are `1e-7`, `1e-10`, and `1e-10`, respectively. A full factorial sweep of 27 combinations varied each criterion over the values listed below.

- Force residual tolerance: `1e-6`, `1e-7`, `1e-8`
- Displacement correction tolerance: `1e-8`, `1e-10`, `1e-12`
- State-transition tolerance: `1e-8`, `1e-10`, `1e-12`

Maximum changes relative to the adopted case across all 27 combinations:

| Quantity | Maximum relative change or NRMSE |
| --- | ---: |
| Peak top displacement | 0 |
| Peak top acceleration | 0 |
| Peak isolation displacement | 0 |
| Peak base shear | 0 |
| Hysteresis loop area | 0 |
| Top displacement curve NRMSE | 0 |
| Top acceleration curve NRMSE | 0 |
| Isolation displacement curve NRMSE | 0 |
| Base shear curve NRMSE | 0 |

All runs had at most 0 failed steps and required at most 3 Newton iterations.
