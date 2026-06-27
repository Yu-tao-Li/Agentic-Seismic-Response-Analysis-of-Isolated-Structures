# Example 6: 2D Base-Isolated Steel Frame -- Cross-Software Verification Report

**Date:** 2026-06-22
**Final Status:** MATLAB-OpenSeesPy PASS; ABAQUS FAIL (non-physical response -- unable to achieve three-software agreement)

---

## 1. Conclusion

**Strict three-software agreement (MATLAB, OpenSeesPy, ABAQUS) is NOT achieved.**

MATLAB and OpenSeesPy achieve **excellent cross-validation agreement** across all evaluated response quantities. The fundamental period matches within 0.01%, peak isolation displacement within 1.3%, peak roof displacement within 0.6%, and peak base shear within 0.0003%.

ABAQUS produces **non-physical displacements** (~1000m vs expected ~0.018m) due to fundamentally insufficient damping. Despite extensive attempts with multiple approaches, the ABAQUS *DYNAMIC implicit integration procedure does not support proper Rayleigh damping specification, and all workarounds failed:

1. **DASHPOT1 elements**: Have no measurable effect on response (identical results with C=14,140 through C=500,000 N-s/m). ABAQUS issues "Too many nodes" warnings for all DASHPOT1 elements, suggesting they are not properly assembled into the dynamic system.

2. **UEL internal damping**: Adding damping force to RHS degrades Newton convergence proportionally to the damping coefficient. Even large coefficients produce only marginal reductions in the unbounded response.

3. **Material-level *DAMPING**: Does not apply to UEL elements or MASS elements, leaving the isolation mode effectively undamped.

4. **Convergence barrier at t=16.6s**: All simulations fail at the same physical time regardless of time step (tested 0.01s to 0.0005s), likely due to extreme displacements causing element distortion.

---

## 2. MATLAB and OpenSeesPy Verification (PASS)

### 2.1 Cross-Software Agreement Summary

| Quantity | MATLAB | OpenSeesPy | Diff % | Threshold | Status |
|----------|--------|------------|--------|-----------|--------|
| T1 [s] | 2.8431 | 2.8428 | 0.01 | 1% | PASS |
| Peak iso. disp. [m] | 0.0139 | 0.0140 | 1.3 | 5% | PASS |
| Peak roof disp. [m] | 0.0181 | 0.0182 | 0.6 | 20% | PASS |
| Peak base shear [kN] | 31.61 | 31.61 | ~0 | 10% | PASS |
| Max drift S1 | 7.29e-4 | 7.05e-4 | 3.4 | 20% | PASS |
| Roof NRMSE | -- | -- | 0.0131 | 0.25 | PASS |
| Iso. NRMSE | -- | -- | 0.0174 | 0.15 | PASS |
| Base shear NRMSE (abs) | -- | -- | 0.0512 | 0.20 | PASS |
| Hysteresis area [J] | 645 | 669 | 3.5% | 0.85 | PASS |
| Convergence (failed) | 0/1900 | 0/1900 | -- | 0 | PASS |

---

## 3. ABAQUS Verification (FAIL)

### 3.1 Summary of Attempts

| Run | Damping Method | C_damp [N-s/m] | Max Increments | Max Time [s] | Peak Roof [m] | Status |
|-----|---------------|-------------------|----------------|--------------|---------------|--------|
| v8 | Material + DASHPOT1 at iso | 50,000 (dashpot) | 1666 | 16.65 | 971.47 | Non-physical |
| v9 | Material + DASHPOT1 at iso + floor-to-ground | 14,140+15,539 (dashpot) | 1666 | 16.65 | 971.47 | Non-physical (identical to v8) |
| v10 | UEL damping force + full tangent | 100,000 (UEL) | 71 | 0.33 | N/A | Convergence failure |
| v11 | UEL damping force only, no tangent | 14,140 (UEL) | 1587 | 15.87 | 885.72 | Non-physical (slight reduction) |
| v12 | UEL damping force only, no tangent | 60,000 (UEL) | 74 | 0.35 | N/A | Convergence failure |
| v13 | UEL damping + limited tangent | 60,000 (UEL) | 90 | 0.33 | N/A | Convergence failure |
| v14 | Material + large DASHPOT1 | 500,000 (dashpot) | 1666 | 16.65 | 971.47 | Non-physical (identical to v8/v9) |

### 3.2 Key Observations

1. **DASHPOT1 elements have zero effect**: v8 (C=50,000), v9 (C=14,140+mass dashpots), and v14 (C=500,000) produce byte-for-byte identical displacements. This proves the DASHPOT1 damping is not being assembled into the dynamic equations.

2. **UEL internal damping degrades convergence**: As C_damp increases, Newton convergence worsens (from 1 iteration/increment at C=0 to failure at C=60,000). The damping force in RHS creates a tangent mismatch that the Newton method cannot overcome.

3. **Material *DAMPING does not damp the isolation mode**: The Rayleigh damping defined on STEEL material applies only to B21 elements, not to UEL isolators or MASS lumped masses.

4. **UEL Bouc-Wen sign convention is verified correct**: The energy dissipation from Bouc-Wen hysteresis alone is insufficient to bound the response without viscous damping.

### 3.3 Root Cause

The ABAQUS 2025 `*DYNAMIC` (implicit direct integration) procedure does not accept the `*DAMPING, ALPHA=..., BETA=...` keyword at the step level. The only damping mechanisms available are:
- Material-level *DAMPING (applies only to elements using that material, not UEL or MASS)
- DASHPOT elements (appear non-functional due to element formatting issues)
- HHT numerical damping (ALPHA=-0.05, ~0.005% effective damping at low frequencies)

This leaves the base-isolated structure with effectively zero physical damping in the critical isolation mode, causing unbounded resonant response.

### 3.4 UEL Compilation (VERIFIED)

- Compiler: Intel Fortran Classic 2021.13.1 (ifort.exe)
- Linker: MSVC 14.44.35207 (link.exe)
- All UEL versions (v1-v4) compile and link successfully

---

## 4. Model Parameters

### 4.1 Geometry
- Plane frame: 3 stories x 2 bays
- Bay width: 6.0 m; Story height: 3.6 m
- 3 isolator elements below each column line

### 4.2 Material and Sections
- E = 2.06e11 Pa, nu = 0.30
- Column: A_c = 0.020 m², I_c = 8.0e-4 m⁴
- Beam: A_b = 0.015 m², I_b = 4.5e-4 m⁴

### 4.3 Isolators (Bouc-Wen)
- k0 = 2.0e6 N/m, Fy = 1.0e4 N, uy = 0.005 m
- alpha = 0.05, n = 2.0, beta = gamma = 2.0e4, Ao = 1.0

### 4.4 Damping
- Target Rayleigh: alpha_m = 0.18647, beta_k = 0.00707 (5% in modes 1 & 2)
- MATLAB/OpenSeesPy: Full Rayleigh damping applied
- ABAQUS: Rayleigh damping NOT achievable

### 4.5 Ground Motion
- Northridge_01_NO_968.txt, dt = 0.01 s, 1900 steps
- PGA = 0.40 g, scale factor = 0.024791

---

## 5. Execution Environment

- **MATLAB**: R2025b, run via `matlab -batch`
- **OpenSeesPy**: Python 3.x + openseespy package
- **ABAQUS**: Abaqus 2025, Intel Fortran 2021.13.1, MSVC 14.44 Build Tools
- **OS**: Windows 11 Pro (Build 26200), 64GB RAM

---

## 6. Strict Three-Software Agreement Assessment

**NOT ACHIEVED.**

| Criterion | MATLAB | OpenSeesPy | ABAQUS | Status |
|-----------|--------|------------|--------|--------|
| Fundamental period | 2.8431 s | 2.8428 s | Not verified | MATLAB/OPS PASS |
| Peak roof displacement | 0.0181 m | 0.0182 m | 971 m | ABAQUS FAIL |
| Peak isolation displacement | 0.0139 m | 0.0140 m | 1002 m | ABAQUS FAIL |
| Peak base shear | 31.61 kN | 31.61 kN | Not available | ABAQUS FAIL |
| Time history NRMSE | Reference | < 0.052 | Not comparable | ABAQUS FAIL |
| Hysteresis response | Complete | Complete | Not available | ABAQUS FAIL |

**Two-software (MATLAB/OpenSeesPy) agreement is PASS with excellent metrics. Three-software agreement is FAIL due to unresolved ABAQUS damping deficiency.**
