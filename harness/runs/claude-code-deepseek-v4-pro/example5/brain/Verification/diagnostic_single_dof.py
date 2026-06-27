"""
Single-DOF diagnostic: compare OpenSees Steel01 vs exact bilinear (MATLAB algorithm)
for the isolation story (story 1) parameters under cyclic loading.
"""
import numpy as np

# Story 1 parameters
k0 = 50e6      # N/m
fy = 500e3     # N
b = 0.05       # hardening ratio
k_sh = b * k0  # post-yield stiffness = 2.5e6 N/m
fy_eff = fy * (1 - b)  # effective yield half-width = 475000 N

def matlab_return_mapping(last_shear, last_disp, disp, k0, fy, b):
    """Replicate MATLAB computeRsStiffness for a single story."""
    k_sh = b * k0
    fy_minus_b = fy * (1 - b)
    delta_disp = disp - last_disp
    c = last_shear + k0 * delta_disp  # elastic predictor
    shear = max(k_sh * disp - fy_minus_b,
                min(k_sh * disp + fy_minus_b, c))
    if abs(shear - c) < 1e-3:
        stiffness = k0
    else:
        stiffness = k_sh
    return shear, stiffness

def compute_monotonic_response(disp_history):
    """Compute MATLAB-style response for a displacement history."""
    shear = np.zeros(len(disp_history))
    stiff = np.zeros(len(disp_history))
    last_shear = 0.0
    last_disp = 0.0
    for i, d in enumerate(disp_history):
        shear[i], stiff[i] = matlab_return_mapping(last_shear, last_disp, d, k0, fy, b)
        last_shear = shear[i]
        last_disp = d
    return shear, stiff

# Create a cyclic displacement history
t = np.linspace(0, 10, 1000)
dt = t[1] - t[0]
freq = 1.0  # Hz
amplitude = 0.03  # m (large enough to yield)

# Sine wave with increasing amplitude
disp = amplitude * np.sin(2 * np.pi * freq * t) * (t / t[-1])

shear_m, stiff_m = compute_monotonic_response(disp)

# Save MATLAB-style result
np.savetxt('diagnostic_matlab_hyst.txt',
           np.column_stack([t, disp, shear_m, stiff_m]),
           fmt='%.10e',
           header='time disp shear stiffness')

print("=== MATLAB hysteresis diagnostic ===")
print(f"Parameters: k0={k0:.1f}, fy={fy:.1f}, b={b}")
print(f"Post-yield k: {k_sh:.1f}")
print(f"Max disp: {np.max(np.abs(disp)):.6f} m")
print(f"Max shear: {np.max(np.abs(shear_m)):.2f} N")
print(f"Yield disp: {fy/k0:.6f} m")

# Count yield events
n_yield = np.sum(stiff_m < k0 * 0.99)
print(f"Yield steps: {n_yield}/{len(disp)}")

# Now also save a simple 3-point test to verify the algorithm
print("\n=== 3-point verification ===")
test_disps = [0.0, 0.005, 0.02]  # elastic, near-yield, post-yield
last_s = 0.0
last_d = 0.0
for d in test_disps:
    s, st = matlab_return_mapping(last_s, last_d, d, k0, fy, b)
    print(f"  disp={d:.6f}: shear={s:.2f}, stiff={st:.1f}")
    last_s = s
    last_d = d
