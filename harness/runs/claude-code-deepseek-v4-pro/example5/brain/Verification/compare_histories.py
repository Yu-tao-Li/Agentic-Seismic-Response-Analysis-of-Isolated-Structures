"""
Compare MATLAB and OpenSees time histories to diagnose where divergence occurs.
"""
import numpy as np

# Load MATLAB results
matlab_top_disp = np.loadtxt('matlab/topStoDisIso2.txt')
matlab_iso_disp = np.loadtxt('matlab/iso_displacement.txt')
matlab_base_shear = np.loadtxt('matlab/base_shear.txt')
matlab_top_accel = np.loadtxt('matlab/topStoAccIso2.txt')
matlab_hyst = np.loadtxt('matlab/hysteresis_layer1.txt')

# Load OpenSees results
ops_top_disp = np.loadtxt('openseespy/topStoDisIso2.txt')
ops_iso_disp = np.loadtxt('openseespy/iso_displacement.txt')
ops_base_shear = np.loadtxt('openseespy/base_shear.txt')
ops_top_accel = np.loadtxt('openseespy/topStoAccIso2.txt')
ops_hyst = np.loadtxt('openseespy/hysteresis_layer1.txt')

# Align lengths
n = min(len(matlab_top_disp), len(ops_top_disp))
time = matlab_top_disp[:n, 0]

mt = matlab_top_disp[:n, 1]
ot = ops_top_disp[:n, 1]
mi = matlab_iso_disp[:n, 1]
oi = ops_iso_disp[:n, 1]
mbs = matlab_base_shear[:n, 1]
obs = ops_base_shear[:n, 1]
ma = matlab_top_accel[:n, 1]
oa = ops_top_accel[:n, 1]

# Compute differences
top_disp_diff = ot - mt
iso_disp_diff = oi - mi
bs_diff = obs - mbs
acc_diff = oa - ma

# Find where large differences start
print("=== Top Displacement Difference Analysis ===")
print(f"Max abs diff: {np.max(np.abs(top_disp_diff)):.6f} m")
print(f"RMS diff: {np.sqrt(np.mean(top_disp_diff**2)):.6f} m")
print(f"Max abs MATLAB: {np.max(np.abs(mt)):.6f} m")
print(f"Max abs OpenSees: {np.max(np.abs(ot)):.6f} m")

# Find time of max difference
idx_max = np.argmax(np.abs(top_disp_diff))
print(f"\nMax difference at t={time[idx_max]:.3f}s")
print(f"  MATLAB: {mt[idx_max]:.6f}, OpenSees: {ot[idx_max]:.6f}")
print(f"  Diff: {top_disp_diff[idx_max]:.6f}")

# Find where difference exceeds 1% of peak
threshold = 0.01 * np.max(np.abs(mt))
first_significant = np.where(np.abs(top_disp_diff) > threshold)[0]
if len(first_significant) > 0:
    idx_first = first_significant[0]
    print(f"\nFirst time diff > 1% of peak ({threshold:.6f} m):")
    print(f"  t={time[idx_first]:.3f}s")
    print(f"  MATLAB disp: {mt[idx_first]:.6f}, OpenSees: {ot[idx_first]:.6f}")

# Check hysteresis accumulation
print("\n=== Hysteresis Comparison ===")
mh_force = matlab_hyst[:n, 2]  # force
oh_force = ops_hyst[:n, 2]
mh_disp = matlab_hyst[:n, 1]   # disp
oh_disp = ops_hyst[:n, 1]

# Compute hysteresis area (energy dissipated)
# Cumulative sum of force * delta_disp
m_energy = np.cumsum(np.abs(mh_force[1:] * np.diff(mh_disp)))
o_energy = np.cumsum(np.abs(oh_force[1:] * np.diff(oh_disp)))
print(f"MATLAB total hysteretic energy: {m_energy[-1]:.2f} J")
print(f"OpenSees total hysteretic energy: {o_energy[-1]:.2f} J")
print(f"Energy ratio: {o_energy[-1]/m_energy[-1]:.4f}")

# Look at yield flags
matlab_yield = np.loadtxt('matlab/yield_flags.txt')
print(f"\n=== Yielding Summary ===")
print(f"MATLAB yielding steps:")
for i in range(4):
    n_yield_matlab = np.sum(matlab_yield[:n, i+1])
    print(f"  Story {i+1}: {n_yield_matlab} steps ({n_yield_matlab/n*100:.1f}%)")

# Check initial stiffness period
print("\n=== Initial Response (first 50 steps) ===")
idx_early = 50
print(f"Top disp diff at t={time[idx_early]:.2f}s: {top_disp_diff[idx_early]:.8f} m")
print(f"Iso disp diff at t={time[idx_early]:.2f}s: {iso_disp_diff[idx_early]:.8f} m")

# Look at the response growth pattern
# Segment the time history and compute RMS difference per segment
segment_size = 200
n_segments = n // segment_size
print("\n=== RMS Difference by Segment ===")
for s in range(n_segments):
    start = s * segment_size
    end = start + segment_size
    rms_top = np.sqrt(np.mean(top_disp_diff[start:end]**2))
    rms_iso = np.sqrt(np.mean(iso_disp_diff[start:end]**2))
    t_start = time[start]
    t_end = time[end-1]
    print(f"  t=[{t_start:.1f},{t_end:.1f}]: top_rms={rms_top:.6f}, iso_rms={rms_iso:.6f}")
