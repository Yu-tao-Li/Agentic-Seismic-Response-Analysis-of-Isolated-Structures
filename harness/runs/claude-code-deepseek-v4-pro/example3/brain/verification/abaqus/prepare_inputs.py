"""Prepare ABAQUS auxiliary input files and compute model parameters."""
import numpy as np

# Read and normalize ground motion
gm = np.loadtxt('../../input_data/GM1.txt')
dt = 0.02
g = 9.81
target_pga = 0.40 * g
max_abs_gm = np.max(np.abs(gm))
scale = target_pga / max_abs_gm
ag = gm * scale
N = len(ag)
t = np.arange(N) * dt

print(f"GM1: N={N}, max_abs={max_abs_gm:.8f}, scale={scale:.6f}, target_PGA={target_pga:.4f}")

# Write ABAQUS amplitude file (time, value pairs)
# ABAQUS expects amplitude at each time point
# Use DEFINE= option to specify how to interpret
with open('GM_normalized.inp', 'w') as f:
    for i in range(N):
        f.write(f'{t[i]:.10e}, {ag[i]:.10e}\n')
print(f"Written GM_normalized.inp with {N} rows")

# Compute Rayleigh damping coefficients
M = np.diag([500., 1000., 1000., 1000.])
k_iso = 200e3
k_story = 500e3
K = np.array([
    [k_iso+k_story, -k_story,       0.,        0.],
    [-k_story,      2*k_story, -k_story,        0.],
    [0.,           -k_story,  2*k_story, -k_story],
    [0.,              0.,     -k_story,   k_story]
])

eigvals, _ = np.linalg.eig(np.linalg.solve(M, K))
omega = np.sqrt(np.sort(eigvals.real))
print(f"Natural frequencies (rad/s): {omega}")

omega1, omega2 = omega[0], omega[1]
xi = 0.05
alpha = 2*xi*omega1*omega2/(omega1+omega2)
beta_ray = 2*xi/(omega1+omega2)
print(f"Rayleigh: alpha={alpha:.8f}, beta={beta_ray:.8e}")
print(f"Use in ABAQUS: *DAMPING, ALPHA={alpha:.6f}, BETA={beta_ray:.6e}")

# Write parameters file
with open('model_params.txt', 'w') as f:
    f.write(f"alpha={alpha:.10f}\n")
    f.write(f"beta={beta_ray:.10e}\n")
    f.write(f"scale={scale:.10f}\n")
    f.write(f"target_pga={target_pga:.10f}\n")
    f.write(f"dt={dt}\n")
    f.write(f"n_samples={N}\n")
    f.write(f"omega={omega.tolist()}\n")
