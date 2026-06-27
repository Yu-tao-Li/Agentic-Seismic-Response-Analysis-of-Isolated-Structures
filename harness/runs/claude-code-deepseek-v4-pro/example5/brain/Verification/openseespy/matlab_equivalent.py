"""
Python implementation of the exact MATLAB nonlinear MDOF algorithm.
Used to verify OpenSees results and diagnose discrepancies.
"""
import numpy as np

# Parameters
m = np.array([250, 270, 270, 180]) * 1e3
k0 = np.array([50, 245, 195, 98]) * 1e6
fy = np.array([500, 1225, 975, 490]) * 1e3
b_hard = np.array([0.05, 0.0, 0.0, 0.0])
n_stories = 4

dt = 0.01
g = 9.81
PGA = 0.40 * g

# Load ground motion
data_path = '../../input_data/Northridge_01_NO_968.txt'
acc_raw = np.loadtxt(data_path)
acc_data = acc_raw / np.max(np.abs(acc_raw)) * PGA
n_steps = len(acc_data)
time = np.arange(n_steps) * dt

# Mass matrix
M = np.diag(m)

# Initial stiffness matrix
K0 = np.zeros((4, 4))
K0[0, 0] = k0[0] + k0[1]
K0[0, 1] = -k0[1]
K0[1, 0] = -k0[1]
K0[1, 1] = k0[1] + k0[2]
K0[1, 2] = -k0[2]
K0[2, 1] = -k0[2]
K0[2, 2] = k0[2] + k0[3]
K0[2, 3] = -k0[3]
K0[3, 2] = -k0[3]
K0[3, 3] = k0[3]

# Rayleigh damping
damping_ratio = 0.05
eigvals, _ = np.linalg.eig(np.linalg.solve(M, K0))
omega = np.sqrt(np.sort(np.real(eigvals)))
w1, w2 = omega[0], omega[1]
alpha = 2 * damping_ratio * w1 * w2 / (w1 + w2)
beta_val = 2 * damping_ratio / (w1 + w2)
C = alpha * M + beta_val * K0
print(f"alpha={alpha:.6f}, beta={beta_val:.6f}")

# Newmark parameters
beta_n = 0.25
gamma = 0.5

def compute_stiffness_matrix(story_stiffness):
    """MATLAB computeStiffnessMatrix."""
    K = np.zeros((4, 4))
    for i in range(4):
        if i < 3:
            K[i, i] = story_stiffness[i] + story_stiffness[i + 1]
            K[i, i + 1] = -story_stiffness[i + 1]
            K[i + 1, i] = -story_stiffness[i + 1]
        else:
            K[i, i] = story_stiffness[i]
    return K

def compute_rs_stiffness(last_Rs, last_disp, disp, fy, k0, b):
    """Exact MATLAB computeRsStiffness."""
    n = len(last_Rs)
    last_shear = np.zeros(n)
    last_rel_disp = np.zeros(n)
    rel_disp = np.zeros(n)
    shear = np.zeros(n)
    stiffness = np.zeros(n)

    for i in range(n):
        last_shear[i] = np.sum(last_Rs[i:n])

    last_rel_disp[0] = last_disp[0]
    rel_disp[0] = disp[0]
    for i in range(1, n):
        last_rel_disp[i] = last_disp[i] - last_disp[i-1]
        rel_disp[i] = disp[i] - disp[i-1]

    for i in range(n):
        fy_minus_b = fy[i] * (1 - b[i])
        k_sh = b[i] * k0[i]
        delta_disp = rel_disp[i] - last_rel_disp[i]
        c = last_shear[i] + k0[i] * delta_disp
        shear[i] = max(k_sh * rel_disp[i] - fy_minus_b,
                       min(k_sh * rel_disp[i] + fy_minus_b, c))
        if abs(shear[i] - c) < 1e-3:
            stiffness[i] = k0[i]
        else:
            stiffness[i] = k_sh

    Rs = np.zeros(n)
    Rs[3] = shear[3]
    for i in range(2, -1, -1):
        Rs[i] = shear[i] - shear[i+1]
    return Rs, stiffness

# Initial conditions
u = np.zeros(4)
u_dot = np.zeros(4)
ones = np.ones(4)
u_dot_dot = np.linalg.solve(M, -acc_data[0] * M @ ones - C @ u_dot)
Rs = np.zeros(4)

# Storage
top_disp = np.zeros(n_steps)
top_accel = np.zeros(n_steps)
iso_disp = np.zeros(n_steps)
base_shear = np.zeros(n_steps)

top_disp[0] = u[3]
top_accel[0] = u_dot_dot[3]
iso_disp[0] = u[0]
base_shear[0] = np.sum(Rs)

# Initial effective stiffness
K_eff = K0 + (1.0/(beta_n*dt*dt))*M + (gamma/(beta_n*dt))*C
k_story = k0.copy()
k_trial = k0.copy()

max_iter = 100
tol = 1e-5
conv_fails = 0

print(f"Starting analysis: {n_steps} steps...")

for j in range(n_steps - 1):
    # External load increment
    delta_P = -(acc_data[j+1] - acc_data[j]) * M @ ones

    # Effective load increment
    delta_P_eff = (delta_P
        + ((1.0/(beta_n*dt))*M + (gamma/beta_n)*C) @ u_dot
        + ((1.0/(2*beta_n))*M + dt*(gamma/(2*beta_n) - 1)*C) @ u_dot_dot)

    # Newton-Raphson
    u_trial = u.copy()
    Rs_trial = Rs.copy()
    delta_R_trial = delta_P_eff.copy()

    converged = False
    for iter in range(max_iter):
        if np.sum(np.abs(delta_R_trial)) < tol:
            converged = True
            break

        delta_u = np.linalg.solve(K_eff, delta_R_trial)
        u_trial_n = u_trial + delta_u

        Rs_trial_n, k_trial = compute_rs_stiffness(
            Rs_trial, u_trial, u_trial_n, fy, k0, b_hard)
        K_trial = compute_stiffness_matrix(k_trial)
        K_eff_trial = K_trial + (1.0/(beta_n*dt*dt))*M + (gamma/(beta_n*dt))*C

        delta_F = (Rs_trial_n - Rs_trial
            + (gamma/(beta_n*dt))*C @ delta_u
            + (1.0/(beta_n*dt*dt))*M @ delta_u)
        delta_R_trial = delta_R_trial - delta_F

        u_trial = u_trial_n
        Rs_trial = Rs_trial_n
        K_eff = K_eff_trial

    if not converged:
        conv_fails += 1

    delta_u = u_trial - u
    delta_u_dot = ((gamma/(beta_n*dt))*delta_u
        - (gamma/beta_n)*u_dot
        + (1 - gamma/(2*beta_n))*dt*u_dot_dot)
    delta_u_dot_dot = ((1.0/(beta_n*dt*dt))*delta_u
        - (1.0/(beta_n*dt))*u_dot
        - (1.0/(2*beta_n))*u_dot_dot)

    u = u_trial
    u_dot = u_dot + delta_u_dot
    u_dot_dot = u_dot_dot + delta_u_dot_dot
    Rs = Rs_trial

    top_disp[j+1] = u[3]
    top_accel[j+1] = u_dot_dot[3]
    iso_disp[j+1] = u[0]
    base_shear[j+1] = np.sum(Rs)

    if (j+1) % 200 == 0:
        print(f"  Step {j+1}: top_d={u[3]:.6f}")

print(f"\n=== Python MATLAB-equivalent Results ===")
print(f"Peak top displacement: {np.max(np.abs(top_disp)):.8f} m")
print(f"Peak top acceleration: {np.max(np.abs(top_accel)):.8f} m/s^2")
print(f"Peak iso displacement: {np.max(np.abs(iso_disp)):.8f} m")
print(f"Peak base shear: {np.max(np.abs(base_shear)):.4f} N")
print(f"Convergence failures: {conv_fails}/{n_steps}")

# Compare with MATLAB results
matlab_peak_top_disp = 0.0281488123
matlab_peak_top_accel = 5.079157
matlab_peak_iso_disp = 0.0231487
matlab_peak_base_shear = 532869.17

print(f"\n=== Comparison with MATLAB ===")
print(f"Top disp: py={np.max(np.abs(top_disp)):.8f} vs matlab={matlab_peak_top_disp:.8f} "
      f"(diff={abs(np.max(np.abs(top_disp))-matlab_peak_top_disp)/matlab_peak_top_disp*100:.4f}%)")
print(f"Top accel: py={np.max(np.abs(top_accel)):.8f} vs matlab={matlab_peak_top_accel:.8f} "
      f"(diff={abs(np.max(np.abs(top_accel))-matlab_peak_top_accel)/matlab_peak_top_accel*100:.4f}%)")
print(f"Iso disp: py={np.max(np.abs(iso_disp)):.8f} vs matlab={matlab_peak_iso_disp:.8f} "
      f"(diff={abs(np.max(np.abs(iso_disp))-matlab_peak_iso_disp)/matlab_peak_iso_disp*100:.4f}%)")
print(f"Base shear: py={np.max(np.abs(base_shear)):.4f} vs matlab={matlab_peak_base_shear:.4f} "
      f"(diff={abs(np.max(np.abs(base_shear))-matlab_peak_base_shear)/matlab_peak_base_shear*100:.4f}%)")

# Compare with OpenSees
ops_top_disp = 0.0246977
ops_iso_disp = 0.0188663
print(f"\n=== Python MATLAB-eq vs OpenSees ===")
print(f"Top disp: py={np.max(np.abs(top_disp)):.8f} vs ops={ops_top_disp:.8f} "
      f"(diff={abs(np.max(np.abs(top_disp))-ops_top_disp)/max(np.max(np.abs(top_disp)),1e-10)*100:.4f}%)")
print(f"Iso disp: py={np.max(np.abs(iso_disp)):.8f} vs ops={ops_iso_disp:.8f} "
      f"(diff={abs(np.max(np.abs(iso_disp))-ops_iso_disp)/max(np.max(np.abs(iso_disp)),1e-10)*100:.4f}%)")

# Save results
np.savetxt('topStoDisIso2_py.txt',
           np.column_stack([time, top_disp]), fmt='%.10e')
np.savetxt('topStoAccIso2_py.txt',
           np.column_stack([time, top_accel]), fmt='%.10e')
print("\nPython results saved.")
