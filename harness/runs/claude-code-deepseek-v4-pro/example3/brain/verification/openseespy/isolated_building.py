import numpy as np
import json
import openseespy.opensees as ops

# ============================================================
# Load and normalize ground motion
# ============================================================
gm = np.loadtxt('../../input_data/GM1.txt')
dt = 0.02
g = 9.81
target_pga = 0.40 * g
max_abs_gm = np.max(np.abs(gm))
scale = target_pga / max_abs_gm
ag = gm * scale
N = len(ag)
t = np.arange(N) * dt

print(f"GM1: N={N}, max_abs={max_abs_gm:.6f}, scale={scale:.6f}, target_PGA={target_pga:.4f}")

# ============================================================
# Compute Rayleigh damping coefficients (same as MATLAB)
# ============================================================
M_mat = np.diag([500., 1000., 1000., 1000.])
k_iso = 200e3
k_story = 500e3
K_mat = np.array([
    [k_iso+k_story, -k_story,       0.,        0.],
    [-k_story,      2*k_story, -k_story,        0.],
    [0.,           -k_story,  2*k_story, -k_story],
    [0.,              0.,     -k_story,   k_story]
])

# Solve generalized eigenvalue problem
eigvals, eigvecs = np.linalg.eig(np.linalg.solve(M_mat, K_mat))
omega = np.sqrt(np.sort(eigvals.real))
omega1, omega2 = omega[0], omega[1]

xi = 0.05
alpha_ray = 2*xi*omega1*omega2/(omega1+omega2)
beta_ray = 2*xi/(omega1+omega2)
c_iso = 2000.0

print(f"Natural frequencies (rad/s): {omega}")
print(f"Rayleigh: alpha={alpha_ray:.6f}, beta={beta_ray:.6e}")
print(f"Added isolation damping: c_iso={c_iso}")

# ============================================================
# Build OpenSeesPy model (1D shear building)
# ============================================================
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Nodes: 1=base(fixed), 2=isolation, 3=story2, 4=story3, 5=story4(top)
for i in range(1, 6):
    ops.node(i, 0.0)

# Fix base
ops.fix(1, 1)

# Masses
ops.mass(2, 500.0)   # isolation layer
ops.mass(3, 1000.0)  # story 2
ops.mass(4, 1000.0)  # story 3
ops.mass(5, 1000.0)  # story 4 (top)

# Elastic springs between adjacent floors
ops.uniaxialMaterial('Elastic', 1, k_iso)    # base-isolation
ops.uniaxialMaterial('Elastic', 2, k_story)  # iso-story2
ops.uniaxialMaterial('Elastic', 3, k_story)  # story2-story3
ops.uniaxialMaterial('Elastic', 4, k_story)  # story3-story4

ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1, '-doRayleigh', 1)
ops.element('zeroLength', 2, 2, 3, '-mat', 2, '-dir', 1, '-doRayleigh', 1)
ops.element('zeroLength', 3, 3, 4, '-mat', 3, '-dir', 1, '-doRayleigh', 1)
ops.element('zeroLength', 4, 4, 5, '-mat', 4, '-dir', 1, '-doRayleigh', 1)

# Additional viscous damper at isolation level (C(1,1) += c_iso)
ops.uniaxialMaterial('Viscous', 10, c_iso, 1.0)  # C, alpha=1.0 for linear
ops.element('zeroLength', 10, 1, 2, '-mat', 10, '-dir', 1)

# Apply Rayleigh damping
ops.rayleigh(alpha_ray, beta_ray, 0.0, 0.0)

# ============================================================
# Define ground motion time series and load pattern
# ============================================================
time_series_tag = 1
ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *ag.tolist(), '-factor', 1.0)

# Equivalent forces: F_i = -m_i * a_g applied at each floor
# The ground motion is already in the time series, so use magnitude -mass
ops.pattern('Plain', 1, time_series_tag)
ops.load(2, -500.0)   # -m1 * a_g
ops.load(3, -1000.0)  # -m2 * a_g
ops.load(4, -1000.0)  # -m3 * a_g
ops.load(5, -1000.0)  # -m4 * a_g

# ============================================================
# Analysis setup
# ============================================================
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('ProfileSPD')
ops.test('NormDispIncr', 1.0e-8, 100, 0)
ops.algorithm('Linear')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

# ============================================================
# Run analysis and record results
# ============================================================
top_dis = np.zeros(N)
top_vel = np.zeros(N)
top_acc_rel = np.zeros(N)
top_acc_abs = np.zeros(N)
iso_dis = np.zeros(N)

# Store the initial state at t = 0 so the OpenSeesPy output uses the
# same time convention as the MATLAB and ABAQUS files.
top_dis[0] = ops.nodeDisp(5, 1)
top_vel[0] = ops.nodeVel(5, 1)
top_acc_rel[0] = -ag[0]
top_acc_abs[0] = top_acc_rel[0] + ag[0]
iso_dis[0] = ops.nodeDisp(2, 1)

for i in range(1, N):
    ok = ops.analyze(1, dt)
    if ok != 0:
        print(f"Analysis failed at step {i}")
        break
    top_dis[i] = ops.nodeDisp(5, 1)
    top_vel[i] = ops.nodeVel(5, 1)
    top_acc_rel[i] = ops.nodeAccel(5, 1)        # relative acceleration
    top_acc_abs[i] = top_acc_rel[i] + ag[i]      # absolute = relative + ground
    iso_dis[i] = ops.nodeDisp(2, 1)

ops.wipe()

# ============================================================
# Output results
# ============================================================
print(f"Peak top displacement: {np.max(np.abs(top_dis)):.6f} m")
print(f"Peak top rel acceleration: {np.max(np.abs(top_acc_rel)):.6f} m/s^2")
print(f"Peak top abs acceleration: {np.max(np.abs(top_acc_abs)):.6f} m/s^2 ({np.max(np.abs(top_acc_abs))/g:.4f} g)")
print(f"Peak isolation displacement: {np.max(np.abs(iso_dis)):.6f} m")

# topStoDis.txt
np.savetxt('topStoDis.txt', np.column_stack([t, top_dis]), fmt='%.10e')

# topStoAcc.txt
np.savetxt('topStoAcc.txt', np.column_stack([t, top_acc_rel]), fmt='%.10e')

# isoDis.txt
np.savetxt('isoDis.txt', np.column_stack([t, iso_dis]), fmt='%.10e')

# Response JSON
resp = {
    'time': t.tolist(),
    'top_displacement': top_dis.tolist(),
    'top_rel_acceleration': top_acc_rel.tolist(),
    'top_abs_acceleration': top_acc_abs.tolist(),
    'isolation_displacement': iso_dis.tolist(),
    'input_pga': float(target_pga),
    'dt': dt,
    'n_samples': N,
    'peak_top_dis': float(np.max(np.abs(top_dis))),
    'peak_top_rel_acc': float(np.max(np.abs(top_acc_rel))),
    'peak_top_abs_acc': float(np.max(np.abs(top_acc_abs))),
    'peak_iso_dis': float(np.max(np.abs(iso_dis))),
    'omega': omega.tolist(),
    'alpha': float(alpha_ray),
    'beta': float(beta_ray),
    'c_iso_added': c_iso
}

with open('response.json', 'w') as f:
    json.dump(resp, f, indent=2)

print('OpenSeesPy analysis complete. Output files written.')
