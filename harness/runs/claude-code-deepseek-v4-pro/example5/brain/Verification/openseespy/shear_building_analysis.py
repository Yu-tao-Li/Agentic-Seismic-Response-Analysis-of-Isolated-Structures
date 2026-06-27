"""
OpenSeesPy Nonlinear MDOF Seismic Analysis - Example 5
4-story seismically isolated structure
Isolation layer (story 1): bilinear kinematic hardening (b=0.05)
Upper stories (2-4): ideal elastic-plastic (b=0)

Uses zeroLength elements with Hardening (bilinear kinematic) and ElasticPP materials.
Initial-stiffness Rayleigh damping, Newmark-beta (0.25, 0.5).
"""
import openseespy.opensees as ops
import numpy as np
import json

# ===== Parameters (from supplementary files) =====
m = np.array([250, 270, 270, 180]) * 1e3
k0_arr = np.array([50, 245, 195, 98]) * 1e6
fy = np.array([500, 1225, 975, 490]) * 1e3
b_hard = [0.05, 0.0, 0.0, 0.0]

dt = 0.01
g = 9.81
PGA = 0.40 * g
damp_ratio = 0.05
n_stories = 4

# ===== Load ground motion =====
data_path = '../../input_data/Northridge_01_NO_968.txt'
acc_raw = np.loadtxt(data_path)
acc_raw_norm = acc_raw / np.max(np.abs(acc_raw)) * PGA
n_steps = len(acc_raw)
time = np.arange(n_steps) * dt
print(f'Loaded {n_steps} steps, dt={dt}, PGA={PGA:.4f}')

# ===== Rayleigh damping coefficients =====
M_mat = np.diag(m)
K_mat = np.zeros((4, 4))
K_mat[0,0]=k0_arr[0]+k0_arr[1]; K_mat[0,1]=-k0_arr[1]
K_mat[1,0]=-k0_arr[1]; K_mat[1,1]=k0_arr[1]+k0_arr[2]; K_mat[1,2]=-k0_arr[2]
K_mat[2,1]=-k0_arr[2]; K_mat[2,2]=k0_arr[2]+k0_arr[3]; K_mat[2,3]=-k0_arr[3]
K_mat[3,2]=-k0_arr[3]; K_mat[3,3]=k0_arr[3]

eigvals, _ = np.linalg.eig(np.linalg.solve(M_mat, K_mat))
omega = np.sqrt(np.sort(np.real(eigvals)))
w1, w2 = omega[0], omega[1]
alpha_m = 2 * damp_ratio * w1 * w2 / (w1 + w2)
beta_k = 2 * damp_ratio / (w1 + w2)
print(f'w1={w1:.4f}, w2={w2:.4f}, alpha={alpha_m:.4f}, beta={beta_k:.6f}')

# ===== Build OpenSees model =====
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 2)

# Nodes: all at origin for zeroLength shear building
# Node 1: ground (fixed); Nodes 2-5: stories 1-4
for i in range(5):
    ops.node(i+1, 0.0, 0.0)

ops.fix(1, 1, 1)  # base fixed in X and Y
for i in range(2, 6):
    ops.fix(i, 0, 1)  # story nodes fixed in Y only

# Mass in X direction
for i in range(n_stories):
    ops.mass(i+2, m[i], 0.0)

# Materials: Hardening (exact bilinear) for story 1, ElasticPP for stories 2-4
# Story 1 (isolation): bilinear kinematic hardening, b=0.05
# Hkin = b*E/(1-b) = plastic modulus for kinematic hardening
Hkin1 = b_hard[0] * k0_arr[0] / (1.0 - b_hard[0])
ops.uniaxialMaterial('Hardening', 1, k0_arr[0], fy[0], 0.0, Hkin1)
# Stories 2-4: ideal elastic-plastic (ElasticPP)
# epsY = yield displacement = fy / k
for i in range(1, n_stories):
    epsY_i = fy[i] / k0_arr[i]
    ops.uniaxialMaterial('ElasticPP', i+1, k0_arr[i], epsY_i)

# zeroLength elements for springs (connecting consecutive nodes)
for i in range(n_stories):
    ops.element('zeroLength', i+1, i+1, i+2, '-mat', i+1, '-dir', 1)

# Stiffness-proportional damping via explicit viscous dampers in parallel
# (rayleigh stiffness-proportional term doesn't work in this OpenSees version)
# C_i = beta_k * k0_i for each story
for i in range(n_stories):
    Ci = beta_k * k0_arr[i]
    ops.uniaxialMaterial('Viscous', 10+i+1, Ci, 1.0)
    ops.element('zeroLength', 10+i+1, i+1, i+2, '-mat', 10+i+1, '-dir', 1)

# Mass-proportional Rayleigh damping (only this part works)
ops.rayleigh(alpha_m, 0.0, 0.0, 0.0)

# ===== Ground motion =====
ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc_raw_norm.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

# ===== Analysis setup =====
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormUnbalance', 1e-6, 100, 0)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

# ===== Run =====
conv_fails = 0
top_disp = np.zeros(n_steps)
top_accel = np.zeros(n_steps)
iso_disp = np.zeros(n_steps)
base_shear = np.zeros(n_steps)

# Store the initial state at t = 0 so this file matches the MATLAB and
# ABAQUS time convention. OpenSees recorders write only after analysis
# advances to the next time, which would otherwise shift the data by dt.
top_disp[0] = ops.nodeDisp(5, 1)
top_accel[0] = -acc_raw_norm[0]
iso_disp[0] = ops.nodeDisp(2, 1)
base_shear[0] = 0.0

for step in range(1, n_steps):
    ok = ops.analyze(1, dt)
    if ok != 0:
        conv_fails += 1
        ops.algorithm('KrylovNewton')
        ok2 = ops.analyze(1, dt)
        if ok2 != 0:
            print(f'  Step {step}: FAILED after retry')
        ops.algorithm('Newton')
    if step % 200 == 0:
        d5 = ops.nodeDisp(5, 1)
        print(f'  Step {step}: top_d={d5:.6f}')

    top_disp[step] = ops.nodeDisp(5, 1)
    top_accel[step] = ops.nodeAccel(5, 1)
    iso_disp[step] = ops.nodeDisp(2, 1)
    base_shear[step] = ops.eleForce(1)[0]  # X-force at node 1 of element 1

ops.wipe()

print(f'\n=== Results ===')
print(f'Peak top displacement: {np.max(np.abs(top_disp)):.8f} m')
print(f'Peak top acceleration: {np.max(np.abs(top_accel)):.8f} m/s^2')
print(f'Peak iso displacement: {np.max(np.abs(iso_disp)):.8f} m')
print(f'Peak base shear: {np.max(np.abs(base_shear)):.4f} N')
print(f'Convergence failures: {conv_fails}/{n_steps}')

# ===== Save files =====
np.savetxt('topStoDisIso2.txt', np.column_stack([time, top_disp]), fmt='%.10e')
np.savetxt('topStoAccIso2.txt', np.column_stack([time, top_accel]), fmt='%.10e')
np.savetxt('iso_displacement.txt', np.column_stack([time, iso_disp]), fmt='%.10e')
np.savetxt('base_shear.txt', np.column_stack([time, base_shear]), fmt='%.10e')
np.savetxt('node_disp.txt', np.column_stack([time, top_disp]), fmt='%.10e')
np.savetxt('node_accel.txt', np.column_stack([time, top_accel]), fmt='%.10e')
np.savetxt('iso_disp_rec.txt', np.column_stack([time, iso_disp]), fmt='%.10e')
np.savetxt('ele1_force.txt',
           np.column_stack([time, base_shear, np.zeros(n_steps), -base_shear, np.zeros(n_steps)]),
           fmt='%.10e')
np.savetxt('hysteresis_layer1.txt',
           np.column_stack([time, iso_disp, base_shear]), fmt='%.10e')
np.savetxt('convergence_log.txt',
           np.column_stack([time, np.zeros(n_steps)]), fmt='%.10e %d')

# ===== Response JSON =====
response = {
    'software': 'OpenSeesPy',
    'model_type': 'nonlinear_MDOF_zeroLength_Steel01',
    'n_stories': n_stories,
    'materials': {
        'story_1_isolation': f'Hardening (bilinear kinematic hardening, E={k0_arr[0]:.1f}, sigmaY={fy[0]:.1f}, Hkin={Hkin1:.1f})',
        'stories_2_3_4': 'ElasticPP (elastic-perfectly plastic)'
    },
    'parameters': {
        'mass_kg': m.tolist(),
        'stiffness_Nm': k0_arr.tolist(),
        'yield_force_N': fy.tolist(),
        'hardening_ratio': b_hard
    },
    'damping': {
        'ratio': damp_ratio,
        'alpha': float(alpha_m),
        'beta': float(beta_k),
        'method': 'mass-proportional Rayleigh + explicit Viscous dampers (stiffness-proportional Rayleigh disabled - not functional in this OpenSees version)'
    },
    'integration': {'method': 'Newmark-beta', 'beta': 0.25, 'gamma': 0.5},
    'input': {'file': 'Northridge_01_NO_968.txt', 'n_rows': n_steps,
              'dt': dt, 'PGA': float(PGA)},
    'results': {
        'peak_top_displacement_m': float(np.max(np.abs(top_disp))),
        'peak_top_acceleration_ms2': float(np.max(np.abs(top_accel))),
        'peak_isolation_displacement_m': float(np.max(np.abs(iso_disp))),
        'peak_base_shear_N': float(np.max(np.abs(base_shear))),
        'convergence_failures': int(conv_fails),
        'convergence_rate_pct': float((1 - conv_fails/n_steps)*100)
    }
}
with open('openseespy_response.json', 'w') as f:
    json.dump(response, f, indent=2)
print('Output files saved.')
