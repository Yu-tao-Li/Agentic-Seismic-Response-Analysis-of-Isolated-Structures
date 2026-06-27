"""
OpenSeesPy Nonlinear MDOF - Example 5 (v2 - equivalent force approach)
Uses equivalent nodal forces (like MATLAB) instead of UniformExcitation.
1D model, Hardening (b=0.05) + ElasticPP (b=0) materials.
"""
import openseespy.opensees as ops
import numpy as np
import json

# Parameters
m = np.array([250, 270, 270, 180]) * 1e3
k0_arr = np.array([50, 245, 195, 98]) * 1e6
fy = np.array([500, 1225, 975, 490]) * 1e3
b_hard = [0.05, 0.0, 0.0, 0.0]

dt = 0.01
g = 9.81
PGA = 0.40 * g
damp_ratio = 0.05
n_stories = 4

# Load ground motion
data_path = '../../input_data/Northridge_01_NO_968.txt'
acc_raw = np.loadtxt(data_path)
acc_norm = acc_raw / np.max(np.abs(acc_raw)) * PGA
n_steps = len(acc_raw)
time = np.arange(n_steps) * dt

# Rayleigh coefficients
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
print(f'w1={w1:.4f}, w2={w2:.4f}, alpha={alpha_m:.6f}, beta={beta_k:.6f}')

# Build 1D model (ndm=1, ndf=1)
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Nodes
for i in range(5):
    ops.node(i+1, 0.0)
ops.fix(1, 1)

# Mass
for i in range(n_stories):
    ops.mass(i+2, m[i])

# Materials
# Story 1: Hardening with bilinear kinematic hardening
Hkin1 = b_hard[0] * k0_arr[0] / (1.0 - b_hard[0])
ops.uniaxialMaterial('Hardening', 1, k0_arr[0], fy[0], 0.0, Hkin1)
# Stories 2-4: ElasticPP
for i in range(1, n_stories):
    epsY_i = fy[i] / k0_arr[i]
    ops.uniaxialMaterial('ElasticPP', i+1, k0_arr[i], epsY_i)

# zeroLength elements
for i in range(n_stories):
    ops.element('zeroLength', i+1, i+1, i+2, '-mat', i+1, '-dir', 1)

# Rayleigh damping (initial stiffness)
ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)

# Equivalent nodal forces: F_i(t) = -m_i * a_g(t)
ref_acc = -1.0  # reference value for plain load pattern
ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc_norm.tolist(),
               '-factor', 1.0)
ops.pattern('Plain', 1, 1)
for i in range(n_stories):
    ops.load(i+2, ref_acc * m[i])  # reference load = -m_i, timeSeries = a_g(t)

# Recorders
ops.recorder('Node', '-file', 'node_disp_v2.txt', '-time', '-node', 5, '-dof', 1, 'disp')
ops.recorder('Node', '-file', 'node_accel_v2.txt', '-time', '-node', 5, '-dof', 1, 'accel')
ops.recorder('Node', '-file', 'iso_disp_v2.txt', '-time', '-node', 2, '-dof', 1, 'disp')
ops.recorder('Element', '-file', 'ele_force_v2.txt', '-time', '-ele', 1, 'force')

# Analysis
ops.constraints('Transformation')
ops.numberer('Plain')
ops.system('UmfPack')
ops.test('NormUnbalance', 1e-5, 100, 0)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

conv_fails = 0
for step in range(n_steps):
    ok = ops.analyze(1, dt)
    if ok != 0:
        conv_fails += 1
        ops.algorithm('KrylovNewton')
        ok2 = ops.analyze(1, dt)
        if ok2 != 0:
            print(f'  Step {step}: FAILED')
        ops.algorithm('Newton')
    if step % 200 == 0:
        d5 = ops.nodeDisp(5, 1)
        print(f'  Step {step}: top_d={d5:.6f}')

ops.wipe()

# Read results
disp_d = np.loadtxt('node_disp_v2.txt')
accel_d = np.loadtxt('node_accel_v2.txt')
iso_d = np.loadtxt('iso_disp_v2.txt')
ele1_f = np.loadtxt('ele_force_v2.txt')

top_disp = disp_d[:, 1]
top_accel = accel_d[:, 1]
iso_disp = iso_d[:, 1]
base_shear = ele1_f[:, 1]  # force at node 1 of element 1

print(f'\n=== Results (v2 - equivalent forces) ===')
print(f'Peak top displacement: {np.max(np.abs(top_disp)):.8f} m')
print(f'Peak top acceleration: {np.max(np.abs(top_accel)):.8f} m/s^2')
print(f'Peak iso displacement: {np.max(np.abs(iso_disp)):.8f} m')
print(f'Peak base shear: {np.max(np.abs(base_shear)):.4f} N')
print(f'Convergence failures: {conv_fails}/{n_steps}')

# Compare
matlab_top = 0.0281488123
matlab_iso = 0.0231487
matlab_accel = 5.079157
matlab_bs = 532869.17
print(f'\nvs MATLAB top disp: diff={abs(np.max(np.abs(top_disp))-matlab_top)/matlab_top*100:.2f}%')
print(f'vs MATLAB iso disp: diff={abs(np.max(np.abs(iso_disp))-matlab_iso)/matlab_iso*100:.2f}%')
print(f'vs MATLAB top accel: diff={abs(np.max(np.abs(top_accel))-matlab_accel)/matlab_accel*100:.2f}%')
print(f'vs MATLAB base shear: diff={abs(np.max(np.abs(base_shear))-matlab_bs)/matlab_bs*100:.2f}%')

print('\nDone.')
