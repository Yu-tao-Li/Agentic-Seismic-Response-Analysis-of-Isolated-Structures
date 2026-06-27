#!/usr/bin/env python
"""Diagnose OpenSees Rayleigh damping implementation with a simple 1-DOF test
and then the full 4-DOF system step-by-step."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

# ============ 1-DOF TEST ============
print("=== 1-DOF Rayleigh Damping Test ===")
m1 = 1.0
k1 = 100.0
omega_n = np.sqrt(k1/m1)  # = 10 rad/s
alpha = 2*0.05*omega_n     # mass-proportional: 2*xi*omega for single mode = alpha only
beta = 2*0.05/omega_n      # stiffness-proportional: 2*xi/omega

print(f"omega_n = {omega_n:.4f}")
print(f"alpha = {alpha:.6f}, beta = {beta:.6f}")

# MATLAB-style: C = alpha*M + beta*K = 1.0 + 0.01 = 1.01 (for each DOF)
c_matlab_style = alpha*m1 + beta*k1
print(f"Expected C = {c_matlab_style:.6f}")

dt = 0.01
F0 = 100.0  # step force

# Test 1a: OpenSees with Rayleigh
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
ops.node(1, 0.0); ops.node(2, 0.0)
ops.fix(1, 1)
ops.mass(2, m1)
ops.uniaxialMaterial('Elastic', 1, k1)
ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.rayleigh(alpha, 0.0, beta, 0.0)  # init stiffness
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, F0)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

t_end = 2.0
n_steps = int(t_end/dt)
t_os = []; d_os = []; v_os = []; a_os = []
for i in range(n_steps):
    ops.analyze(1, dt)
    t_os.append(ops.getTime())
    d_os.append(ops.nodeDisp(2, 1))
    v_os.append(ops.nodeVel(2, 1))
    a_os.append(ops.nodeAccel(2, 1))
ops.wipe()

# Test 1b: Python Newmark with same C
d_py = [0.0]; v_py = [0.0]; a_py = [0.0]
c = c_matlab_style

beta_nm = 0.25; gamma_nm = 0.5
c0 = 1/(beta_nm*dt**2); c1 = gamma_nm/(beta_nm*dt); c2 = 1/(beta_nm*dt)
c3 = 1/(2*beta_nm)-1; c4 = gamma_nm/beta_nm-1; c5 = dt*(gamma_nm/(2*beta_nm)-1)

K_eff = k1 + c1*c + c0*m1

for i in range(n_steps):
    u_n, v_n, a_n = d_py[-1], v_py[-1], a_py[-1]
    F_eff = F0 + m1*(c0*u_n + c2*v_n + c3*a_n) + c*(c1*u_n + c4*v_n + c5*a_n)
    u_new = F_eff / K_eff
    a_new = c0*(u_new - u_n) - c2*v_n - c3*a_n
    v_new = v_n + dt*((1-gamma_nm)*a_n + gamma_nm*a_new)
    d_py.append(u_new); v_py.append(v_new); a_py.append(a_new)

d_py = d_py[1:]; v_py = v_py[1:]; a_py = a_py[1:]

# Compare
d_diff = np.max(np.abs(np.array(d_os) - np.array(d_py)))
v_diff = np.max(np.abs(np.array(v_os) - np.array(v_py)))
print(f"Max displacement diff: {d_diff:.8f}")
print(f"Max velocity diff:     {v_diff:.8f}")

if d_diff < 1e-6:
    print("1-DOF: PASS - OpenSees Rayleigh matches Python reference")
else:
    print(f"1-DOF: FAIL - OpenSees Rayleigh differs from Python reference!")
    # Show first few steps
    print("\nFirst 5 steps comparison:")
    print("Step  OS_disp    PY_disp    OS_vel     PY_vel")
    for i in range(5):
        print(f"{i}: {d_os[i]:.8f} {d_py[i]:.8f} {v_os[i]:.8f} {v_py[i]:.8f}")

# ============ 4-DOF TEST - STEP-BY-STEP ============
print("\n=== 4-DOF Rayleigh Damping Step-by-Step ===")
mass_iso = 250e3
mass_st = [270e3, 270e3, 180e3]
k_st = [245e6, 195e6, 98e6]
keq = 50e6
c_iso = 1000e3  # NOT used in Rayleigh-only test

alpha_m = 0.505929
beta_k = 0.003477

M = np.diag([mass_iso] + mass_st)
K = np.zeros((4, 4))
K[0,0] = keq + k_st[0]; K[0,1] = -k_st[0]; K[1,0] = -k_st[0]
for i in range(1, 3):
    K[i,i] = k_st[i-1] + k_st[i]
    K[i,i-1] = -k_st[i-1]; K[i,i+1] = -k_st[i]
K[3,3] = k_st[2]; K[3,2] = -k_st[2]

C_matlab = alpha_m * M + beta_k * K  # NO c_iso

g = 9.81; PGA = 0.4 * g

acc = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc = acc / np.max(np.abs(acc)) * PGA

# OpenSees: Rayleigh only, NO viscous, NO Steel01 (all Elastic), first 20 steps
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
for i in range(1, 6):
    ops.node(i, 0.0)
ops.fix(1, 1)
ops.mass(2, mass_iso); ops.mass(3, mass_st[0]); ops.mass(4, mass_st[1]); ops.mass(5, mass_st[2])
ops.uniaxialMaterial('Elastic', 1, keq)
for i in range(3):
    ops.uniaxialMaterial('Elastic', i+2, k_st[i])
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.element('twoNodeLink', 2, 2, 3, '-mat', 2, '-dir', 1)
ops.element('twoNodeLink', 3, 3, 4, '-mat', 3, '-dir', 1)
ops.element('twoNodeLink', 4, 4, 5, '-mat', 4, '-dir', 1)
ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)
ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

os_d = np.zeros((20, 4)); os_v = np.zeros((20, 4)); os_a = np.zeros((20, 4))
for s in range(20):
    ops.analyze(1, dt)
    for dof in range(4):
        os_d[s, dof] = ops.nodeDisp(dof+2, 1)
        os_v[s, dof] = ops.nodeVel(dof+2, 1)
        os_a[s, dof] = ops.nodeAccel(dof+2, 1)
ops.wipe()

# Python Newmark with MATLAB formulation
beta_nm = 0.25; gamma_nm = 0.5
c0 = 1/(beta_nm*dt**2); c1 = gamma_nm/(beta_nm*dt); c2 = 1/(beta_nm*dt)
c3 = 1/(2*beta_nm)-1; c4 = gamma_nm/beta_nm-1; c5 = dt*(gamma_nm/(2*beta_nm)-1)

K_eff_py = K + c1*C_matlab + c0*M
K_eff_py_inv = np.linalg.inv(K_eff_py)

u = np.zeros(4); v = np.zeros(4); a = -np.ones(4) * acc[0]
ones_v = np.ones(4)

py_d = np.zeros((20, 4)); py_v = np.zeros((20, 4)); py_a = np.zeros((20, 4))
for i in range(20):
    ag = acc[i]
    F_eff = -M @ ones_v * ag + M @ (c0*u + c2*v + c3*a) + C_matlab @ (c1*u + c4*v + c5*a)
    u_new = K_eff_py_inv @ F_eff
    a_new = c0*(u_new - u) - c2*v - c3*a
    v_new = v + dt*((1-gamma_nm)*a + gamma_nm*a_new)
    py_d[i] = u_new; py_v[i] = v_new; py_a[i] = a_new
    u, v, a = u_new, v_new, a_new

# Compare step by step
print("Step-by-step displacement comparison (top DOF, first 20 steps):")
print("Step  OS_disp    PY_disp    diff")
max_diff = 0
max_diff_step = 0
for i in range(20):
    diff = abs(os_d[i, 3] - py_d[i, 3])
    if diff > max_diff:
        max_diff = diff
        max_diff_step = i
    print(f"{i:3d}: {os_d[i,3]:.8f} {py_d[i,3]:.8f} {diff:.2e}")

print(f"\nMax diff at step {max_diff_step}: {max_diff:.6e}")

# Also compare K_eff: compute what OpenSees uses vs MATLAB
# We can't directly query OpenSees's K_eff, but we can compute
# the effective stiffness from the first step's response
# Check: does first step match?
print(f"\nFirst step displacement (all DOFs):")
print(f"  OpenSees: {os_d[0]}")
print(f"  Python:   {py_d[0]}")

# Check Rayleigh matrix components
print(f"\nRayleigh C matrix diagonal (MATLAB):")
print(f"  diag(C_matlab) = {np.diag(C_matlab)}")
print(f"  c1*C(1,1) = {c1*C_matlab[0,0]:.2f}")
print(f"  keq = {keq:.2f}")
print(f"  c1*beta_k*keq = {c1*beta_k*keq:.2f}")

# Check: what if OpenSees uses c1 different?
# Let's find what c1 makes OpenSees match MATLAB
# by computing the effective stiffness from the first step
# u0_os * K_eff_os = F_eff_0
# u0_os = os_d[0], F_eff_0 = -M*ones*acc[0] + ...
# Actually we can't easily invert this for general 4x4

# Let's try a different approach: what if OpenSees Rayleigh uses K_committed?
# Actually, the initial and committed stiffness should be the same at step 1
# since no yielding. So this shouldn't matter.
print("\nDone.")
