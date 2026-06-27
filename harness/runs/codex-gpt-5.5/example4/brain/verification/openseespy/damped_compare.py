#!/usr/bin/env python
"""Step-by-step damping comparison."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

mass_iso = 250e3
mass_st = [270e3, 270e3, 180e3]
k_st = [245e6, 195e6, 98e6]
keq = 50e6
c_iso = 1000e3
dt = 0.01
g = 9.81; PGA = 0.4 * g

alpha_m = 0.505929
beta_k = 0.003477

acc = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc = acc / np.max(np.abs(acc)) * PGA
n = len(acc)
ones_v = np.ones(4)

# MATLAB equivalent damping matrix
M = np.diag([mass_iso] + mass_st)
K = np.zeros((4, 4))
K[0,0] = keq + k_st[0]; K[0,1] = -k_st[0]; K[1,0] = -k_st[0]
for i in range(1, 3):
    K[i,i] = k_st[i-1] + k_st[i]
    K[i,i-1] = -k_st[i-1]; K[i,i+1] = -k_st[i]
K[3,3] = k_st[2]; K[3,2] = -k_st[2]

C_matlab = alpha_m * M + beta_k * K
C_matlab[0,0] += c_iso

beta_n = 0.25; gamma_n = 0.5
c0 = 1/(beta_n*dt**2); c1 = gamma_n/(beta_n*dt); c2 = 1/(beta_n*dt)
c3 = 1/(2*beta_n)-1; c4 = gamma_n/beta_n-1; c5 = dt*(gamma_n/(2*beta_n)-1)

# Test 1: Rayleigh only (MATLAB)
K_eff_m1 = K + c1*(alpha_m*M + beta_k*K) + c0*M  # no c_iso
u = np.zeros(4); v = np.zeros(4); a = -ones_v*acc[0]
dmax1 = 0; amax1 = 0
for i in range(n):
    ag = acc[i]
    F_eff = -M@ones_v*ag + M@(c0*u + c2*v + c3*a) + (alpha_m*M + beta_k*K)@(c1*u + c4*v + c5*a)
    u_new = np.linalg.solve(K_eff_m1, F_eff)
    a_new = c0*(u_new-u) - c2*v - c3*a
    v_new = v + dt*((1-gamma_n)*a + gamma_n*a_new)
    u, v, a = u_new, v_new, a_new
    dmax1 = max(dmax1, abs(u[3])); amax1 = max(amax1, abs(a[3]))
print(f"MATLAB (Rayleigh only):    d_max={dmax1:.6f}, a_max={amax1:.6f}")

# Test 2: Rayleigh only (OpenSees)
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
for i in range(1, 6): ops.node(i, 0.0)
ops.fix(1, 1)
ops.mass(2, mass_iso); ops.mass(3, mass_st[0]); ops.mass(4, mass_st[1]); ops.mass(5, mass_st[2])
ops.uniaxialMaterial('Elastic', 1, keq)
for i in range(3): ops.uniaxialMaterial('Elastic', i+2, k_st[i])
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.element('twoNodeLink', 2, 2, 3, '-mat', 2, '-dir', 1)
ops.element('twoNodeLink', 3, 3, 4, '-mat', 3, '-dir', 1)
ops.element('twoNodeLink', 4, 4, 5, '-mat', 4, '-dir', 1)
# ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)  # init stiff prop
ops.rayleigh(alpha_m, beta_k, 0.0, 0.0)    # tangent stiff prop
ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')
dmax2 = 0; amax2 = 0
for s in range(n):
    ops.analyze(1, dt)
    d = ops.nodeDisp(5, 1); a = ops.nodeAccel(5, 1)
    dmax2 = max(dmax2, abs(d)); amax2 = max(amax2, abs(a))
ops.wipe()
print(f"OpenSees (tangent stiff damp): d_max={dmax2:.6f}, a_max={amax2:.6f}")

# Test 3: Rayleigh only (OpenSees) - init stiffness
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
for i in range(1, 6): ops.node(i, 0.0)
ops.fix(1, 1)
ops.mass(2, mass_iso); ops.mass(3, mass_st[0]); ops.mass(4, mass_st[1]); ops.mass(5, mass_st[2])
ops.uniaxialMaterial('Elastic', 1, keq)
for i in range(3): ops.uniaxialMaterial('Elastic', i+2, k_st[i])
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.element('twoNodeLink', 2, 2, 3, '-mat', 2, '-dir', 1)
ops.element('twoNodeLink', 3, 3, 4, '-mat', 3, '-dir', 1)
ops.element('twoNodeLink', 4, 4, 5, '-mat', 4, '-dir', 1)
ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)  # init stiff prop
ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')
dmax3 = 0; amax3 = 0
for s in range(n):
    ops.analyze(1, dt)
    d = ops.nodeDisp(5, 1); a = ops.nodeAccel(5, 1)
    dmax3 = max(dmax3, abs(d)); amax3 = max(amax3, abs(a))
ops.wipe()
print(f"OpenSees (init stiff damp):   d_max={dmax3:.6f}, a_max={amax3:.6f}")

# Test 4: Full damping MATLAB
K_eff_m2 = K + c1*C_matlab + c0*M
u = np.zeros(4); v = np.zeros(4); a = -ones_v*acc[0]
dmax4 = 0; amax4 = 0
for i in range(n):
    ag = acc[i]
    F_eff = -M@ones_v*ag + M@(c0*u + c2*v + c3*a) + C_matlab@(c1*u + c4*v + c5*a)
    u_new = np.linalg.solve(K_eff_m2, F_eff)
    a_new = c0*(u_new-u) - c2*v - c3*a
    v_new = v + dt*((1-gamma_n)*a + gamma_n*a_new)
    u, v, a = u_new, v_new, a_new
    dmax4 = max(dmax4, abs(u[3])); amax4 = max(amax4, abs(a[3]))
print(f"MATLAB (full damp):          d_max={dmax4:.6f}, a_max={amax4:.6f}")

# Test 5: OpenSees with c_iso via equivalent viscous stiffness in elastic element
# Instead of Viscous material, use the fact that for Newmark, a dashpot acts like a spring
# with effective stiffness c_iso * c1. Add this as an additional elastic spring.
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
for i in range(1, 6): ops.node(i, 0.0)
ops.fix(1, 1)
ops.mass(2, mass_iso); ops.mass(3, mass_st[0]); ops.mass(4, mass_st[1]); ops.mass(5, mass_st[2])

# Isolation: elastic spring + equivalent dashpot stiffness
ops.uniaxialMaterial('Elastic', 1, keq + c_iso*c1)  # keq + c_iso*c1 in stiffness

# For damping force: add it as effective load
# Actually, we need: C_matlab * v_new = C_rayleigh*v_new + c_iso*v1_new
# In effective force formulation: c_iso * (c1*u_new + c4*v_n + c5*a_n)
# The c_iso*c1*u_new part goes to K_eff
# The c_iso*(c4*v_n + c5*a_n) part goes to F_eff

# Hmm, this is getting complicated. Let me try another approach:
# Use the fact that in MATLAB, F_eff_damp = C_matlab * (c1*u_n + c4*v_n + c5*a_n)
# The c1*C_matlab*u_n part is included in effective force...

# Actually, the cleanest approach is to directly implement the MATLAB formulation
# using a load pattern instead of relying on OpenSees internal damping handling.

# For now, let me just see the difference between these tests
print(f"\nComparison:")
print(f"  MATLAB Rayleigh-only:        {dmax1:.6f}, {amax1:.6f}")
print(f"  OpenSees tangent-stiff damp: {dmax2:.6f}, {amax2:.6f}")
print(f"  OpenSees init-stiff damp:    {dmax3:.6f}, {amax3:.6f}")
print(f"  MATLAB full damp:            {dmax4:.6f}, {amax4:.6f}")
