#!/usr/bin/env python
"""Undamped comparison: OpenSees vs MATLAB theory."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

mass_iso = 250e3
mass_st = [270e3, 270e3, 180e3]
k_st = [245e6, 195e6, 98e6]
keq = 50e6
dt = 0.01
g = 9.81; PGA = 0.4 * g

acc = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc = acc / np.max(np.abs(acc)) * PGA
n = len(acc)

# --- OpenSees: NO damping ---
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

# NO rayleigh, NO viscous

ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

dmax_ops = 0; amax_ops = 0
for s in range(n):
    ops.analyze(1, dt)
    d = ops.nodeDisp(5, 1); a = ops.nodeAccel(5, 1)
    dmax_ops = max(dmax_ops, abs(d)); amax_ops = max(amax_ops, abs(a))
ops.wipe()

# MATLAB equivalent (compute in Python)
M = np.diag([mass_iso] + mass_st)
K = np.zeros((4, 4))
K[0,0] = keq + k_st[0]; K[0,1] = -k_st[0]; K[1,0] = -k_st[0]
for i in range(1, 3):
    K[i,i] = k_st[i-1] + k_st[i]
    K[i,i-1] = -k_st[i-1]; K[i,i+1] = -k_st[i]
K[3,3] = k_st[2]; K[3,2] = -k_st[2]

beta_n = 0.25; gamma_n = 0.5
c0 = 1/(beta_n*dt**2)
c1 = gamma_n/(beta_n*dt)
c2 = 1/(beta_n*dt)
c3 = 1/(2*beta_n) - 1

K_eff_mat = K + c0*M  # no damping
dmax_mat = 0; amax_mat = 0

u = np.zeros(4); v = np.zeros(4); a = -np.ones(4) * acc[0];
ones_v = np.ones(4)

for i in range(n):
    ag = acc[i]
    F_eff = -M @ ones_v * ag + M @ (c0*u + c2*v + c3*a)
    u_new = np.linalg.solve(K_eff_mat, F_eff)
    a_new = c0*(u_new - u) - c2*v - c3*a
    v_new = v + dt*((1-gamma_n)*a + gamma_n*a_new)
    u, v, a = u_new, v_new, a_new
    dmax_mat = max(dmax_mat, abs(u[3]))
    amax_mat = max(amax_mat, abs(a[3]))

print(f"OpenSees (undamped): d_max={dmax_ops:.6f}, a_max={amax_ops:.6f}")
print(f"MATLAB (undamped):   d_max={dmax_mat:.6f}, a_max={amax_mat:.6f}")
