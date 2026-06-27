#!/usr/bin/env python
"""Debug: compare OpenSees and MATLAB linear responses step by step."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

# Parameters
mass_iso = 250e3
mass_st = [270e3, 270e3, 180e3]
k_st = [245e6, 195e6, 98e6]
keq = 50e6
c_iso = 1000e3
dr = 0.05
dt = 0.01
g = 9.81
PGA = 0.4 * g

# Rayleigh coefficients (from MATLAB)
alpha_m = 0.505929
beta_k = 0.003477

acc = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc = acc / np.max(np.abs(acc)) * PGA
n = len(acc)

# --- OpenSees: Full model ---
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

ops.uniaxialMaterial('Viscous', 10, c_iso, 1.0)
ops.element('twoNodeLink', 5, 1, 2, '-mat', 10, '-dir', 1)

ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)

ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

# Record first 20 steps for comparison
ops_disp = np.zeros((20, 4))
for s in range(20):
    ops.analyze(1, dt)
    for dof in range(4):
        ops_disp[s, dof] = ops.nodeDisp(dof+2, 1)
ops.wipe()

# Save for MATLAB comparison
np.savetxt('debug_ops_disp.txt', ops_disp, fmt='%.12e')
print("First 20 steps saved to debug_ops_disp.txt")
