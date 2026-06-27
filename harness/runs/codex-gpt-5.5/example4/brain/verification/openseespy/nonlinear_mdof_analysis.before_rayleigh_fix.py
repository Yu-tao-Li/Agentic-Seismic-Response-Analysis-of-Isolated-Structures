#!/usr/bin/env python
"""OpenSeesPy nonlinear MDOF time-history analysis for 4-story isolated structure.
Upper 3 stories: ideal elastic-plastic (Steel01, b=0)
Isolation layer: linear (keq + dashpot c)
Rayleigh damping: 5%
Note: OpenSees Rayleigh damping implementation differs slightly from MATLAB.
"""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

# Parameters from supplementary/Example4_2.m
mass_isolation = 250e3    # kg
mass_stories = [270e3, 270e3, 180e3]  # kg
k_stories = [245e6, 195e6, 98e6]  # N/m (k2, k3, k4)
fy_stories = [1225e3, 975e3, 490e3]  # N (fy2, fy3, fy4)
keq = 50e6  # N/m (isolation stiffness)
c_iso = 1000e3  # N/(m/s) (isolation damping)
damping_ratio = 0.05
dt = 0.01
g = 9.81
PGA = 0.40 * g

# Load ground motion
accel_raw = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
accel_raw = accel_raw / np.max(np.abs(accel_raw)) * PGA
n_steps = len(accel_raw)

# Eigenvalues for Rayleigh coefficients (computed from M, K0 matrices)
n_dof = 4
M_mat = np.diag([mass_isolation] + mass_stories)
K_mat = np.zeros((n_dof, n_dof))
K_mat[0,0] = keq + k_stories[0]
K_mat[0,1] = -k_stories[0]; K_mat[1,0] = -k_stories[0]
for i in range(1, n_dof-1):
    K_mat[i,i] = k_stories[i-1] + k_stories[i]
    K_mat[i,i-1] = -k_stories[i-1]; K_mat[i,i+1] = -k_stories[i]
K_mat[-1,-1] = k_stories[-1]; K_mat[-1,-2] = -k_stories[-1]

eigvals, _ = np.linalg.eig(np.linalg.solve(M_mat, K_mat))
omega_all = np.sqrt(np.sort(np.abs(eigvals)))
omega1, omega2 = omega_all[0], omega_all[1]

alpha_m = 2*omega1*omega2*(damping_ratio*omega2 - damping_ratio*omega1)/(omega2**2 - omega1**2)
beta_k = 2*(damping_ratio*omega2 - damping_ratio*omega1)/(omega2**2 - omega1**2)

print(f"omega1={omega1:.4f}, omega2={omega2:.4f}")
print(f"alpha_m={alpha_m:.6f}, beta_k={beta_k:.6f}")

# Build OpenSees model
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Nodes: 1=base(fixed), 2=isolation, 3=story2, 4=story3, 5=story4(top)
for i in range(1, 6):
    ops.node(i, 0.0)
ops.fix(1, 1)

# Masses (kg)
ops.mass(2, mass_isolation)
ops.mass(3, mass_stories[0])
ops.mass(4, mass_stories[1])
ops.mass(5, mass_stories[2])

# Materials
# Isolation: elastic spring keq
ops.uniaxialMaterial('Elastic', 1, keq)

# Upper stories: Steel01 (ideal elastic-plastic, b=0)
# Steel01 tag Fy E0 b
for i in range(3):
    ops.uniaxialMaterial('Steel01', i+2, fy_stories[i], k_stories[i], 0.0)

# Elements: twoNodeLink between consecutive floors
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)  # isolation spring
ops.element('twoNodeLink', 2, 2, 3, '-mat', 2, '-dir', 1)  # story 2
ops.element('twoNodeLink', 3, 3, 4, '-mat', 3, '-dir', 1)  # story 3
ops.element('twoNodeLink', 4, 4, 5, '-mat', 4, '-dir', 1)  # story 4

# Isolation dashpot (parallel to isolation spring)
ops.uniaxialMaterial('Viscous', 10, c_iso, 1.0)
ops.element('twoNodeLink', 5, 1, 2, '-mat', 10, '-dir', 1)

# Rayleigh damping: initial-stiffness-proportional
ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)

# Ground motion
ops.timeSeries('Path', 1, '-dt', dt, '-values', *accel_raw.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

# Analysis settings
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-8, 50, 0)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

# Run and record
rec_time, rec_disp_top, rec_acc_top, rec_disp_iso = [], [], [], []
rec_reaction, rec_yield = [], []
failed_steps = 0

for step in range(n_steps):
    ok = ops.analyze(1, dt)
    if ok != 0:
        failed_steps += 1

    t = ops.getTime()
    rec_time.append(t)
    rec_disp_top.append(ops.nodeDisp(5, 1))
    rec_acc_top.append(ops.nodeAccel(5, 1))
    rec_disp_iso.append(ops.nodeDisp(2, 1))
    rec_reaction.append(ops.nodeReaction(1, 1))

rec_time = np.array(rec_time)
rec_disp_top = np.array(rec_disp_top)
rec_acc_top = np.array(rec_acc_top)
rec_disp_iso = np.array(rec_disp_iso)
rec_reaction = np.array(rec_reaction)

# Save results
np.savetxt('topStoDisIso1.txt', np.column_stack([rec_time, rec_disp_top]), fmt='%.10e')
np.savetxt('topStoAccIso1.txt', np.column_stack([rec_time, rec_acc_top]), fmt='%.10e')
np.savetxt('iso_displacement.txt', np.column_stack([rec_time, rec_disp_iso]), fmt='%.10e')
np.savetxt('hysteresis_layer1.txt',
           np.column_stack([rec_time, rec_reaction, rec_disp_iso]), fmt='%.10e')

print(f"\nAnalysis complete.")
print(f"Steps: {n_steps}, Failed: {failed_steps}")
print(f"Max top displacement: {np.max(np.abs(rec_disp_top)):.6f} m")
print(f"Max top acceleration: {np.max(np.abs(rec_acc_top)):.6f} m/s^2")
print(f"Max isolation displacement: {np.max(np.abs(rec_disp_iso)):.6f} m")

ops.wipe()
