#!/usr/bin/env python
"""Verify: mass Rayleigh + Viscous elements = correct damping for 4-DOF system."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

mass_iso = 250e3; mass_st = [270e3, 270e3, 180e3]
k_st = [245e6, 195e6, 98e6]; keq = 50e6
alpha_m = 0.505929; beta_k = 0.003477
dt = 0.01; g = 9.81; PGA = 0.4 * g

acc = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc = acc / np.max(np.abs(acc)) * PGA
n = len(acc)

# MATLAB reference
M = np.diag([mass_iso] + mass_st)
K = np.zeros((4, 4))
K[0,0] = keq + k_st[0]; K[0,1] = -k_st[0]; K[1,0] = -k_st[0]
for i in range(1, 3):
    K[i,i] = k_st[i-1] + k_st[i]
    K[i,i-1] = -k_st[i-1]; K[i,i+1] = -k_st[i]
K[3,3] = k_st[2]; K[3,2] = -k_st[2]

# Compute Rayleigh coefficients from eigenvalues (for verification)
eigvals, _ = np.linalg.eig(np.linalg.solve(M, K))
omega = np.sqrt(np.sort(np.abs(eigvals)))
omega1, omega2 = omega[0], omega[1]
alpha_calc = 2*omega1*omega2*(0.05*omega2 - 0.05*omega1)/(omega2**2 - omega1**2)
beta_calc = 2*(0.05*omega2 - 0.05*omega1)/(omega2**2 - omega1**2)
print(f"omega1={omega1:.4f}, omega2={omega2:.4f}")
print(f"alpha_calc={alpha_calc:.6f}, beta_calc={beta_calc:.6f}")
print(f"(Precomputed: alpha={alpha_m:.6f}, beta={beta_k:.6f})")

# Use the computed values to ensure consistency
C_matlab = alpha_calc*M + beta_calc*K

beta_nm = 0.25; gamma_nm = 0.5
c0 = 1/(beta_nm*dt**2); c1 = gamma_nm/(beta_nm*dt); c2 = 1/(beta_nm*dt)
c3 = 1/(2*beta_nm)-1; c4 = gamma_nm/beta_nm-1; c5 = dt*(gamma_nm/(2*beta_nm)-1)

K_eff = K + c1*C_matlab + c0*M
K_eff_inv = np.linalg.inv(K_eff)

u = np.zeros(4); v = np.zeros(4); a = -np.ones(4) * acc[0]
ones_v = np.ones(4)

dmax_mat = 0; amax_mat = 0
for i in range(n):
    ag = acc[i]
    F_eff = -M@ones_v*ag + M@(c0*u + c2*v + c3*a) + C_matlab@(c1*u + c4*v + c5*a)
    u_new = K_eff_inv @ F_eff
    a_new = c0*(u_new - u) - c2*v - c3*a
    v_new = v + dt*((1-gamma_nm)*a + gamma_nm*a_new)
    u, v, a = u_new, v_new, a_new
    dmax_mat = max(dmax_mat, abs(u[3]))
    amax_mat = max(amax_mat, abs(a[3]))

print(f"\nMATLAB reference: d_max={dmax_mat:.6f}, a_max={amax_mat:.6f}")

# OpenSees: mass Rayleigh ONLY + Viscous elements for stiffness-proportional
# Viscous coefficients
c_vis_iso = beta_calc * keq
c_vis_s2 = beta_calc * k_st[0]
c_vis_s3 = beta_calc * k_st[1]
c_vis_s4 = beta_calc * k_st[2]

ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
for i in range(1, 6): ops.node(i, 0.0)
ops.fix(1, 1)
ops.mass(2, mass_iso); ops.mass(3, mass_st[0]); ops.mass(4, mass_st[1]); ops.mass(5, mass_st[2])

# Springs
ops.uniaxialMaterial('Elastic', 1, keq)
for i in range(3): ops.uniaxialMaterial('Elastic', i+2, k_st[i])
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.element('twoNodeLink', 2, 2, 3, '-mat', 2, '-dir', 1)
ops.element('twoNodeLink', 3, 3, 4, '-mat', 3, '-dir', 1)
ops.element('twoNodeLink', 4, 4, 5, '-mat', 4, '-dir', 1)

# Stiffness-proportional damping via Viscous elements (parallel to each spring)
ops.uniaxialMaterial('Viscous', 11, c_vis_iso, 1.0)
ops.element('twoNodeLink', 6, 1, 2, '-mat', 11, '-dir', 1)
ops.uniaxialMaterial('Viscous', 12, c_vis_s2, 1.0)
ops.element('twoNodeLink', 7, 2, 3, '-mat', 12, '-dir', 1)
ops.uniaxialMaterial('Viscous', 13, c_vis_s3, 1.0)
ops.element('twoNodeLink', 8, 3, 4, '-mat', 13, '-dir', 1)
ops.uniaxialMaterial('Viscous', 14, c_vis_s4, 1.0)
ops.element('twoNodeLink', 9, 4, 5, '-mat', 14, '-dir', 1)

# Mass-proportional Rayleigh ONLY
ops.rayleigh(alpha_calc, 0.0, 0.0, 0.0)

ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

dmax_os = 0; amax_os = 0
for i in range(n):
    ops.analyze(1, dt)
    d = ops.nodeDisp(5, 1); a = ops.nodeAccel(5, 1)
    dmax_os = max(dmax_os, abs(d)); amax_os = max(amax_os, abs(a))
ops.wipe()

print(f"OpenSees (fixed damp): d_max={dmax_os:.6f}, a_max={amax_os:.6f}")

err_d = abs(dmax_os - dmax_mat) / dmax_mat * 100
err_a = abs(amax_os - amax_mat) / amax_mat * 100
print(f"Displacement error: {err_d:.4f}%")
print(f"Acceleration error: {err_a:.4f}%")

if err_d < 0.1 and err_a < 0.1:
    print("PASS: Viscous + mass Rayleigh approach works for 4-DOF system!")
else:
    print("Still a discrepancy, further investigation needed.")
