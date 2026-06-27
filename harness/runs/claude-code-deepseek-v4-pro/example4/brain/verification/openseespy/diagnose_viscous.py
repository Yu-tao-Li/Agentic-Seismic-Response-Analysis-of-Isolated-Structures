#!/usr/bin/env python
"""Test if Viscous materials can substitute for stiffness-proportional Rayleigh damping."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import openseespy.opensees as ops

m1 = 1.0; k1 = 100.0
omega_n = np.sqrt(k1/m1)
alpha = 2*0.05*omega_n      # = 1.0
beta = 2*0.05/omega_n        # = 0.01
dt = 0.01; F0 = 100.0
n_steps = 20

beta_nm = 0.25; gamma_nm = 0.5
c0 = 1/(beta_nm*dt**2); c1 = gamma_nm/(beta_nm*dt); c2 = 1/(beta_nm*dt)
c3 = 1/(2*beta_nm)-1; c4 = gamma_nm/beta_nm-1; c5 = dt*(gamma_nm/(2*beta_nm)-1)

def test_run(desc, setup_fn):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.node(2, 0.0)
    ops.fix(1, 1)
    setup_fn()
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, F0)
    ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

    d = np.zeros(n_steps)
    for i in range(n_steps):
        ops.analyze(1, dt)
        d[i] = ops.nodeDisp(2, 1)
    ops.wipe()

    K_eff = F0 / d[0]
    C_eff = (K_eff - k1 - c0*m1) / c1
    print(f"  {desc}: u(1)={d[0]:.8f}, C_eff={C_eff:.4f}")
    return d

# Reference: MATLAB full damping (C=2.0)
C_full = alpha*m1 + beta*k1
K_eff_full = k1 + c1*C_full + c0*m1
u_ref = F0 / K_eff_full
print(f"Reference (C=2.0): u(1)={u_ref:.8f}\n")

# Test 1: mass rayleigh only (C_mass=1.0)
def setup_1():
    ops.mass(2, m1)
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.rayleigh(alpha, 0.0, 0.0, 0.0)
test_run("Mass Rayleigh only     ", setup_1)

# Test 2: mass rayleigh + viscous element in parallel (C_stiff via viscous)
def setup_2():
    c_stiff_eqv = beta * k1  # = 1.0
    ops.mass(2, m1)
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.uniaxialMaterial('Viscous', 10, c_stiff_eqv, 1.0)
    ops.element('zeroLength', 2, 1, 2, '-mat', 10, '-dir', 1)
    ops.rayleigh(alpha, 0.0, 0.0, 0.0)  # mass only
test_run("Mass Rayleigh + Viscous", setup_2)

# Test 3: full damping via mass rayleigh + Viscous for stiff part
# Use twoNodeLink elements instead
def setup_3():
    c_stiff_eqv = beta * k1
    ops.mass(2, m1)
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.uniaxialMaterial('Viscous', 10, c_stiff_eqv, 1.0)
    ops.element('twoNodeLink', 2, 1, 2, '-mat', 10, '-dir', 1)
    ops.rayleigh(alpha, 0.0, 0.0, 0.0)  # mass only
test_run("TNL: Mass Ray + Vis  ", setup_3)

# Test 4: Only Viscous, no mass Rayleigh
def setup_4():
    c_stiff_eqv = beta * k1
    ops.mass(2, m1)
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.uniaxialMaterial('Viscous', 10, c_stiff_eqv, 1.0)
    ops.element('twoNodeLink', 2, 1, 2, '-mat', 10, '-dir', 1)
    # NO rayleigh at all
test_run("TNL: Viscous only     ", setup_4)

# Test 5: Use Node damping instead of Rayleigh mass-proportional
# In OpenSees, we can set damping at nodes
def setup_5():
    c_stiff_eqv = beta * k1  # = 1.0
    c_mass_eqv = alpha * m1  # = 1.0
    ops.mass(2, m1)
    # Combined damping in one Viscous element: C = alpha*m1 + beta*k1 = 2.0
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.uniaxialMaterial('Viscous', 10, c_stiff_eqv + c_mass_eqv, 1.0)
    ops.element('twoNodeLink', 2, 1, 2, '-mat', 10, '-dir', 1)
    # NO rayleigh
test_run("TNL: Vis(C=2) no Rayl", setup_5)

# Test 6: Full equivalent with viscous, full time history comparison
print("\n=== Full 1-DOF time history comparison ===")
# Python reference
d_ref = [0.0]; v_ref = [0.0]; a_ref = [0.0]
C = 2.0
K_eff_py = k1 + c1*C + c0*m1
for i in range(n_steps):
    u_n, v_n, a_n = d_ref[-1], v_ref[-1], a_ref[-1]
    F_eff = F0 + m1*(c0*u_n + c2*v_n + c3*a_n) + C*(c1*u_n + c4*v_n + c5*a_n)
    u_new = F_eff / K_eff_py
    a_new = c0*(u_new - u_n) - c2*v_n - c3*a_n
    v_new = v_n + dt*((1-gamma_nm)*a_n + gamma_nm*a_new)
    d_ref.append(u_new); v_ref.append(v_new); a_ref.append(a_new)
d_ref = d_ref[1:]

# OpenSees: Viscous(C=2.0) for full damping
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
ops.node(1, 0.0); ops.node(2, 0.0)
ops.fix(1, 1)
ops.mass(2, m1)
ops.uniaxialMaterial('Elastic', 1, k1)
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.uniaxialMaterial('Viscous', 10, 2.0, 1.0)  # C=2.0
ops.element('twoNodeLink', 2, 1, 2, '-mat', 10, '-dir', 1)
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, F0)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

d_os = np.zeros(n_steps)
for i in range(n_steps):
    ops.analyze(1, dt)
    d_os[i] = ops.nodeDisp(2, 1)
ops.wipe()

max_err = np.max(np.abs(d_os - np.array(d_ref)))
print(f"Max disp error (Viscous C=2 vs Ref): {max_err:.8f}")
if max_err < 1e-6:
    print("PASS: Viscous elements correctly implement damping!")
else:
    # Check how big error actually is
    max_d = np.max(np.abs(d_ref))
    print(f"Relative error: {max_err/max_d*100:.4f}%")
    # Maybe Viscous doesn't contribute to K_eff. Check.
    K_eff_os = F0 / d_os[0]
    C_eff_os = (K_eff_os - k1 - c0*m1) / c1
    print(f"K_eff from OpenSees: {K_eff_os:.2f}, C_eff: {C_eff_os:.4f}")

# Check what happens if I use ONLY rayleigh mass-proportional + viscous for stiff
# This is the most practical approach
print("\n=== Approach: mass Rayleigh + Viscous per element ===")
# For the 4-DOF system, this would be:
# - ops.rayleigh(alpha_m, 0.0, 0.0, 0.0) for mass-proportional
# - For each element i: Viscous material with C = beta_k * k_i, parallel to the spring
print("This approach should give correct total damping.")
