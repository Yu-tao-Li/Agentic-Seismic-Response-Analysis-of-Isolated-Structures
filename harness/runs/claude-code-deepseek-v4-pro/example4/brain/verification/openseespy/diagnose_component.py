#!/usr/bin/env python
"""Isolate which Rayleigh component is broken in OpenSees."""
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

def run_os(description, am, bk, bki, bkc):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.mass(2, m1)
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    if am != 0 or bk != 0 or bki != 0 or bkc != 0:
        ops.rayleigh(am, bk, bki, bkc)
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

    # Compute effective C from first step
    K_eff_first = F0 / d[0]
    C_eff = (K_eff_first - k1 - c0*m1) / c1
    print(f"{description}:")
    print(f"  u(1)={d[0]:.8f}, K_eff={K_eff_first:.2f}, C_eff={C_eff:.4f}")
    return d, C_eff

# Reference: MATLAB-style with full damping
C_full = alpha*m1 + beta*k1  # = 2.0
K_eff_full = k1 + c1*C_full + c0*m1
u1_full = F0 / K_eff_full
print(f"MATLAB full damping: C={C_full:.1f}, K_eff={K_eff_full:.2f}, u(1)={u1_full:.8f}\n")

# Test 1: Mass-proportional only
run_os("OS mass-only  ", alpha, 0.0, 0.0, 0.0)

# Test 2: Stiffness-proportional (current tangent)
run_os("OS stiff-curr ", 0.0, beta, 0.0, 0.0)

# Test 3: Stiffness-proportional (initial)
run_os("OS stiff-init ", 0.0, 0.0, beta, 0.0)

# Test 4: Stiffness-proportional (committed)
run_os("OS stiff-comm ", 0.0, 0.0, 0.0, beta)

# Test 5: Both mass + stiff-init
run_os("OS mass+stiff ", alpha, 0.0, beta, 0.0)

# Test 6: Both mass + stiff-curr
run_os("OS mass+curr  ", alpha, beta, 0.0, 0.0)

# Test 7: No rayleigh (baseline)
run_os("OS no damping ", 0.0, 0.0, 0.0, 0.0)

# Also test with twoNodeLink instead of zeroLength
print("\n=== twoNodeLink test ===")
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)
ops.node(1, 0.0); ops.node(2, 0.0)
ops.fix(1, 1)
ops.mass(2, m1)
ops.uniaxialMaterial('Elastic', 1, k1)
ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
ops.rayleigh(alpha, 0.0, beta, 0.0)
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, F0)
ops.constraints('Transformation'); ops.numberer('RCM'); ops.system('FullGeneral')
ops.test('NormDispIncr', 1e-12, 50, 0); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

d_tnl = np.zeros(n_steps)
for i in range(n_steps):
    ops.analyze(1, dt)
    d_tnl[i] = ops.nodeDisp(2, 1)
ops.wipe()

K_eff_tnl = F0 / d_tnl[0]
C_tnl = (K_eff_tnl - k1 - c0*m1) / c1
print(f"twoNodeLink mass+stiff: u(1)={d_tnl[0]:.8f}, K_eff={K_eff_tnl:.2f}, C_eff={C_tnl:.4f}")

# Try different rayleigh argument combinations
print("\n=== Alternative rayleigh args ===")
for desc, args in [
    ("am, bk, bki=0, bkc=0      ", (alpha, beta, 0.0, 0.0)),
    ("am, bk=0, bki, bkc=0       ", (alpha, 0.0, beta, 0.0)),
    ("am=0, bk, bki, bkc=0       ", (0.0, beta, beta, 0.0)),
    ("am=0, bk=0, bki, bkc       ", (0.0, 0.0, beta, beta)),
    ("am, bk=0, bki=0, bkc       ", (alpha, 0.0, 0.0, beta)),
    ("am, bk, bki, bkc           ", (alpha, beta, beta, beta)),
]:
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.mass(2, m1)
    ops.uniaxialMaterial('Elastic', 1, k1)
    ops.element('twoNodeLink', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.rayleigh(*args)
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
