"""
OpenSeesPy single-DOF diagnostic test - compare Steel01 vs MATLAB bilinear.
Uses transient analysis with a slow sine wave (quasi-static).
"""
import openseespy.opensees as ops
import numpy as np

k0 = 50e6
fy = 500e3
b = 0.05
k_sh = b * k0  # 2.5e6

def matlab_force(disp):
    """MATLAB-style return mapping (single story)."""
    shear = np.zeros(len(disp))
    last_s = 0.0
    last_d = 0.0
    for i, d in enumerate(disp):
        delta_d = d - last_d
        c = last_s + k0 * delta_d
        fy_eff = fy * (1 - b)
        shear[i] = max(k_sh * d - fy_eff, min(k_sh * d + fy_eff, c))
        last_s = shear[i]
        last_d = d
    return shear

def run_opensees_transient(material_tag, mat_name, mat_args):
    """Transient analysis with imposed displacement via large-stiffness penalty."""
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.node(3, 0.0)  # reference node for displacement control
    ops.fix(1, 1)
    ops.fix(3, 1)

    ops.mass(2, 250000.0)

    # Material for zeroLength
    ops.uniaxialMaterial(mat_name, material_tag, *mat_args)
    ops.element('zeroLength', 1, 1, 2, '-mat', material_tag, '-dir', 1)

    # Displacement control via stiff elastic element
    ops.uniaxialMaterial('Elastic', 99, 1e15)  # very stiff
    ops.element('zeroLength', 99, 3, 2, '-mat', 99, '-dir', 1)

    # Apply displacement to node 3
    dt = 0.01
    t_final = 10.0
    n_steps = int(t_final / dt) + 1
    t = np.linspace(0, t_final, n_steps)

    freq = 1.0
    amp = 0.02  # enough to yield but not extreme
    ramp = np.minimum(t / 1.0, 1.0)
    disp = amp * np.sin(2 * np.pi * freq * t) * ramp

    ops.timeSeries('Path', 1, '-dt', dt, '-values', *disp.tolist())
    ops.pattern('Plain', 1, 1)
    ops.sp(3, 1, 1.0)

    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('UmfPack')
    ops.test('NormUnbalance', 1e-8, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', dt)
    ops.analysis('Static')

    forces = np.zeros(n_steps)
    for i in range(n_steps):
        ok = ops.analyze(1)
        if ok != 0:
            print(f"  Step {i} FAILED")
            forces[i] = forces[i-1] if i > 0 else 0.0
        else:
            forces[i] = -ops.eleResponse(1, 'force')[0]  # element 1 force
    ops.wipe()
    return forces, disp

# Run all material tests
print("=== Single-DOF Material Diagnostic ===\n")

# 1. Steel01
print("1. Steel01 (default)")
f_steel, d = run_opensees_transient(1, 'Steel01', [fy, k0, b])
f_matlab = matlab_force(d)
print(f"   Max |f|: {np.max(np.abs(f_steel)):.2f} N (MATLAB: {np.max(np.abs(f_matlab)):.2f})")
print(f"   Max diff: {np.max(np.abs(f_steel - f_matlab)):.2f} N ({np.max(np.abs(f_steel - f_matlab))/fy*100:.4f}% fy)")
print(f"   RMS diff: {np.sqrt(np.mean((f_steel - f_matlab)**2)):.2f} N")

# 2. Hardening
print("\n2. Hardening")
Hkin = b * k0 / (1.0 - b)  # plastic modulus
f_hard, d = run_opensees_transient(2, 'Hardening', [k0, fy, 0.0, Hkin])
f_matlab = matlab_force(d)
print(f"   Max |f|: {np.max(np.abs(f_hard)):.2f} N (MATLAB: {np.max(np.abs(f_matlab)):.2f})")
print(f"   Max diff: {np.max(np.abs(f_hard - f_matlab)):.2f} N ({np.max(np.abs(f_hard - f_matlab))/fy*100:.4f}% fy)")
print(f"   RMS diff: {np.sqrt(np.mean((f_hard - f_matlab)**2)):.2f} N")

# 3. ElasticPP (b=0 case for upper stories)
print("\n3. ElasticPP (for b=0 upper stories)")
f_pp, d = run_opensees_transient(3, 'ElasticPP', [k0, fy/k0])
f_matlab = matlab_force(d)
print(f"   Max |f|: {np.max(np.abs(f_pp)):.2f} N (MATLAB: {np.max(np.abs(f_matlab)):.2f})")
print(f"   Max diff: {np.max(np.abs(f_pp - f_matlab)):.2f} N ({np.max(np.abs(f_pp - f_matlab))/fy*100:.4f}% fy)")
print(f"   RMS diff: {np.sqrt(np.mean((f_pp - f_matlab)**2)):.2f} N")

# Now test EPP (b=0) case more carefully
print("\n=== EPP (b=0) test ===")
def matlab_force_epp(disp):
    shear = np.zeros(len(disp))
    last_s = 0.0
    last_d = 0.0
    for i, d in enumerate(disp):
        delta_d = d - last_d
        c = last_s + k0 * delta_d
        shear[i] = max(-fy, min(fy, c))  # clip to [-fy, fy]
        last_s = shear[i]
        last_d = d
    return shear

f_epp, d = run_opensees_transient(4, 'ElasticPP', [k0, fy/k0])
f_matlab_epp = matlab_force_epp(d)
print(f"   Max |f|: {np.max(np.abs(f_epp)):.2f} N (MATLAB: {np.max(np.abs(f_matlab_epp)):.2f})")
print(f"   Max diff: {np.max(np.abs(f_epp - f_matlab_epp)):.2f} N")

# Summary
print("\n=== SUMMARY ===")
print("Hardening material should match MATLAB bilinear kinematic hardening exactly.")
print("ElasticPP with b=0 should match MATLAB EPP exactly.")
print("Steel01 has rounded yield transition, expected to differ slightly near yield.")
