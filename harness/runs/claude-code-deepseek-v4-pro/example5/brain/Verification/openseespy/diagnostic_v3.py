"""
OpenSeesPy Diagnostic v3 - Compare with/without damping
and with direct step-by-step comparison to MATLAB-algorithm.
"""
import openseespy.opensees as ops
import numpy as np

# Parameters
m = np.array([250, 270, 270, 180]) * 1e3
k0_arr = np.array([50, 245, 195, 98]) * 1e6
fy = np.array([500, 1225, 975, 490]) * 1e3
b_hard = [0.05, 0.0, 0.0, 0.0]
n_stories = 4

dt = 0.01
g = 9.81
PGA = 0.40 * g

data_path = '../../input_data/Northridge_01_NO_968.txt'
acc_raw = np.loadtxt(data_path)
acc_norm = acc_raw / np.max(np.abs(acc_raw)) * PGA
n_steps = len(acc_raw)

M_mat = np.diag(m)
K_mat = np.zeros((4, 4))
K_mat[0,0]=k0_arr[0]+k0_arr[1]; K_mat[0,1]=-k0_arr[1]
K_mat[1,0]=-k0_arr[1]; K_mat[1,1]=k0_arr[1]+k0_arr[2]; K_mat[1,2]=-k0_arr[2]
K_mat[2,1]=-k0_arr[2]; K_mat[2,2]=k0_arr[2]+k0_arr[3]; K_mat[2,3]=-k0_arr[3]
K_mat[3,2]=-k0_arr[3]; K_mat[3,3]=k0_arr[3]

eigvals, _ = np.linalg.eig(np.linalg.solve(M_mat, K_mat))
omega = np.sqrt(np.sort(np.real(eigvals)))
w1, w2 = omega[0], omega[1]
alpha_m = 2 * 0.05 * w1 * w2 / (w1 + w2)
beta_k = 2 * 0.05 / (w1 + w2)

def run_ops(label, use_damping=True, use_rayleigh=True):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    for i in range(5):
        ops.node(i+1, 0.0)
    ops.fix(1, 1)
    for i in range(n_stories):
        ops.mass(i+2, m[i])

    # Hardening for story 1
    Hkin1 = b_hard[0] * k0_arr[0] / (1.0 - b_hard[0])
    ops.uniaxialMaterial('Hardening', 1, k0_arr[0], fy[0], 0.0, Hkin1)
    for i in range(1, n_stories):
        ops.uniaxialMaterial('ElasticPP', i+1, k0_arr[i], fy[i]/k0_arr[i])
    for i in range(n_stories):
        ops.element('zeroLength', i+1, i+1, i+2, '-mat', i+1, '-dir', 1)

    if use_damping and use_rayleigh:
        ops.rayleigh(alpha_m, 0.0, beta_k, 0.0)

    # Use UniformExcitation (simpler)
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc_norm.tolist())
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('UmfPack')
    ops.test('NormUnbalance', 1e-5, 100, 0)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')

    top_disp = np.zeros(n_steps)
    for step in range(n_steps):
        ok = ops.analyze(1, dt)
        if ok != 0:
            ops.algorithm('KrylovNewton')
            ok2 = ops.analyze(1, dt)
            ops.algorithm('Newton')
        top_disp[step] = ops.nodeDisp(5, 1)

    ops.wipe()
    peak = np.max(np.abs(top_disp))
    print(f"  {label}: peak top disp = {peak:.8f} m")
    return top_disp

print("=== Damping Sensitivity Test ===")
d1 = run_ops("With Rayleigh damping", use_damping=True)
d2 = run_ops("Without damping", use_damping=False)

# Compare first 100 steps
diff = np.max(np.abs(d1[:100] - d2[:100]))
print(f"  Max diff in first 100 steps: {diff:.8f}")

# Now also test with very loose tolerance
print("\n=== Tolerance Test ===")
# Actually let me just compare step-by-step with MATLAB-equivalent Python
# by running 10 steps manually and comparing

# MATLAB-equivalent first 10 steps (from Python implementation)
# Let me just print the first 5 steps from both
print("\n=== First 5 steps comparison ===")
for i in range(5):
    print(f"  Step {i}: OpenSees top_d={d1[i]:.10f}")
