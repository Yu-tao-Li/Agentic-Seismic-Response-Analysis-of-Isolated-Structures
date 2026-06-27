"""
Three-story linear shear building - OpenSeesPy transient analysis.
Rayleigh damping (xi=0.05), GM1.txt normalized to 0.40g.
Uses 1D model with zeroLength shear springs.
"""
import numpy as np
import openseespy.opensees as ops
import json, os

# Parameters (SI)
nStories = 3
m_story  = 1000.0
k_story  = 500e3
xi       = 0.05
dt       = 0.02
g        = 9.81
targetPGA = 0.40 * g

# Load and normalize ground motion
gm_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'input_data', 'GM1.txt')
rawGM = np.loadtxt(gm_path)
rawPGA = np.max(np.abs(rawGM))
scale = targetPGA / rawPGA
ugddot = rawGM * scale
nSteps = len(ugddot)
print(f"N={nSteps}  rawPGA={rawPGA:.6f}  scale={scale:.6f}  PGA={np.max(np.abs(ugddot)):.6f}")

# Analytical eigenvalues
k_over_m = k_story / m_story
N = nStories
omega = np.array([2.0 * np.sqrt(k_over_m) * np.sin(np.pi/2 * (2*j-1)/(2*N+1)) for j in range(1, N+1)])
omega1, omega2, omega3 = omega[0], omega[1], omega[2]
alphaM = 2.0 * xi * omega1 * omega2 / (omega1 + omega2)
betaK  = 2.0 * xi / (omega1 + omega2)
print(f"omega1={omega1:.6f}  omega2={omega2:.6f}  omega3={omega3:.6f}")
print(f"alphaM={alphaM:.6f}  betaK={betaK:.8f}")

# Build OpenSees model (1D: ndm=1, ndf=1)
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Nodes: 1=base, 2=story1, 3=story2, 4=roof
for n in range(1, 5):
    ops.node(n, 0.0)

ops.fix(1, 1)

for n in range(2, 5):
    ops.mass(n, m_story)

# Zero-length shear springs
ops.uniaxialMaterial('Elastic', 1, k_story)
ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1, '-doRayleigh', 1)
ops.element('zeroLength', 2, 2, 3, '-mat', 1, '-dir', 1, '-doRayleigh', 1)
ops.element('zeroLength', 3, 3, 4, '-mat', 1, '-dir', 1, '-doRayleigh', 1)

# Rayleigh damping
ops.rayleigh(alphaM, betaK, 0.0, 0.0)

# Uniform excitation
ops.timeSeries('Path', 1, '-dt', dt, '-values', *ugddot.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

# Analysis setup
ops.constraints('Plain')
ops.numberer('Plain')
ops.system('UmfPack')
ops.test('NormUnbalance', 1.0e-12, 50, 0)
ops.algorithm('Linear')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

# Run
nAnalysis = nSteps - 1
time   = np.zeros(nAnalysis)
topDis = np.zeros(nAnalysis)
topAcc = np.zeros(nAnalysis)

for i in range(nAnalysis):
    ok = ops.analyze(1, dt)
    if ok != 0:
        print(f"Analysis failed at step {i}")
        break
    time[i]   = ops.getTime()
    topDis[i] = ops.nodeDisp(4, 1)
    topAcc[i] = ops.nodeAccel(4, 1)

print(f"Analysis complete. {nAnalysis} steps.")

# Write outputs
out_dir = os.path.dirname(os.path.abspath(__file__))
np.savetxt(os.path.join(out_dir, 'topStoDis.txt'), np.column_stack([time, topDis]), fmt='%.10e')
np.savetxt(os.path.join(out_dir, 'topStoAcc.txt'), np.column_stack([time, topAcc]), fmt='%.10e')

response_data = {
    "software": "OpenSeesPy",
    "time_step": dt,
    "n_samples": nAnalysis,
    "input_pga_normalized_m_s2": targetPGA,
    "raw_pga_m_s2": float(rawPGA),
    "scale_factor": float(scale),
    "top_displacement_peak_m": float(np.max(np.abs(topDis))),
    "top_acceleration_peak_m_s2": float(np.max(np.abs(topAcc))),
    "omega1_rad_s": float(omega1),
    "omega2_rad_s": float(omega2),
    "omega3_rad_s": float(omega3),
    "rayleigh_alpha": float(alphaM),
    "rayleigh_beta": float(betaK)
}

with open(os.path.join(out_dir, 'response_output.json'), 'w') as f:
    json.dump(response_data, f, indent=2)

print(f"Peak top displacement: {np.max(np.abs(topDis)):.6e} m")
print(f"Peak top acceleration: {np.max(np.abs(topAcc)):.6e} m/s^2")
ops.wipe()
