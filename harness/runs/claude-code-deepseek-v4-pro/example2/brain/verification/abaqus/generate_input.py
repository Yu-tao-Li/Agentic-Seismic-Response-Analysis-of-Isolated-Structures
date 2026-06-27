"""
Generate ABAQUS INP and amplitude files for three-story shear building.
Rayleigh damping via DASHPOT1 (alpha*M) + CONNECTOR DAMPING (beta*K).
*A(0)=0 in amplitude to match MATLAB initial acceleration convention.
"""
import numpy as np
import os

nStories = 3
m_story  = 1000.0
k_story  = 500e3
xi       = 0.05
dt       = 0.02
g        = 9.81
targetPGA = 0.40 * g

script_dir = os.path.dirname(os.path.abspath(__file__))
gm_path = os.path.join(script_dir, '..', '..', 'input_data', 'GM1.txt')
rawGM = np.loadtxt(gm_path)
rawPGA = np.max(np.abs(rawGM))
scale = targetPGA / rawPGA
ugddot = rawGM * scale
nSteps = len(ugddot)
total_time = (nSteps - 1) * dt

print(f"N={nSteps}  rawPGA={rawPGA:.6f}  scale={scale:.6f}  PGA={targetPGA:.6f}")
print(f"Total time: {total_time:.2f}s  dt={dt:.2f}s")

k_over_m = k_story / m_story
N = nStories
omega = np.array([2.0 * np.sqrt(k_over_m) * np.sin(np.pi/2 * (2*j-1)/(2*N+1)) for j in range(1, N+1)])
omega1, omega2 = omega[0], omega[1]
alphaM = 2.0 * xi * omega1 * omega2 / (omega1 + omega2)
betaK  = 2.0 * xi / (omega1 + omega2)
c_alpha = alphaM * m_story    # dashpot coefficient for mass-proportional damping
c_beta  = betaK * k_story     # connector damping coefficient for stiffness-proportional damping

print(f"omega1={omega1:.6f}  omega2={omega2:.6f}  omega3={omega[2]:.6f}")
print(f"alphaM={alphaM:.6f}  betaK={betaK:.8f}")
print(f"DASHPOT1 coefficient (alpha*m): {c_alpha:.3f} N*s/m")
print(f"CONNECTOR DAMPING (beta*k): {c_beta:.3f} N*s/m")

# Amplitude file: first value=0 so that initial acceleration a_0=0 (matches MATLAB convention)
ugddot_adj = ugddot.copy()
ugddot_adj[0] = 0.0

amp_file = os.path.join(script_dir, 'eq_amplitude.txt')
np.savetxt(amp_file, np.column_stack([np.arange(nSteps)*dt, ugddot_adj]),
           fmt=['%.6f', '%.12e'], delimiter=', ', header='')
print(f"Amplitude written: {amp_file} ({nSteps} rows, first=0)")

inp_path = os.path.join(script_dir, 'three_story_shear.inp')
with open(inp_path, 'w') as f:
    f.write('*HEADING\n')
    f.write(f'Three-story shear building, DASHPOT1+CONNECTOR DAMPING, dt={dt}, PGA={targetPGA:.3f}m/s2\n')
    f.write('*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO\n')
    f.write('**\n')
    f.write('*NODE\n')
    f.write(' 1, 0.0, 0.0, 0.0\n')
    f.write(' 2, 1.0, 0.0, 0.0\n')
    f.write(' 3, 2.0, 0.0, 0.0\n')
    f.write(' 4, 3.0, 0.0, 0.0\n')
    f.write('**\n')
    f.write('** Connector elements with stiffness + beta*K damping\n')
    f.write('*ELEMENT, TYPE=CONN3D2, ELSET=CONN\n')
    f.write('1, 1, 2\n')
    f.write('2, 2, 3\n')
    f.write('3, 3, 4\n')
    f.write('*CONNECTOR BEHAVIOR, NAME=CONN_BEHAVIOR\n')
    f.write('*CONNECTOR ELASTICITY, COMPONENT=1\n')
    f.write(f'{k_story:.1f},\n')
    f.write('*CONNECTOR DAMPING, COMPONENT=1\n')
    f.write(f'{c_beta:.3f},\n')
    f.write('*CONNECTOR SECTION, ELSET=CONN, BEHAVIOR=CONN_BEHAVIOR\n')
    f.write(' CARTESIAN,\n')
    f.write(' 1, 0, 0\n')
    f.write('**\n')
    f.write('** Dashpot elements for alpha*M damping (to ground)\n')
    f.write('*ELEMENT, TYPE=DASHPOT1, ELSET=DASH\n')
    f.write('11, 2\n')
    f.write('12, 3\n')
    f.write('13, 4\n')
    f.write('*DASHPOT, ELSET=DASH\n')
    f.write('1\n')
    f.write(f'{c_alpha:.3f},\n')
    f.write('**\n')
    f.write('** Lumped masses\n')
    f.write('*ELEMENT, TYPE=MASS, ELSET=MAS\n')
    f.write('101, 2\n')
    f.write('102, 3\n')
    f.write('103, 4\n')
    f.write('*MASS, ELSET=MAS\n')
    f.write(f'{m_story:.1f},\n')
    f.write('**\n')
    f.write('*NSET, NSET=ROOF\n')
    f.write('4\n')
    f.write('**\n')
    f.write('*BOUNDARY\n')
    f.write('1, 1, 3\n')
    f.write('2, 2, 3\n')
    f.write('3, 2, 3\n')
    f.write('4, 2, 3\n')
    f.write('**\n')
    f.write('*AMPLITUDE, NAME=EQ_AMP, INPUT=eq_amplitude.txt\n')
    f.write('**\n')
    f.write(f'*STEP, NAME=DYNAMIC, INC={nSteps}\n')
    f.write('*DYNAMIC, ALPHA=0.0, DIRECT\n')
    f.write(f'{dt}, {total_time:.4f}, {dt}, {dt}\n')
    f.write('**\n')
    f.write('*CLOAD, AMPLITUDE=EQ_AMP\n')
    f.write(f'2, 1, {-m_story:.1f}\n')
    f.write(f'3, 1, {-m_story:.1f}\n')
    f.write(f'4, 1, {-m_story:.1f}\n')
    f.write('**\n')
    f.write('*OUTPUT, FIELD, NUMBER INTERVAL=53\n')
    f.write('*NODE OUTPUT\n')
    f.write('U, A\n')
    f.write('*OUTPUT, HISTORY, FREQUENCY=1\n')
    f.write('*NODE OUTPUT, NSET=ROOF\n')
    f.write('U1, A1\n')
    f.write('*END STEP\n')

print(f"INP written: {inp_path}")
print("Run: abaqus job=three_story_shear input=three_story_shear.inp")
print("Extract: abaqus python extract_results.py")
