"""
Generate ABAQUS input file for Example 5:
4-story seismically isolated structure.
CONN3D2 connector elements for nonlinear shear springs.
Ground motion via CLOAD on mass nodes.
Rayleigh damping via CONNECTOR DAMPING (beta*K) + DASHPOT1 (alpha*M).
"""
import numpy as np

m = np.array([250, 270, 270, 180]) * 1e3
k = np.array([50, 245, 195, 98]) * 1e6
fy = np.array([500, 1225, 975, 490]) * 1e3
b_hard = [0.05, 0.0, 0.0, 0.0]
dt = 0.01
g = 9.81
PGA = 0.40 * g
n_stories = 4

acc_raw = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc_norm = acc_raw / np.max(np.abs(acc_raw)) * PGA
n_steps = len(acc_norm)
total_time = (n_steps - 1) * dt

M_mat = np.diag(m)
K0_mat = np.zeros((4, 4))
K0_mat[0,0]=k[0]+k[1]; K0_mat[0,1]=-k[1]
K0_mat[1,0]=-k[1]; K0_mat[1,1]=k[1]+k[2]; K0_mat[1,2]=-k[2]
K0_mat[2,1]=-k[2]; K0_mat[2,2]=k[2]+k[3]; K0_mat[2,3]=-k[3]
K0_mat[3,2]=-k[3]; K0_mat[3,3]=k[3]
eigvals, _ = np.linalg.eig(np.linalg.solve(M_mat, K0_mat))
omega = np.sqrt(np.sort(np.real(eigvals)))
w1, w2 = omega[0], omega[1]
alpha_m = 2 * 0.05 * w1 * w2 / (w1 + w2)
beta_k = 2 * 0.05 / (w1 + w2)
print(f'Damping: alpha={alpha_m:.6f}, beta={beta_k:.6f}')

lines = []
def L(s=''):
    lines.append(s)

L('*HEADING')
L('Example 5: 4-story seismically isolated structure')
L('*PREPRINT,ECHO=NO,MODEL=NO,HISTORY=NO')

# Nodes (co-located for shear building)
L('*NODE')
for i in range(5):
    L(f'{i+1}, 0.0, 0.0, 0.0')
L('*NSET,NSET=BASE')
L('1')
L('*NSET,NSET=TOP')
L('5')

# CONN3D2 elements
for i in range(4):
    L(f'*ELEMENT,TYPE=CONN3D2,ELSET=E{i+1}')
    L(f'{i+1}, {i+1}, {i+2}')

# Mass elements
for i in range(4):
    L(f'*ELEMENT,TYPE=MASS,ELSET=MASS{i+1}')
    L(f'{100+i+1}, {i+2}')
    L(f'*MASS,ELSET=MASS{i+1}')
    L(f'{m[i]:.1f},')

# Combined mass elset for loading
L('*ELSET,ELSET=MASSEL')
L('MASS1,MASS2,MASS3,MASS4')

# Mass-proportional damping: DASHPOT1 elements (alpha*m_i to ground in X-dir)
for i in range(4):
    L(f'*ELEMENT,TYPE=DASHPOT1,ELSET=DAMP_M{i+1}')
    L(f'{200+i+1}, {i+2}')
    L(f'*DASHPOT,ELSET=DAMP_M{i+1}')
    L('1')
    L(f'{alpha_m * m[i]:.6f}')

# Connector behaviors and sections
for i in range(4):
    L(f'*CONNECTOR BEHAVIOR,NAME=BHV{i+1}')
    L('*CONNECTOR ELASTICITY,COMPONENT=1')
    L(f'{k[i]:.1f},')
    L('*CONNECTOR PLASTICITY,COMPONENT=1')
    L(f'{fy[i]:.1f},')
    # Kinematic hardening (required by ABAQUS when plasticity is defined)
    if b_hard[i] > 1e-6:
        Hkin = b_hard[i] * k[i] / (1.0 - b_hard[i])
        upl_ref = 0.05
        L('*CONNECTOR HARDENING, TYPE=KINEMATIC')
        L(f'{fy[i]:.1f}, 0.0')
        L(f'{fy[i] + Hkin * upl_ref:.1f}, {upl_ref:.6f}')
    else:
        L('*CONNECTOR HARDENING, TYPE=KINEMATIC')
        L(f'{fy[i]:.1f}, 0.0')
        L(f'{fy[i]:.1f}, 1.0')
    # Stiffness-proportional damping: connector dashpot C_i = beta * k_i
    L('*CONNECTOR DAMPING, COMPONENT=1')
    L(f'{beta_k * k[i]:.6f},')
    # Components 2,3 (Y,Z): rigid
    for comp in [2, 3]:
        L(f'*CONNECTOR ELASTICITY,COMPONENT={comp}')
        L('1.0E12,')
    L(f'*CONNECTOR SECTION,ELSET=E{i+1},BEHAVIOR=BHV{i+1}')
    L('CARTESIAN,')

# Boundary conditions (base fixed)
L('*BOUNDARY')
L('BASE, 1, 3')

# Ground motion amplitude: unit-normalized acceleration (-1 to 1)
acc_unit = acc_raw / np.max(np.abs(acc_raw))
L('*AMPLITUDE,NAME=GM,TIME=TOTAL TIME')
for i, a in enumerate(acc_unit):
    t = i * dt
    L(f'{t:.6f}, {a:.6e}')

# Step: implicit dynamic analysis (pure Newmark, ALPHA=0)
L(f'*STEP,NAME=DYNAMIC,NLGEOM=NO,INC={n_steps+10}')
L('*DYNAMIC,ALPHA=0.0,DIRECT')
L(f'{dt:.6f},{total_time:.6f},{dt/100:.6f},{dt:.6f}')

# Ground motion: concentrated forces at each mass node
L('*CLOAD, AMPLITUDE=GM')
for i in range(4):
    peak_force = -m[i] * PGA
    L(f'{i+2}, 1, {peak_force:.6f}')

# Output at every increment
L('*OUTPUT, FIELD, FREQUENCY=1')
L('*NODE OUTPUT')
L('U,')
L('A,')
L('V,')
L('*ELEMENT OUTPUT')
L('CTF,')
L('*END STEP')

with open('example5_model.inp', 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f'ABAQUS input file: example5_model.inp ({n_steps} steps)')
