"""
Generate ABAQUS input file for Example 5.

The MATLAB/OpenSeesPy reference uses bilinear kinematic hysteresis for each
story spring and Rayleigh damping C = alpha*M + beta*K0.  Connector plasticity
in Abaqus gave a different cyclic rule for the isolation layer, so this model
uses axial T2D2 truss elements with 1D kinematic plasticity for the restoring
springs and explicit DASHPOT2 elements for beta*K0 damping.
"""
import os
import numpy as np

damping_mode = os.environ.get('ABAQUS_DAMPING_MODE', 'dashpot').lower()
if damping_mode not in ('dashpot', 'global_beta'):
    raise ValueError('ABAQUS_DAMPING_MODE must be dashpot or global_beta')
alpha_damping_mode = os.environ.get('ABAQUS_ALPHA_DAMPING_MODE', 'explicit_dashpot').lower()
if alpha_damping_mode not in ('explicit_dashpot', 'global'):
    raise ValueError('ABAQUS_ALPHA_DAMPING_MODE must be explicit_dashpot or global')
plasticity_mode = os.environ.get('ABAQUS_PLASTICITY_MODE', 'combined_parameters').lower()
if plasticity_mode not in ('kinematic_table', 'combined_parameters'):
    raise ValueError('ABAQUS_PLASTICITY_MODE must be kinematic_table or combined_parameters')

m = np.array([250, 270, 270, 180]) * 1e3
k = np.array([50, 245, 195, 98]) * 1e6
fy = np.array([500, 1225, 975, 490]) * 1e3
b_hard = np.array([0.05, 0.0, 0.0, 0.0])
dt = 0.01
substeps_per_dt = int(os.environ.get('ABAQUS_SUBSTEPS_PER_DT', '1'))
if substeps_per_dt < 1:
    raise ValueError('ABAQUS_SUBSTEPS_PER_DT must be a positive integer')
analysis_dt = dt / substeps_per_dt
g = 9.81
PGA = 0.40 * g
n_stories = 4

acc_raw = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc_norm = acc_raw / np.max(np.abs(acc_raw)) * PGA
n_steps = len(acc_norm)
total_time = (n_steps - 1) * dt

M_mat = np.diag(m)
K0_mat = np.zeros((4, 4))
K0_mat[0, 0] = k[0] + k[1]
K0_mat[0, 1] = -k[1]
K0_mat[1, 0] = -k[1]
K0_mat[1, 1] = k[1] + k[2]
K0_mat[1, 2] = -k[2]
K0_mat[2, 1] = -k[2]
K0_mat[2, 2] = k[2] + k[3]
K0_mat[2, 3] = -k[3]
K0_mat[3, 2] = -k[3]
K0_mat[3, 3] = k[3]
eigvals, _ = np.linalg.eig(np.linalg.solve(M_mat, K0_mat))
omega = np.sqrt(np.sort(np.real(eigvals)))
w1, w2 = omega[0], omega[1]
alpha_m = 2 * 0.05 * w1 * w2 / (w1 + w2)
beta_k = 2 * 0.05 / (w1 + w2)
alpha_scale = float(os.environ.get('ABAQUS_ALPHA_SCALE', '1.0'))
beta_scale = float(os.environ.get('ABAQUS_BETA_SCALE', '1.0'))
alpha_m *= alpha_scale
beta_k *= beta_scale
print(f'Damping: alpha={alpha_m:.6f}, beta={beta_k:.6f}')

lines = []


def L(s=''):
    lines.append(s)


L('*HEADING')
L('Example 5: 4-story isolated shear building, truss plasticity')
L('*PREPRINT,ECHO=NO,MODEL=NO,HISTORY=NO')

# Nodes are spaced one meter apart in X so axial truss strain equals relative
# story displacement. NLGEOM=NO keeps this a small-displacement shear model.
L('*NODE')
for i in range(5):
    L(f'{i + 1}, {float(i):.1f}, 0.0')
L('*NSET,NSET=BASE')
L('1')
L('*NSET,NSET=TOP')
L('5')

# Nonlinear restoring springs.
for i in range(n_stories):
    L(f'*ELEMENT,TYPE=T2D2,ELSET=SPRING{i + 1}')
    L(f'{i + 1}, {i + 1}, {i + 2}')

if damping_mode == 'dashpot':
    # Linear dashpots for the stiffness-proportional Rayleigh term beta*K0.
    for i in range(n_stories):
        eid = 201 + i
        L(f'*ELEMENT,TYPE=DASHPOT2,ELSET=DASHPOT{i + 1}')
        L(f'{eid}, {i + 1}, {i + 2}')
        L(f'*DASHPOT,ELSET=DASHPOT{i + 1}')
        L('1, 1')
        L(f'{beta_k * k[i]:.10e},')

if alpha_damping_mode == 'explicit_dashpot':
    # Abaqus global alpha damping did not affect this point-mass/truss model
    # consistently, so alpha*M is represented by nodal dashpots to the fixed
    # base in the relative-coordinate equation of motion.
    for i in range(n_stories):
        eid = 301 + i
        node = i + 2
        L(f'*ELEMENT,TYPE=DASHPOT2,ELSET=ALPHA_DASHPOT{i + 1}')
        L(f'{eid}, 1, {node}')
        L(f'*DASHPOT,ELSET=ALPHA_DASHPOT{i + 1}')
        L('1, 1')
        L(f'{alpha_m * m[i]:.10e},')

# Lumped masses at floor nodes.
for i in range(n_stories):
    L(f'*ELEMENT,TYPE=MASS,ELSET=MASS{i + 1}')
    L(f'{101 + i}, {i + 2}')
    L(f'*MASS,ELSET=MASS{i + 1}')
    L(f'{m[i]:.1f},')

# Material mapping: A=L=1, so force = stress and stiffness = E.
for i in range(n_stories):
    L(f'*MATERIAL,NAME=MAT{i + 1}')
    L('*DENSITY')
    L('1.0e-12,')
    L('*ELASTIC')
    L(f'{k[i]:.10e}, 0.0')
    if b_hard[i] > 0:
        H = b_hard[i] * k[i] / (1.0 - b_hard[i])
    else:
        H = 0.0
    if plasticity_mode == 'combined_parameters':
        L('*PLASTIC,HARDENING=COMBINED,DATA TYPE=PARAMETERS,NUMBER BACKSTRESSES=1')
        L(f'{fy[i]:.10e}, {H:.10e}, 0.0')
    else:
        L('*PLASTIC,HARDENING=KINEMATIC')
        L(f'{fy[i]:.10e}, 0.0')
        L(f'{fy[i] + H:.10e}, 1.0')
    L(f'*SOLID SECTION,ELSET=SPRING{i + 1},MATERIAL=MAT{i + 1}')
    L('1.0,')

# Boundary conditions: only horizontal X motion is active at floor nodes.
L('*BOUNDARY')
L('BASE, 1, 2')
for node in range(2, 6):
    L(f'{node}, 2, 2')

acc_unit = acc_raw / np.max(np.abs(acc_raw))
L('*AMPLITUDE,NAME=GM,TIME=TOTAL TIME')
for i, a in enumerate(acc_unit):
    L(f'{i * dt:.6f}, {a:.10e}')

max_increments = (n_steps - 1) * substeps_per_dt + 20
L(f'*STEP,NAME=DYNAMIC,NLGEOM=NO,INC={max_increments}')
L('*DYNAMIC,ALPHA=0.0,DIRECT')
L(f'{analysis_dt:.6f},{total_time:.6f},{analysis_dt / 100:.6f},{analysis_dt:.6f}')

if damping_mode == 'dashpot':
    # beta*K0 is the explicit story dashpot chain above.  alpha*M is either
    # represented explicitly by nodal dashpots or delegated to Abaqus global
    # damping for diagnostic control runs.
    alpha_global = alpha_m if alpha_damping_mode == 'global' else 0.0
    L(f'*GLOBAL DAMPING, ALPHA={alpha_global:.10e}, BETA=0.0')
else:
    # Control run: let Abaqus apply stiffness-proportional damping globally.
    L(f'*GLOBAL DAMPING, ALPHA={alpha_m:.10e}, BETA={beta_k:.10e}')

L('*CLOAD, AMPLITUDE=GM')
for i in range(n_stories):
    L(f'{i + 2}, 1, {-m[i] * PGA:.10e}')

L('*OUTPUT, FIELD, FREQUENCY=1')
L('*NODE OUTPUT')
L('U, V, A, RF')
for i in range(n_stories):
    L(f'*ELEMENT OUTPUT, ELSET=SPRING{i + 1}')
    L('S, E, PE, PEEQ')
L('*END STEP')

with open('example5_model.inp', 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(
    f'ABAQUS input file: example5_model.inp '
    f'({n_steps} input points, analysis_dt={analysis_dt}, damping={damping_mode}, '
    f'plasticity={plasticity_mode})'
)
