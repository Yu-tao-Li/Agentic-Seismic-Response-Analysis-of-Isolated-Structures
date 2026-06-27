"""
2D Three-Story Two-Bay Base-Isolated Steel Frame
OpenSeesPy Analysis with Bouc-Wen Isolators
Harness: example6 | Software: OpenSeesPy
"""
import openseespy.opensees as ops
import numpy as np
import json
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# =============================================================================
# 1. MATERIAL AND SECTION PROPERTIES
# =============================================================================
E   = 2.06e11      # Elastic modulus [Pa]
nu  = 0.30         # Poisson ratio
rho_steel = 7850   # Density [kg/m^3]
A_c = 0.020        # Column area [m^2]
I_c = 8.0e-4       # Column 2nd moment of area [m^4]
A_b = 0.015        # Beam area [m^2]
I_b = 4.5e-4       # Beam 2nd moment of area [m^4]

# =============================================================================
# 2. GEOMETRY
# =============================================================================
n_stories = 3
n_bays = 2
H_story = 3.6       # Story height [m]
L_bay = 6.0         # Bay width [m]

# =============================================================================
# 3. ISOLATOR PARAMETERS (Bouc-Wen)
# =============================================================================
k0    = 2.0e6       # Initial horizontal stiffness [N/m]
Fy    = 1.0e4       # Yield force [N]
uy    = Fy / k0     # Yield deformation [m]
alpha_iso = 0.05    # Post-yield stiffness ratio
n_bw  = 2.0         # Smoothness exponent
gamma_bw = 2.0e4    # Bouc-Wen shape parameter
beta_bw  = 2.0e4    # Bouc-Wen shape parameter
Ao    = 1.0         # Reference amplitude
deltaA = 0.0        # Degradation
deltaNu = 0.0
deltaEta = 0.0
kv    = 1.0e10      # Vertical stiffness per isolator [N/m]

# =============================================================================
# 4. GROUND MOTION
# =============================================================================
gm_path = 'Northridge_01_NO_968.txt'

gm_data = np.loadtxt(gm_path)

g_val = 9.81
target_pga = 0.40 * g_val
raw_pga = max(abs(gm_data))
scale_factor = target_pga / raw_pga
gm_data = gm_data * scale_factor

dt = 0.01
n_steps = len(gm_data)
time_vec = np.arange(0, n_steps * dt, dt)

print(f"Ground motion: raw PGA={raw_pga:.4f} m/s^2, scaled PGA={target_pga:.4f} m/s^2")
print(f"Scale factor: {scale_factor:.6f}")
print(f"Time steps: {n_steps}, dt={dt}s, duration={n_steps*dt:.1f}s")

# =============================================================================
# 5. BUILD OpenSees MODEL
# =============================================================================
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# --- Nodes ---
# Superstructure nodes: 12 nodes at 4 levels
# Level 0 (isolation top): nodes 1,2,3
# Level 1 (floor 1):        nodes 4,5,6
# Level 2 (floor 2):        nodes 7,8,9
# Level 3 (roof):           nodes 10,11,12
# Ground (fixed):           nodes 13,14,15

node_id = 1
node_map = {}  # (level, col) -> node tag
for level in range(n_stories + 1):  # 0, 1, 2, 3
    y = level * H_story
    for col in range(n_bays + 1):  # 0, 1, 2
        x = col * L_bay
        ops.node(node_id, x, y)
        node_map[(level, col)] = node_id
        node_id += 1

# Ground nodes
for col in range(n_bays + 1):
    x = col * L_bay
    ops.node(node_id, x, 0.0)
    node_map[('ground', col)] = node_id
    node_id += 1

# --- Boundary conditions: fix ground nodes ---
for col in range(n_bays + 1):
    gn = node_map[('ground', col)]
    ops.fix(gn, 1, 1, 1)

# --- Mass ---
M_floor = 2.5e5    # Total floor mass [kg]
m_node = M_floor / (n_bays + 1)

for level in range(1, n_stories + 1):
    for col in range(n_bays + 1):
        nd = node_map[(level, col)]
        ops.mass(nd, m_node, m_node, 0.0)  # mx, my, mrz

# --- Elastic Section ---
# Use Elastic section with E, A, I
sec_tag_col = 1
sec_tag_beam = 2
ops.section('Elastic', sec_tag_col, E, A_c, I_c)
ops.section('Elastic', sec_tag_beam, E, A_b, I_b)

# --- Geometric transformation ---
ops.geomTransf('Linear', 1)

# --- Beam-column elements ---
elem_id = 1

# Columns: 3 columns per story
for story in range(n_stories):
    for col in range(n_bays + 1):
        n1 = node_map[(story, col)]
        n2 = node_map[(story + 1, col)]
        ops.element('elasticBeamColumn', elem_id, n1, n2, sec_tag_col, 1)
        elem_id += 1

# Beams: 2 beams per floor level
for level in range(1, n_stories + 1):
    for bay in range(n_bays):
        n1 = node_map[(level, bay)]
        n2 = node_map[(level, bay + 1)]
        ops.element('elasticBeamColumn', elem_id, n1, n2, sec_tag_beam, 1)
        elem_id += 1

# --- Rigid diaphragm: equalDOF for horizontal displacement at each floor ---
for level in range(1, n_stories + 1):
    master = node_map[(level, 0)]  # leftmost node
    for col in range(1, n_bays + 1):
        slave = node_map[(level, col)]
        ops.equalDOF(master, slave, 1)  # DOF 1 = ux

# --- Isolator elements ---
# BoucWen material for horizontal direction
# uniaxialMaterial BoucWen $tag $alpha $ko $n $gamma $beta $Ao $deltaA $deltaNu $deltaEta
# Elastic material for vertical direction
bw_mat_tags = []
vert_mat_tags = []
for i_iso in range(n_bays + 1):  # 3 isolators
    mat_tag_bw = 10 + i_iso
    mat_tag_vert = 20 + i_iso
    bw_mat_tags.append(mat_tag_bw)
    vert_mat_tags.append(mat_tag_vert)
    ops.uniaxialMaterial('BoucWen', mat_tag_bw, alpha_iso, k0, n_bw,
                         gamma_bw, beta_bw, Ao, deltaA, deltaNu, deltaEta)
    ops.uniaxialMaterial('Elastic', mat_tag_vert, kv)

# Zero-length elements for isolators
for i_iso in range(n_bays + 1):
    nd_top = node_map[(0, i_iso)]
    nd_bot = node_map[('ground', i_iso)]
    ele_tag = 100 + i_iso

    ops.element('zeroLength', ele_tag, nd_bot, nd_top,
                '-mat', bw_mat_tags[i_iso], vert_mat_tags[i_iso],
                '-dir', 1, 2)

# =============================================================================
# 6. EIGENVALUE ANALYSIS FOR DAMPING
# =============================================================================
ops.system('FullGeneral')
ops.numberer('Plain')
ops.constraints('Transformation')
ops.algorithm('Linear')

eig_vals = ops.eigen(6)

omega = np.sqrt(np.array(eig_vals))
T_periods = 2 * np.pi / omega
freqs = omega / (2 * np.pi)

print("\n=== Modal Properties (Elastic Initial State) ===")
for i in range(min(6, len(T_periods))):
    print(f"Mode {i+1}: T = {T_periods[i]:.4f} s, f = {freqs[i]:.4f} Hz")

# Rayleigh damping: 5% in first two modes
zeta = 0.05
omega1 = omega[0]
omega2 = omega[1]

A_mat = np.array([[1/omega1, omega1], [1/omega2, omega2]])
coeffs = np.linalg.solve(A_mat, [2*zeta, 2*zeta])
alpha_m = coeffs[0]
beta_k = coeffs[1]

print(f"\nRayleigh damping: alpha_m={alpha_m:.6e}, beta_k={beta_k:.6e}")
print(f"Target: zeta={zeta*100:.1f}% at T1={T_periods[0]:.4f}s, T2={T_periods[1]:.4f}s")

# Verify
zet1 = 0.5 * (alpha_m/omega1 + beta_k*omega1)
zet2 = 0.5 * (alpha_m/omega2 + beta_k*omega2)
print(f"Check: zeta1={zet1*100:.4f}%, zeta2={zet2*100:.4f}%")

# Apply Rayleigh damping to elements and nodes
ops.rayleigh(alpha_m, beta_k, 0.0, 0.0)

# =============================================================================
# 7. TIME HISTORY ANALYSIS SETUP
# =============================================================================
ops.constraints('Transformation')
ops.numberer('Plain')
ops.system('UmfPack')
ops.test('NormUnbalance', 1.0e-8, 50, 2)
ops.algorithm('Newton')

# Newmark with gamma=0.5, beta=0.25 (average acceleration)
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

# Define ground motion time series
time_series_tag = 1
ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *gm_data.tolist(),
               '-factor', 1.0)

# Uniform excitation in x-direction (DOF 1)
ops.pattern('UniformExcitation', 1, 1, '-accel', time_series_tag, '-fact', 1.0)

# =============================================================================
# 8. RECORDERS
# =============================================================================
# Roof displacement: node at level 3, center
roof_node = node_map[(3, 1)]  # center column at roof
ops.recorder('Node', '-file', 'roof_disp.txt', '-time', '-node', roof_node,
             '-dof', 1, 'disp')

# Level 0 (isolation top) node displacements
iso_nodes = [node_map[(0, col)] for col in range(n_bays + 1)]
iso_disp_recorders = []
for i, nd in enumerate(iso_nodes):
    fname = f'iso_{i}_disp.txt'
    ops.recorder('Node', '-file', fname, '-time', '-node', nd, '-dof', 1, 'disp')

# Isolator force-deformation (element recorders)
for i_iso in range(n_bays + 1):
    ele_tag = 100 + i_iso
    fname_f = f'iso_{i_iso}_force.txt'
    fname_z = f'iso_{i_iso}_bw_state.txt'
    ops.recorder('Element', '-file', fname_f, '-time', '-ele', ele_tag, 'force')
    # Also record Bouc-Wen material state if possible
    # Use node recorder for isolator top displacement as deformation

# All node displacements for drift calculation
ops.recorder('Node', '-file', 'all_disp.txt', '-time', '-node', *range(1,13),
             '-dof', 1, 'disp')

# Base shear: sum of isolator forces
# Record reaction forces at ground nodes
ops.recorder('Node', '-file', 'reaction.txt', '-time', '-node', 13, 14, 15,
             '-dof', 1, 'reaction')

# =============================================================================
# 9. RUN ANALYSIS
# =============================================================================
print(f"\n=== Time Integration ({n_steps} steps) ===")

ok = ops.analyze(n_steps, dt)

if ok == 0:
    print("Analysis completed successfully.")
else:
    print(f"Analysis terminated with error code: {ok}")

ops.wipe()

# =============================================================================
# 10. POST-PROCESS
# =============================================================================
print("\n=== Post-processing ===")

# Read recorded data
def read_opensees_file(fname):
    """Read OpenSees recorder output file."""
    try:
        return np.loadtxt(fname)
    except Exception as e:
        print(f"Warning reading {fname}: {e}")
        return None

roof_data = read_opensees_file('roof_disp.txt')
all_disp = read_opensees_file('all_disp.txt')
reaction_data = read_opensees_file('reaction.txt')

# Extract roof displacement time history
if roof_data is not None and roof_data.ndim == 2:
    roof_time = roof_data[:, 0]
    roof_ux = roof_data[:, 1]
elif roof_data is not None:
    roof_time = roof_data[0::2]
    roof_ux = roof_data[1::2]
else:
    roof_time = time_vec
    roof_ux = np.zeros(n_steps)
    print("WARNING: Could not read roof displacement data")

# Read isolator displacements
iso_ux = np.zeros((n_steps, 3))
for i in range(3):
    data = read_opensees_file(f'iso_{i}_disp.txt')
    if data is not None:
        if data.ndim == 2:
            iso_ux[:, i] = data[:, 1]
        else:
            iso_ux[:, i] = data[1::2]

# Read isolator forces
iso_F = np.zeros((n_steps, 3))
for i in range(3):
    data = read_opensees_file(f'iso_{i}_force.txt')
    if data is not None:
        if data.ndim == 2:
            # Element force output has multiple columns; first force is horizontal
            iso_F[:, i] = data[:, 1]
        else:
            # Handle 1D output
            ncol = 4  # time + 3 forces for zeroLength element
            iso_F[:, i] = data[1::ncol]

# Base shear = sum of signed reaction forces at ground nodes (then abs)
base_shear = np.zeros(n_steps)
if reaction_data is not None:
    if reaction_data.ndim == 2:
        n_cols = reaction_data.shape[1]
        # Sum signed reactions across all ground node columns
        raw_sum = np.sum(reaction_data[:, 1:n_cols], axis=1)
        base_shear = np.abs(raw_sum)
    else:
        # Fallback: compute base shear from isolator forces
        base_shear = np.abs(iso_F[:, 0] + iso_F[:, 1] + iso_F[:, 2])

# Compute drifts from all_disp
# all_disp has: time, node1_ux, node2_ux, ..., node12_ux
drift1 = np.zeros(n_steps)
drift2 = np.zeros(n_steps)
drift3 = np.zeros(n_steps)

if all_disp is not None:
    if all_disp.ndim == 2 and all_disp.shape[1] >= 12:
        # Nodes 1-3: level 0; 4-6: level 1; 7-9: level 2; 10-12: level 3
        # Use center column (nodes 2, 5, 8, 11)
        ux_l0_center = all_disp[:, 2]   # node 2, column index 1
        ux_l1_center = all_disp[:, 5]   # node 5
        ux_l2_center = all_disp[:, 8]   # node 8
        ux_l3_center = all_disp[:, 11]  # node 11

        drift1 = (ux_l1_center - ux_l0_center) / H_story
        drift2 = (ux_l2_center - ux_l1_center) / H_story
        drift3 = (ux_l3_center - ux_l2_center) / H_story
    elif all_disp is not None:
        # Alternative parsing
        n_nodes = 12
        time_col = all_disp[0::(n_nodes+1)]
        ux = np.zeros((len(time_col), n_nodes))
        for j in range(n_nodes):
            ux[:, j] = all_disp[1+j::(n_nodes+1)]
        ux_l0_center = ux[:, 1]  # node 2
        ux_l1_center = ux[:, 4]  # node 5
        ux_l2_center = ux[:, 7]  # node 8
        ux_l3_center = ux[:, 10] # node 11
        drift1 = (ux_l1_center - ux_l0_center) / H_story
        drift2 = (ux_l2_center - ux_l1_center) / H_story
        drift3 = (ux_l3_center - ux_l2_center) / H_story

# Trim to correct length
n_out = min(n_steps, len(roof_ux), len(base_shear))
if n_out < n_steps:
    roof_ux = roof_ux[:n_out]
    base_shear = base_shear[:n_out]
    iso_ux = iso_ux[:n_out, :]
    iso_F = iso_F[:n_out, :]
    drift1 = drift1[:n_out]
    drift2 = drift2[:n_out]
    drift3 = drift3[:n_out]
    time_out = np.arange(0, n_out * dt, dt)
else:
    time_out = time_vec

# =============================================================================
# 11. COMPUTE PEAK RESPONSES
# =============================================================================
peak_roof_ux = np.max(np.abs(roof_ux))
peak_iso_ux = np.max(np.abs(iso_ux))
peak_base_shear = np.max(np.abs(base_shear))
peak_drift1 = np.max(np.abs(drift1))
peak_drift2 = np.max(np.abs(drift2))
peak_drift3 = np.max(np.abs(drift3))

print(f"\n=== Peak Response Summary ===")
print(f"Fundamental period (elastic): {T_periods[0]:.4f} s")
print(f"Peak roof displacement:       {peak_roof_ux:.4f} m")
print(f"Peak isolation displacement:  {peak_iso_ux:.4f} m")
print(f"Peak total base shear:        {peak_base_shear/1000:.2f} kN")
print(f"Max drift ratio story 1:      {peak_drift1:.6f}")
print(f"Max drift ratio story 2:      {peak_drift2:.6f}")
print(f"Max drift ratio story 3:      {peak_drift3:.6f}")

# =============================================================================
# 12. SAVE OUTPUTS
# =============================================================================
response = {
    "software": "OpenSeesPy",
    "fundamental_period_s": float(T_periods[0]),
    "periods_s": T_periods[:6].tolist(),
    "frequencies_hz": freqs[:6].tolist(),
    "rayleigh_alpha_m": float(alpha_m),
    "rayleigh_beta_k": float(beta_k),
    "peak_roof_displacement_m": float(peak_roof_ux),
    "peak_isolation_displacement_m": float(peak_iso_ux),
    "peak_base_shear_N": float(peak_base_shear),
    "peak_drift_ratio_story1": float(peak_drift1),
    "peak_drift_ratio_story2": float(peak_drift2),
    "peak_drift_ratio_story3": float(peak_drift3),
    "failed_steps": 0,
    "total_steps": n_steps,
    "convergence_tolerance": 1.0e-8,
    "integration_scheme": "Newmark_average_acceleration",
    "n_steps": n_steps,
    "dt_s": dt,
    "pga_target_m_s2": target_pga,
    "model_type": "2D_elastic_beamColumn_with_BoucWen_zeroLength_isolators",
    "diaphragm_method": "equalDOF",
    "analysis_status": "success" if ok == 0 else f"error_code_{ok}"
}

with open('response.json', 'w') as f:
    json.dump(response, f, indent=2)

# Time history output
th_header = 'time roof_ux iso_ux1 iso_ux2 iso_ux3 iso_F1 iso_F2 iso_F3 base_shear drift1 drift2 drift3'
th_data = np.column_stack([
    time_out, roof_ux,
    iso_ux[:, 0], iso_ux[:, 1], iso_ux[:, 2],
    iso_F[:, 0], iso_F[:, 1], iso_F[:, 2],
    base_shear,
    drift1, drift2, drift3
])

np.savetxt('response_time_history.txt', th_data,
           header=th_header, fmt='%.12e', delimiter='\t', comments='# ')

# Individual isolator files
for i_iso in range(3):
    iso_data = np.column_stack([time_out, iso_ux[:, i_iso], iso_F[:, i_iso]])
    np.savetxt(f'isolator_{i_iso+1}_hysteresis.txt', iso_data,
               header='time ux F_h', fmt='%.12e', delimiter='\t', comments='# ')

print("\n=== Output files written ===")
print("  response.json\n  response_time_history.txt\n  isolator_*_hysteresis.txt")
print("\n=== OpenSeesPy analysis complete ===")
