import numpy as np
import os

base = os.path.dirname(os.path.abspath(__file__))

# MATLAB writes alternating [time, disp, time, acc, ...] on separate lines
mat_dis_raw = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoDis.txt'))
mat_dis = np.column_stack([mat_dis_raw[0::2], mat_dis_raw[1::2]])

mat_acc_raw = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoAcc.txt'))
mat_acc = np.column_stack([mat_acc_raw[0::2], mat_acc_raw[1::2]])

aba_dis = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoDis.txt'))
aba_acc = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoAcc.txt'))

print(f'MATLAB dis shape: {mat_dis.shape}, ABAQUS dis shape: {aba_dis.shape}')
print(f'MATLAB acc shape: {mat_acc.shape}, ABAQUS acc shape: {aba_acc.shape}')
print()

# Interpolate ABAQUS to MATLAB time points for comparison
from scipy.interpolate import interp1d
# Actually, let's just compare at the intersection
n_cmp = min(len(mat_dis), len(aba_dis))
print(f'First 20 displacement points (n_cmp={n_cmp}):')
for i in range(min(20, n_cmp)):
    r = aba_dis[i,1]/mat_dis[i,1] if abs(mat_dis[i,1])>1e-12 else float('inf')
    print(f'  t={mat_dis[i,0]:.4f}  MATLAB={mat_dis[i,1]:.10e}  ABAQUS={aba_dis[i,1]:.10e}  ratio={r:.6f}')

print()
mdp = np.argmax(np.abs(mat_dis[:,1]))
adp = np.argmax(np.abs(aba_dis[:,1]))
print(f'Peak displacement:')
print(f'  MATLAB: {np.max(np.abs(mat_dis[:,1])):.10e} at t={mat_dis[mdp,0]:.4f} (idx={mdp})')
print(f'  ABAQUS: {np.max(np.abs(aba_dis[:,1])):.10e} at t={aba_dis[adp,0]:.4f} (idx={adp})')
print(f'  Ratio: {np.max(np.abs(aba_dis[:,1]))/np.max(np.abs(mat_dis[:,1])):.6f}')

map_ = np.argmax(np.abs(mat_acc[:,1]))
aap = np.argmax(np.abs(aba_acc[:,1]))
print(f'Peak acceleration:')
print(f'  MATLAB: {np.max(np.abs(mat_acc[:,1])):.10e} at t={mat_acc[map_,0]:.4f} (idx={map_})')
print(f'  ABAQUS: {np.max(np.abs(aba_acc[:,1])):.10e} at t={aba_acc[aap,0]:.4f} (idx={aap})')
print(f'  Ratio: {np.max(np.abs(aba_acc[:,1]))/np.max(np.abs(mat_acc[:,1])):.6f}')

# Check time alignment
print()
print('Time alignment check:')
for i in [0, 1, 100, 1000, 2685]:
    if i < len(mat_dis) and i < len(aba_dis):
        print(f'  i={i}: mat_t={mat_dis[i,0]:.6f}  aba_t={aba_dis[i,0]:.6f}  dt_diff={abs(mat_dis[i,0]-aba_dis[i,0]):.6e}')
