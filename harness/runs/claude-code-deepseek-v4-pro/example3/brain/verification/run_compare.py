import numpy as np, json, os
base = r'harness\runs\claude-code-deepseek-v4-pro\example3\brain\verification'

def load(dir_name):
    d = np.loadtxt(os.path.join(base, dir_name, 'topStoDis.txt'))
    a = np.loadtxt(os.path.join(base, dir_name, 'topStoAcc.txt'))
    i = np.loadtxt(os.path.join(base, dir_name, 'isoDis.txt'))
    return d[:,1], a[:,1], i[:,1]

d1,a1,i1 = load('matlab')
d2,a2,i2 = load('openseespy')

for name, m, o in [('Displacement', d1, d2), ('Acceleration', a1, a2), ('Iso Disp', i1, i2)]:
    pm = np.max(np.abs(m)); po = np.max(np.abs(o))
    err = abs(pm-po)/max(pm,po)*100
    rng = np.max(o) - np.min(o)
    nrmse = np.sqrt(np.mean((m-o)**2))/rng
    print(f'{name}: peak_err={err:.6f}%, NRMSE={nrmse:.8f}, MATLAB={pm:.6f}, OPS={po:.6f}')

# Write initial cross-validation (pre-ABAQUS)
report = {
    'matlab_vs_openseespy': {
        'displacement_peak_error_pct': float(abs(np.max(np.abs(d1))-np.max(np.abs(d2)))/max(np.max(np.abs(d1)),np.max(np.abs(d2)))*100),
        'acceleration_peak_error_pct': float(abs(np.max(np.abs(a1))-np.max(np.abs(a2)))/max(np.max(np.abs(a1)),np.max(np.abs(a2)))*100),
        'isolation_peak_error_pct': float(abs(np.max(np.abs(i1))-np.max(np.abs(i2)))/max(np.max(np.abs(i1)),np.max(np.abs(i2)))*100),
        'displacement_nrmse': float(np.sqrt(np.mean((d1-d2)**2))/(np.max(d2)-np.min(d2))),
        'acceleration_nrmse': float(np.sqrt(np.mean((a1-a2)**2))/(np.max(a2)-np.min(a2))),
    }
}
print('\nPreliminary report:', json.dumps(report, indent=2))
