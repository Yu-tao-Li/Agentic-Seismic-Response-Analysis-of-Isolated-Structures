import numpy as np

d1 = np.loadtxt('matlab/topStoDis.txt')
a1 = np.loadtxt('matlab/topStoAcc.txt')
d2 = np.loadtxt('openseespy/topStoDis.txt')
a2 = np.loadtxt('openseespy/topStoAcc.txt')

for name, mat, ops in [('Displacement', d1[:,1], d2[:,1]), ('Acceleration', a1[:,1], a2[:,1])]:
    pr = np.max(np.abs(mat))
    po = np.max(np.abs(ops))
    err = abs(pr-po)/max(pr,po)*100
    rng = np.max(ops) - np.min(ops)
    nrmse = np.sqrt(np.mean((mat-ops)**2))/rng
    print(f'{name}: MATLAB_peak={pr:.8f}, OPS_peak={po:.8f}, err={err:.6f}%, NRMSE={nrmse:.8f}')

# Also check isolation displacement
i1 = np.loadtxt('matlab/isoDis.txt')
i2 = np.loadtxt('openseespy/isoDis.txt')
pr = np.max(np.abs(i1[:,1]))
po = np.max(np.abs(i2[:,1]))
err = abs(pr-po)/max(pr,po)*100
rng = np.max(i2[:,1]) - np.min(i2[:,1])
nrmse = np.sqrt(np.mean((i1[:,1]-i2[:,1])**2))/rng
print(f'Iso Disp: MATLAB_peak={pr:.8f}, OPS_peak={po:.8f}, err={err:.6f}%, NRMSE={nrmse:.8f}')
