"""Check whether the Abaqus ODB satisfies the MATLAB MDOF equation."""
from odbAccess import openOdb
import numpy as np

m = np.array([250, 270, 270, 180]) * 1e3
k = np.array([50, 245, 195, 98]) * 1e6
dt = 0.01
g = 9.81
PGA = 0.40 * g

acc_raw = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
ag = acc_raw / np.max(np.abs(acc_raw)) * PGA

M = np.diag(m)
K0 = np.zeros((4, 4))
K0[0, 0] = k[0] + k[1]
K0[0, 1] = -k[1]
K0[1, 0] = -k[1]
K0[1, 1] = k[1] + k[2]
K0[1, 2] = -k[2]
K0[2, 1] = -k[2]
K0[2, 2] = k[2] + k[3]
K0[2, 3] = -k[3]
K0[3, 2] = -k[3]
K0[3, 3] = k[3]

eigvals, _ = np.linalg.eig(np.linalg.solve(M, K0))
omega = np.sqrt(np.sort(np.real(eigvals)))
alpha = 2 * 0.05 * omega[0] * omega[1] / (omega[0] + omega[1])
beta = 2 * 0.05 / (omega[0] + omega[1])
C = alpha * M + beta * K0

odb = openOdb('example5_model.odb')
step = odb.steps['DYNAMIC']
inst = odb.rootAssembly.instances['PART-1-1']
node_regions = [inst.nodes[i] for i in range(1, 5)]
spring_regions = [inst.elementSets[f'SPRING{i}'] for i in range(1, 5)]

n = len(step.frames)
u = np.zeros((n, 4))
v = np.zeros((n, 4))
a = np.zeros((n, 4))
shear = np.zeros((n, 4))
time = np.zeros(n)

for j, frame in enumerate(step.frames):
    time[j] = frame.frameValue
    for i, node in enumerate(node_regions):
        u[j, i] = frame.fieldOutputs['U'].getSubset(region=node).values[0].data[0]
        v[j, i] = frame.fieldOutputs['V'].getSubset(region=node).values[0].data[0]
        a[j, i] = frame.fieldOutputs['A'].getSubset(region=node).values[0].data[0]
    for i, region in enumerate(spring_regions):
        s_data = frame.fieldOutputs['S'].getSubset(region=region).values[0].data
        shear[j, i] = s_data[0] if hasattr(s_data, '__len__') else s_data

odb.close()

Rs = np.zeros_like(shear)
Rs[:, 3] = shear[:, 3]
Rs[:, 0] = shear[:, 0] - shear[:, 1]
Rs[:, 1] = shear[:, 1] - shear[:, 2]
Rs[:, 2] = shear[:, 2] - shear[:, 3]

ag = ag[:n]
res = np.zeros_like(u)
for j in range(n):
    res[j, :] = M.dot(a[j, :]) + C.dot(v[j, :]) + Rs[j, :] + M.dot(np.ones(4) * ag[j])

ground_force = ag[:, None] * m[None, :]
scale = np.max(np.abs(ground_force))
print('alpha', alpha, 'beta', beta)
print('max_abs_residual', np.max(np.abs(res)))
print('max_abs_residual_over_peak_ground_force', np.max(np.abs(res)) / scale)
for idx in [0, 368, 580, 900, 1173, 1400, 1899]:
    print('t', time[idx], 'res', res[idx, :], 'u', u[idx, :])
