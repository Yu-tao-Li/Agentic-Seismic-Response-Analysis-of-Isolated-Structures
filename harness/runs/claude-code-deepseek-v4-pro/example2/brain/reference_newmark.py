"""Reference Newmark-beta solver for 3-story shear building."""
import numpy as np
import os, json

nStories = 3
m = 1000.0
k_story = 500e3
xi = 0.05
dt = 0.02
g = 9.81
targetPGA = 0.40 * g

M = m * np.eye(nStories)
K = np.zeros((nStories, nStories))
for s in range(nStories):
    if s == 0:
        K[s, s] = 2 * k_story
        K[s, s+1] = -k_story
    elif s == nStories - 1:
        K[s, s-1] = -k_story
        K[s, s] = k_story
    else:
        K[s, s-1] = -k_story
        K[s, s] = 2 * k_story
        K[s, s+1] = -k_story

eigvals, eigvecs = np.linalg.eig(np.linalg.solve(M, K))
omega = np.sort(np.sqrt(np.abs(eigvals)))
omega1, omega2 = omega[0], omega[1]
alphaM = 2 * xi * omega1 * omega2 / (omega1 + omega2)
betaK = 2 * xi / (omega1 + omega2)
C = alphaM * M + betaK * K

print(f"omega={omega}")
print(f"alphaM={alphaM:.6f}, betaK={betaK:.8f}")

base = os.path.dirname(os.path.abspath(__file__))
rawGM = np.loadtxt(os.path.join(base, 'input_data', 'GM1.txt'))
rawPGA = np.max(np.abs(rawGM))
scale = targetPGA / rawPGA
ugddot = rawGM * scale
np_pts = len(ugddot)

# Newmark-beta average acceleration
gamma = 0.5
beta = 0.25

a0 = 1.0 / (beta * dt**2)
a1 = gamma / (beta * dt)
a2 = 1.0 / (beta * dt)
a3 = 1.0/(2*beta) - 1.0
a4 = gamma/beta - 1.0
a5 = dt/2.0 * (gamma/beta - 2.0)
a6 = dt * (1.0 - gamma)
a7 = gamma * dt

Keff = K + a0*M + a1*C

u = np.zeros(nStories)
ud = np.zeros(nStories)
udd = np.zeros(nStories)

infVec = np.ones((nStories, 1))

time = np.zeros(np_pts)
topDis = np.zeros(np_pts)
topAcc = np.zeros(np_pts)

for i in range(np_pts):
    time[i] = (i) * dt
    topDis[i] = u[-1]
    topAcc[i] = udd[-1]
    if i == np_pts - 1:
        break

    # Force at end of step
    P_ext = -m * ugddot[i+1] * infVec.flatten()

    F_eff = P_ext + M @ (a0*u + a2*ud + a3*udd) + C @ (a1*u + a4*ud + a5*udd)
    u_new = np.linalg.solve(Keff, F_eff)
    udd_new = a0*(u_new - u) - a2*ud - a3*udd
    ud_new = ud + a6*udd + a7*udd_new

    u = u_new
    ud = ud_new
    udd = udd_new

out_dir = os.path.join(base, 'verification', 'reference')
os.makedirs(out_dir, exist_ok=True)

np.savetxt(os.path.join(out_dir, 'topStoDis.txt'), np.column_stack([time, topDis]), fmt='%.10e')
np.savetxt(os.path.join(out_dir, 'topStoAcc.txt'), np.column_stack([time, topAcc]), fmt='%.10e')

data = {
    "software": "Python_Reference_Newmark",
    "time_step": dt,
    "n_samples": np_pts,
    "top_displacement_peak_m": float(np.max(np.abs(topDis))),
    "top_acceleration_peak_m_s2": float(np.max(np.abs(topAcc))),
    "omega1": float(omega1), "omega2": float(omega2), "omega3": float(omega[2]),
    "alphaM": float(alphaM), "betaK": float(betaK)
}
with open(os.path.join(out_dir, 'response_output.json'), 'w') as f:
    json.dump(data, f, indent=2)

print(f"Peak displacement: {np.max(np.abs(topDis)):.6e} m")
print(f"Peak acceleration: {np.max(np.abs(topAcc)):.6e} m/s^2")

# Compare with MATLAB
mat_dis_raw = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoDis.txt'))
mat_dis = np.column_stack([mat_dis_raw[0::2], mat_dis_raw[1::2]])
mat_acc_raw = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoAcc.txt'))
mat_acc = np.column_stack([mat_acc_raw[0::2], mat_acc_raw[1::2]])

print(f"\nComparison with MATLAB:")
d_err = np.max(np.abs(topDis - mat_dis[:,1]))
a_err = np.max(np.abs(topAcc - mat_acc[:,1]))
print(f"  Max displacement diff: {d_err:.6e} m")
print(f"  Max acceleration diff: {a_err:.6e} m/s^2")

# Compare with ABAQUS
aba_dis = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoDis.txt'))
aba_acc = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoAcc.txt'))

n_comp = min(len(topDis), len(aba_dis))
d_err_aba = np.max(np.abs(topDis[:n_comp] - aba_dis[:n_comp,1]))
a_err_aba = np.max(np.abs(topAcc[:n_comp] - aba_acc[:n_comp,1]))
print(f"\nComparison with ABAQUS:")
print(f"  Max displacement diff: {d_err_aba:.6e} m")
print(f"  Max acceleration diff: {a_err_aba:.6e} m/s^2")
print(f"  ABAQUS/Reference ratio (peak dis): {np.max(np.abs(aba_dis[:,1]))/np.max(np.abs(topDis)):.6f}")
print(f"  ABAQUS/Reference ratio (peak acc): {np.max(np.abs(aba_acc[:,1]))/np.max(np.abs(topAcc)):.6f}")
