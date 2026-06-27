"""Compare damped vs undamped responses for MATLAB model."""
import numpy as np, os

nStories = 3; m = 1000.0; k_story = 500e3; xi = 0.05; dt = 0.02; g = 9.81; targetPGA = 0.40 * g

M = m * np.eye(nStories)
K = np.zeros((nStories, nStories))
for s in range(nStories):
    if s == 0: K[s,s]=2*k_story; K[s,s+1]=-k_story
    elif s==nStories-1: K[s,s-1]=-k_story; K[s,s]=k_story
    else: K[s,s-1]=-k_story; K[s,s]=2*k_story; K[s,s+1]=-k_story

eigvals = np.linalg.eigvals(np.linalg.solve(M, K))
omega = np.sort(np.sqrt(np.abs(eigvals)))
alphaM = 2*xi*omega[0]*omega[1]/(omega[0]+omega[1])
betaK = 2*xi/(omega[0]+omega[1])
C = alphaM * M + betaK * K

base = os.path.dirname(os.path.abspath(__file__))
rawGM = np.loadtxt(os.path.join(base, 'input_data', 'GM1.txt'))
scale = targetPGA / np.max(np.abs(rawGM))
ugddot = rawGM * scale; np_pts = len(ugddot)

def run_newmark(damped=True):
    C_use = C if damped else np.zeros_like(C)
    gamma=0.5; beta=0.25
    a0=1/(beta*dt**2); a1=gamma/(beta*dt); a2=1/(beta*dt); a3=1/(2*beta)-1
    a4=gamma/beta-1; a5=dt/2*(gamma/beta-2); a6=dt*(1-gamma); a7=gamma*dt
    Keff = K + a0*M + a1*C_use
    infVec = np.ones((nStories, 1))
    u=np.zeros(nStories); ud=np.zeros(nStories); udd=np.zeros(nStories)
    topDis=np.zeros(np_pts); topAcc=np.zeros(np_pts)
    for i in range(np_pts):
        topDis[i]=u[-1]; topAcc[i]=udd[-1]
        if i==np_pts-1: break
        P_ext = -m * ugddot[i+1] * infVec.flatten()
        F_eff = P_ext + M@(a0*u+a2*ud+a3*udd) + C_use@(a1*u+a4*ud+a5*udd)
        u_new = np.linalg.solve(Keff, F_eff)
        udd_new = a0*(u_new-u) - a2*ud - a3*udd
        ud_new = ud + a6*udd + a7*udd_new
        u=u_new; ud=ud_new; udd=udd_new
    return topDis, topAcc

dis_d, acc_d = run_newmark(True)
dis_u, acc_u = run_newmark(False)

print(f"Damped peak disp:   {np.max(np.abs(dis_d)):.10e} at idx={np.argmax(np.abs(dis_d))}")
print(f"Undamped peak disp: {np.max(np.abs(dis_u)):.10e} at idx={np.argmax(np.abs(dis_u))}")
print(f"Ratio undamped/damped: {np.max(np.abs(dis_u))/np.max(np.abs(dis_d)):.6f}")
print()
print(f"Damped peak acc:   {np.max(np.abs(acc_d)):.10e}")
print(f"Undamped peak acc: {np.max(np.abs(acc_u)):.10e}")
print(f"Ratio undamped/damped: {np.max(np.abs(acc_u))/np.max(np.abs(acc_d)):.6f}")
