import numpy as np, os
base = os.path.dirname(os.path.abspath(__file__))

mat_raw = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoDis.txt'))
mat_dis = np.column_stack([mat_raw[0::2], mat_raw[1::2]])
aba_dis = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoDis.txt'))

print("First 20 displacement points:")
for i in range(20):
    r = aba_dis[i,1]/mat_dis[i,1] if abs(mat_dis[i,1])>1e-12 else -1
    print(f't={mat_dis[i,0]:.4f} MATLAB={mat_dis[i,1]:.10e} ABAQUS={aba_dis[i,1]:.10e} ratio={r:.6f}')

print(f"\nMATLAB peak: {np.max(np.abs(mat_dis[:,1])):.10e} at idx={np.argmax(np.abs(mat_dis[:,1]))}")
print(f"ABAQUS peak: {np.max(np.abs(aba_dis[:,1])):.10e} at idx={np.argmax(np.abs(aba_dis[:,1]))}")
