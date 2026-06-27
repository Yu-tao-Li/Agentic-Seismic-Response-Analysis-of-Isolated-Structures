from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import openseespy.opensees as ops

dt = 0.02
g = 9.81
target_pga = 0.40 * g
masses = np.array([1000.0, 1000.0, 1000.0], dtype=float)
story_k = np.array([500000.0, 500000.0, 500000.0], dtype=float)
c_extra = {}
ndof = 3
raw = np.loadtxt(Path(__file__).resolve().parents[2] / "input_data" / "GM1.txt", dtype=float).reshape(-1)
scale = target_pga / np.max(np.abs(raw))
ag = raw * scale
n = ag.size
t = np.arange(n) * dt
M = np.diag(masses)
K = np.zeros((ndof, ndof))
for spring, value in enumerate(story_k):
    if spring == 0:
        K[0, 0] += value
    else:
        left = spring - 1
        right = spring
        K[left, left] += value
        K[right, right] += value
        K[left, right] -= value
        K[right, left] -= value
eigvals = np.linalg.eigvals(np.linalg.solve(M, K))
omega = np.sqrt(np.sort(np.real(eigvals)))
xi = 0.05
alpha = 2.0 * xi * omega[0] * omega[1] / (omega[0] + omega[1])
beta_ray = 2.0 * xi / (omega[0] + omega[1])

ops.wipe()
ops.model("basic", "-ndm", 1, "-ndf", 1)
for node in range(1, ndof + 2):
    ops.node(node, 0.0)
ops.fix(1, 1)
for i, mass in enumerate(masses, start=2):
    ops.mass(i, float(mass))
for i, k in enumerate(story_k, start=1):
    ops.uniaxialMaterial("Elastic", i, float(k))
    ops.element("zeroLength", i, i, i + 1, "-mat", i, "-dir", 1, "-doRayleigh", 1)
for idx, c in c_extra.items():
    tag = 100 + idx
    ops.uniaxialMaterial("Viscous", tag, float(c), 1.0)
    ops.element("zeroLength", tag, 1, idx + 2, "-mat", tag, "-dir", 1)
ops.rayleigh(alpha, beta_ray, 0.0, 0.0)
ops.timeSeries("Path", 1, "-dt", dt, "-values", *ag.tolist())
ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
ops.constraints("Plain")
ops.numberer("Plain")
ops.system("UmfPack")
ops.test("NormUnbalance", 1.0e-10, 30, 0)
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
top_u = np.zeros(n)
top_a = np.zeros(n)
iso_u = np.zeros(n)
for step in range(1, n):
    ok = ops.analyze(1, dt)
    if ok != 0:
        raise RuntimeError(f"OpenSeesPy failed at step {step}")
    top_u[step] = ops.nodeDisp(ndof + 1, 1)
    top_a[step] = ops.nodeAccel(ndof + 1, 1)
    iso_u[step] = ops.nodeDisp(2, 1)
ops.wipe()
out = Path(__file__).resolve().parent
np.savetxt(out / "topStoDis.txt", np.column_stack([t, top_u]), fmt="%.10e")
np.savetxt(out / "topStoAcc.txt", np.column_stack([t, top_a]), fmt="%.10e")
np.savetxt(out / "isolation_displacement.txt", np.column_stack([t, iso_u]), fmt="%.10e")
with open(out / "response_output.json", "w", encoding="utf-8") as f:
    json.dump({
        "software": "OpenSeesPy",
        "time": t.tolist(),
        "top_displacement": top_u.tolist(),
        "top_acceleration": top_a.tolist(),
        "isolation_displacement": iso_u.tolist(),
        "peak_top_displacement_m": float(np.max(np.abs(top_u))),
        "peak_top_acceleration_m_s2": float(np.max(np.abs(top_a))),
        "peak_isolation_displacement_m": float(np.max(np.abs(iso_u))),
        "input_pga_normalized_m_s2": float(target_pga),
        "raw_pga_m_s2": float(np.max(np.abs(raw))),
        "scale_factor": float(scale),
        "dt": dt,
        "n_samples": int(n),
        "rayleigh_alpha": float(alpha),
        "rayleigh_beta": float(beta_ray),
    }, f, indent=2)
