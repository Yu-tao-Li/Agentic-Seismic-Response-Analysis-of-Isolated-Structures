"""Comprehensive three-software comparison for Example 2."""
import numpy as np, json, os

base = os.path.dirname(os.path.abspath(__file__))

# Load MATLAB
mat_raw = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoDis.txt'))
mat_dis = np.column_stack([mat_raw[0::2], mat_raw[1::2]])
mat_raw_a = np.loadtxt(os.path.join(base, 'verification', 'matlab', 'topStoAcc.txt'))
mat_acc = np.column_stack([mat_raw_a[0::2], mat_raw_a[1::2]])

# Load OpenSeesPy
ops_dis = np.loadtxt(os.path.join(base, 'verification', 'openseespy', 'topStoDis.txt'))
ops_acc = np.loadtxt(os.path.join(base, 'verification', 'openseespy', 'topStoAcc.txt'))

# Load ABAQUS
aba_dis = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoDis.txt'))
aba_acc = np.loadtxt(os.path.join(base, 'verification', 'abaqus', 'topStoAcc.txt'))

def compute_metrics(ref_dis, ref_acc, test_dis, test_acc, label):
    n = min(len(ref_dis), len(test_dis))
    rd = ref_dis[:n, 1]; td = test_dis[:n, 1]
    ra = ref_acc[:n, 1]; ta = test_acc[:n, 1]

    peak_d_err = abs(np.max(np.abs(td)) - np.max(np.abs(rd))) / np.max(np.abs(rd))
    peak_a_err = abs(np.max(np.abs(ta)) - np.max(np.abs(ra))) / np.max(np.abs(ra))

    rmse_d = np.sqrt(np.mean((td - rd)**2))
    nrmse_d = rmse_d / (np.max(np.abs(rd)) - np.min(rd)) if np.max(np.abs(rd)) > np.min(rd) else 0

    rmse_a = np.sqrt(np.mean((ta - ra)**2))
    nrmse_a = rmse_a / (np.max(np.abs(ra)) - np.min(ra)) if np.max(np.abs(ra)) > np.min(ra) else 0

    l2_d = np.sqrt(np.sum((td - rd)**2) / np.sum(rd**2))
    l2_a = np.sqrt(np.sum((ta - ra)**2) / np.sum(ra**2))

    print(f"\n{label}:")
    print(f"  Peak displacement: ref={np.max(np.abs(rd)):.10e}  test={np.max(np.abs(td)):.10e}  rel_err={peak_d_err:.6e}")
    print(f"  Peak acceleration: ref={np.max(np.abs(ra)):.10e}  test={np.max(np.abs(ta)):.10e}  rel_err={peak_a_err:.6e}")
    print(f"  Displacement NRMSE: {nrmse_d:.6e}")
    print(f"  Acceleration NRMSE: {nrmse_a:.6e}")
    print(f"  Displacement L2 error: {l2_d:.6e}")
    print(f"  Acceleration L2 error: {l2_a:.6e}")
    print(f"  Aligned samples: {n}")

    return {
        "displacement_peak_relative_error": float(peak_d_err),
        "acceleration_peak_relative_error": float(peak_a_err),
        "displacement_nrmse": float(nrmse_d),
        "acceleration_nrmse": float(nrmse_a),
        "displacement_relative_l2_error": float(l2_d),
        "acceleration_relative_l2_error": float(l2_a),
        "aligned_samples": n
    }

print("="*60)
print("THREE-SOFTWARE CROSS-VALIDATION")
print("="*60)

print(f"\nMATLAB:  {mat_dis.shape} disp, {mat_acc.shape} acc")
print(f"OpenSeesPy: {ops_dis.shape} disp, {ops_acc.shape} acc")
print(f"ABAQUS: {aba_dis.shape} disp, {aba_acc.shape} acc")

# Use MATLAB as reference (reference Python solver confirms it's correct)
mat_vs_ops = compute_metrics(mat_dis, mat_acc, ops_dis, ops_acc, "MATLAB vs OpenSeesPy")
mat_vs_aba = compute_metrics(mat_dis, mat_acc, aba_dis, aba_acc, "MATLAB vs ABAQUS")
ops_vs_aba = compute_metrics(ops_dis, ops_acc, aba_dis, aba_acc, "OpenSeesPy vs ABAQUS")

# Thresholds (from report)
d_peak_thresh = 0.01
a_peak_thresh = 0.02
d_nrmse_thresh = 0.02
a_nrmse_thresh = 0.05

def check_pass(metrics):
    return (metrics["displacement_peak_relative_error"] <= d_peak_thresh and
            metrics["acceleration_peak_relative_error"] <= a_peak_thresh and
            metrics["displacement_nrmse"] <= d_nrmse_thresh and
            metrics["acceleration_nrmse"] <= a_nrmse_thresh)

m_o_pass = check_pass(mat_vs_ops)
m_a_pass = check_pass(mat_vs_aba)
o_a_pass = check_pass(ops_vs_aba)

print(f"\nThresholds: disp_peak<{d_peak_thresh}, acc_peak<{a_peak_thresh}, disp_nrmse<{d_nrmse_thresh}, acc_nrmse<{a_nrmse_thresh}")
print(f"MATLAB vs OpenSeesPy: {'PASS' if m_o_pass else 'FAIL'}")
print(f"MATLAB vs ABAQUS:     {'PASS' if m_a_pass else 'FAIL'}")
print(f"OpenSeesPy vs ABAQUS: {'PASS' if o_a_pass else 'FAIL'}")
print(f"\nOVERALL THREE-SOFTWARE AGREEMENT: {'PASS' if (m_o_pass and m_a_pass and o_a_pass) else 'FAIL'}")

# Save results
result = {
    "comparisons": {
        "matlab_vs_openseespy": mat_vs_ops,
        "matlab_vs_abaqus": mat_vs_aba,
        "openseespy_vs_abaqus": ops_vs_aba
    },
    "thresholds": {
        "displacement_peak_relative_error": d_peak_thresh,
        "acceleration_peak_relative_error": a_peak_thresh,
        "displacement_nrmse": d_nrmse_thresh,
        "acceleration_nrmse": a_nrmse_thresh
    },
    "overall_pass": m_o_pass and m_a_pass and o_a_pass
}

with open(os.path.join(base, 'cross_validation_report.json'), 'w') as f:
    json.dump(result, f, indent=2)
print("\nResults written to cross_validation_report.json")
