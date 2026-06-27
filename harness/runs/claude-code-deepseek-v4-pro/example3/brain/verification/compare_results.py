"""Compare results across MATLAB, OpenSeesPy, and optionally ABAQUS."""
import numpy as np
import json
import os
import sys

base = os.path.dirname(os.path.abspath(__file__))

def load_results(path, label):
    """Load topStoDis.txt and topStoAcc.txt from given directory."""
    dis_file = os.path.join(path, 'topStoDis.txt')
    acc_file = os.path.join(path, 'topStoAcc.txt')
    iso_file = os.path.join(path, 'isoDis.txt')

    if not os.path.exists(dis_file) or not os.path.exists(acc_file):
        return None, f"{label}: missing files"

    dis = np.loadtxt(dis_file)
    acc = np.loadtxt(acc_file)
    iso = None
    if os.path.exists(iso_file):
        iso = np.loadtxt(iso_file)

    t = dis[:, 0]
    d = dis[:, 1]
    a = acc[:, 1]
    iso_d = iso[:, 1] if iso is not None else None

    # Verify time consistency
    dt = t[1] - t[0] if len(t) > 1 else 0.02

    result = {
        'label': label,
        'time': t,
        'top_displacement': d,
        'top_acceleration': a,
        'isolation_displacement': iso_d,
        'n_samples': len(t),
        'dt': dt,
        'peak_dis': float(np.max(np.abs(d))),
        'peak_acc': float(np.max(np.abs(a))),
        'peak_iso': float(np.max(np.abs(iso_d))) if iso_d is not None else None,
    }
    return result, None

def compute_metrics(ref, other, name):
    """Compute peak error and NRMSE between ref and other."""
    if len(ref) != len(other):
        other = other[:len(ref)]

    peak_ref = np.max(np.abs(ref))
    peak_other = np.max(np.abs(other))
    peak_err_pct = abs(peak_ref - peak_other) / max(peak_ref, peak_other) * 100

    rng = np.max(ref) - np.min(ref)
    if rng < 1e-10:
        rng = peak_ref * 2
    nrmse = np.sqrt(np.mean((ref - other) ** 2)) / rng

    return {
        'peak_ref': float(peak_ref),
        'peak_other': float(peak_other),
        'peak_error_pct': float(peak_err_pct),
        'nrmse': float(nrmse),
    }

def main():
    paths = {
        'matlab': os.path.join(base, 'matlab'),
        'openseespy': os.path.join(base, 'openseespy'),
        'abaqus': os.path.join(base, 'abaqus'),
    }

    results = {}
    errors = {}

    for name, path in paths.items():
        r, err = load_results(path, name)
        if r:
            results[name] = r
            print(f"[{name}] N={r['n_samples']}, dt={r['dt']:.3f}s, "
                  f"peak_dis={r['peak_dis']:.6f}, "
                  f"peak_acc={r['peak_acc']:.4f}, "
                  f"peak_iso={r['peak_iso']}")
        else:
            errors[name] = err
            print(f"[{name}] ERROR: {err}")

    # Cross-validation
    report = {
        'input_consistency': {},
        'comparisons': {},
        'errors': errors,
        'consistency_thresholds': {},
        'overall_pass': None,
    }

    # Check input consistency
    n_samples = set()
    dts = set()
    for name, r in results.items():
        n_samples.add(r['n_samples'])
        dts.add(round(r['dt'], 6))
    report['input_consistency'] = {
        'n_samples_match': len(n_samples) <= 1,
        'n_samples_values': list(n_samples),
        'dt_match': len(dts) <= 1,
        'dt_values': list(dts),
        'pga_normalized': 0.40 * 9.81,
    }
    print(f"\nInput consistency: n_samples={n_samples}, dt={dts}")

    # Compare pairs
    pairs = [('matlab', 'openseespy')]
    if 'abaqus' in results:
        pairs.append(('matlab', 'abaqus'))
        pairs.append(('openseespy', 'abaqus'))

    # Adopted consistency thresholds
    thresholds = {
        'displacement_peak_error_pct': 2.0,    # 2% peak displacement error
        'acceleration_peak_error_pct': 2.0,    # 2% peak acceleration error
        'iso_displacement_peak_error_pct': 2.0, # 2% peak isolation displacement error
        'displacement_nrmse': 0.02,             # NRMSE < 0.02
        'acceleration_nrmse': 0.03,             # NRMSE < 0.03 (acceleration is more sensitive to numerical differences)
    }
    report['consistency_thresholds'] = thresholds

    for (a, b) in pairs:
        ra, rb = results[a], results[b]
        key = f"{a}_vs_{b}"

        dis_metrics = compute_metrics(ra['top_displacement'], rb['top_displacement'], 'displacement')
        acc_metrics = compute_metrics(ra['top_acceleration'], rb['top_acceleration'], 'acceleration')
        iso_metrics = compute_metrics(
            ra['isolation_displacement'], rb['isolation_displacement'], 'iso_displacement'
        )

        comp = {
            'displacement': dis_metrics,
            'acceleration': acc_metrics,
            'isolation_displacement': iso_metrics,
        }
        report['comparisons'][key] = comp

        print(f"\n{a} vs {b}:")
        print(f"  Displacement: peak_err={dis_metrics['peak_error_pct']:.4f}%, nrmse={dis_metrics['nrmse']:.6f}")
        print(f"  Acceleration: peak_err={acc_metrics['peak_error_pct']:.4f}%, nrmse={acc_metrics['nrmse']:.6f}")
        print(f"  Iso Disp:     peak_err={iso_metrics['peak_error_pct']:.4f}%, nrmse={iso_metrics['nrmse']:.6f}")

    # Overall pass/fail
    all_pass = True
    for key, comp in report['comparisons'].items():
        for metric, thresh_key in [
            ('displacement', 'displacement_peak_error_pct'),
            ('acceleration', 'acceleration_peak_error_pct'),
            ('isolation_displacement', 'iso_displacement_peak_error_pct'),
        ]:
            err_val = comp[metric]['peak_error_pct']
            threshold = thresholds[thresh_key]
            if err_val > threshold:
                all_pass = False
                print(f"FAIL: {key} {metric} peak_err={err_val:.4f}% > {threshold}%")
        for metric, thresh_key in [
            ('displacement', 'displacement_nrmse'),
            ('acceleration', 'acceleration_nrmse'),
        ]:
            nrmse_val = comp[metric]['nrmse']
            threshold = thresholds[thresh_key]
            if nrmse_val > threshold:
                all_pass = False
                print(f"FAIL: {key} {metric} NRMSE={nrmse_val:.6f} > {threshold}")

    report['overall_pass'] = all_pass
    if all_pass:
        print("\n=== ALL COMPARISONS PASS ===")
    else:
        print("\n=== SOME COMPARISONS FAILED ===")

    out_path = os.path.join(base, '..', 'cross_validation_report.json')
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"Report written to {out_path}")

    return report

if __name__ == '__main__':
    main()
