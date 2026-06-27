"""
Build cross_validation_report.json by comparing MATLAB, OpenSeesPy, and ABAQUS results.
Mode shapes are compared by absolute normalized amplitudes since signs are arbitrary.
"""
import json
import math

def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

matlab = load_json('verification/matlab/modal_results.json')
opensees = load_json('verification/openseespy/modal_results.json')
abaqus = load_json('verification/abaqus/modal_results.json')

def abs_mode_shape(ms):
    """Return mode shape with absolute values."""
    return [[abs(v) for v in mode] for mode in ms]

def freq_relative_error(f1, f2):
    """Relative error between two frequencies."""
    return abs(f1 - f2) / f1 if f1 != 0 else 0.0

def mac(phi1, phi2):
    """Modal Assurance Criterion between two mode shapes."""
    dot = sum(a * b for a, b in zip(phi1, phi2))
    norm1 = math.sqrt(sum(a * a for a in phi1))
    norm2 = math.sqrt(sum(b * b for b in phi2))
    if norm1 == 0 or norm2 == 0:
        return 0.0
    return (dot / (norm1 * norm2)) ** 2

def find_best_match(source_modes, target_modes):
    """For each source mode, find the best matching target mode by absolute MAC."""
    matches = {}
    used = set()
    for i, src in enumerate(source_modes):
        best_mac = -1
        best_j = -1
        for j, tgt in enumerate(target_modes):
            if j in used:
                continue
            m = mac(src, tgt)
            if m > best_mac:
                best_mac = m
                best_j = j
        matches[i] = {'matched_mode': best_j, 'mac': round(best_mac, 10)}
        used.add(best_j)
    return matches

# Compare frequencies
results = {}
software_pairs = [
    ('matlab', 'openseespy', matlab, opensees),
    ('matlab', 'abaqus', matlab, abaqus),
    ('openseespy', 'abaqus', opensees, abaqus),
]

freq_comparisons = {}
for name1, name2, data1, data2 in software_pairs:
    pair_key = f"{name1}_vs_{name2}"
    freq_errs = []
    for i in range(3):
        f1 = data1['circular_frequencies_rad_per_s'][i]
        f2 = data2['circular_frequencies_rad_per_s'][i]
        freq_errs.append(round(freq_relative_error(f1, f2), 10))
    freq_comparisons[pair_key] = {
        'circular_freq_relative_errors': freq_errs,
        'max_relative_error': max(freq_errs)
    }

# Compare mode shapes (absolute values)
ms_matlab = abs_mode_shape(matlab['normalized_mode_shapes'])
ms_opensees = abs_mode_shape(opensees['normalized_mode_shapes'])
ms_abaqus = abs_mode_shape(abaqus['normalized_mode_shapes'])

mode_shape_comparisons = {}
ms_pairs = [
    ('matlab', 'openseespy', ms_matlab, ms_opensees),
    ('matlab', 'abaqus', ms_matlab, ms_abaqus),
    ('openseespy', 'abaqus', ms_opensees, ms_abaqus),
]

for name1, name2, ms1, ms2 in ms_pairs:
    pair_key = f"{name1}_vs_{name2}"
    matches = find_best_match(ms1, ms2)
    mac_values = [matches[i]['mac'] for i in range(3)]
    mode_shape_comparisons[pair_key] = {
        'mode_matches': matches,
        'min_mac': min(mac_values),
        'all_mac_above_099': all(m >= 0.99 for m in mac_values)
    }

# Thresholds and pass/fail
freq_threshold = 0.001  # 0.1% relative error
mac_threshold = 0.99

all_freq_pass = all(
    v['max_relative_error'] < freq_threshold
    for v in freq_comparisons.values()
)
all_mac_pass = all(
    v['min_mac'] >= mac_threshold
    for v in mode_shape_comparisons.values()
)
overall_pass = all_freq_pass and all_mac_pass

report = {
    "analysis_description": "3-story shear building modal analysis cross-validation",
    "model_parameters": {
        "story_mass_kg": 1000.0,
        "story_stiffness_N_per_m": 500000.0,
        "num_stories": 3
    },
    "consistency_thresholds": {
        "max_circular_frequency_relative_error": freq_threshold,
        "min_mode_shape_MAC": mac_threshold,
        "rationale": "0.1% freq tolerance covers floating-point and solver differences; MAC >= 0.99 ensures mode shapes are essentially identical up to sign"
    },
    "frequency_comparisons": freq_comparisons,
    "mode_shape_comparisons": mode_shape_comparisons,
    "pass_fail": {
        "frequency_check": "PASS" if all_freq_pass else "FAIL",
        "mode_shape_check": "PASS" if all_mac_pass else "FAIL",
        "overall": "PASS" if overall_pass else "FAIL"
    }
}

with open('cross_validation_report.json', 'w') as f:
    json.dump(report, f, indent=2)

print("cross_validation_report.json written")
print(f"Overall: {'PASS' if overall_pass else 'FAIL'}")
