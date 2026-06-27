"""
Extract time history results from ABAQUS v9 ODB for comparison.
Run with: abaqus python extract_results_v9.py
"""
import os
import json

odb_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(odb_dir)

try:
    from odbAccess import openOdb
    from abaqusConstants import *

    odb = openOdb('frame_bw_dynamic_v9.odb', readOnly=True)

    step_name = odb.steps.keys()[-1]
    step = odb.steps[step_name]
    print(f"Step: {step_name}")
    print(f"Number of frames: {len(step.frames)}")

    # Extract time history from field output (all frames)
    print("Extracting full time histories...")
    n_frames = len(step.frames)
    times = []
    roof_ux = []
    iso1_ux = []
    iso2_ux = []
    iso3_ux = []
    floor1_ux = []
    floor2_ux = []
    floor3_ux = []

    for i, frame in enumerate(step.frames):
        t = frame.frameValue
        times.append(t)

        u_field = frame.fieldOutputs['U']
        u_dict = {}
        for val in u_field.values:
            u_dict[int(val.nodeLabel)] = float(val.data[0])

        roof_ux.append(u_dict.get(11, 0.0))
        iso1_ux.append(u_dict.get(1, 0.0))
        iso2_ux.append(u_dict.get(2, 0.0))
        iso3_ux.append(u_dict.get(3, 0.0))
        floor1_ux.append(u_dict.get(5, 0.0))
        floor2_ux.append(u_dict.get(8, 0.0))
        floor3_ux.append(u_dict.get(11, 0.0))

        if (i+1) % 5000 == 0:
            print(f"  Frame {i+1}/{n_frames}: t={t:.3f}s")

    # Save time history
    output = {
        'software': 'ABAQUS v9',
        'max_time_s': max(times),
        'n_frames': n_frames,
        'times': times,
        'roof_ux_m': roof_ux,
        'iso1_ux_m': iso1_ux,
        'iso2_ux_m': iso2_ux,
        'iso3_ux_m': iso3_ux,
        'floor1_ux_m': floor1_ux,
        'floor2_ux_m': floor2_ux,
        'floor3_ux_m': floor3_ux,
        'peak_roof_displacement_m': max(abs(r) for r in roof_ux),
        'peak_isolation_displacement_m': max(
            max(abs(v) for v in iso1_ux),
            max(abs(v) for v in iso2_ux),
            max(abs(v) for v in iso3_ux)
        ),
    }

    with open('abaqus_time_history_v9.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nPeak roof displacement (0-{max(times):.1f}s): {output['peak_roof_displacement_m']:.6f} m")
    print(f"Peak isolation displacement (0-{max(times):.1f}s): {output['peak_isolation_displacement_m']:.6f} m")
    print(f"Results saved to abaqus_time_history_v9.json")

    # Also save final displacements
    last_frame = step.frames[-1]
    u_field = last_frame.fieldOutputs['U']
    node_data = {}
    for val in u_field.values:
        node_data[int(val.nodeLabel)] = {
            'U1': float(val.data[0]),
            'U2': float(val.data[1]) if len(val.data) > 1 else 0.0,
        }

    with open('abaqus_final_displacements_v9.json', 'w') as f:
        json.dump(node_data, f, indent=2)

    for node_id in [1, 2, 3, 4, 5, 6, 10, 11, 12]:
        if node_id in node_data:
            print(f"  Node {node_id}: U1={node_data[node_id]['U1']:.6f}, U2={node_data[node_id]['U2']:.6e}")

    odb.close()

except ImportError as e:
    print(f"ERROR: {e}")
    import sys
    sys.exit(1)
