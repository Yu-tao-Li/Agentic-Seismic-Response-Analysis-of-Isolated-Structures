"""
Extract time history results from ABAQUS ODB for comparison.
Run with: abaqus python extract_results.py
"""
import sys
import os
import json

# Change to the ODB directory
odb_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(odb_dir)

try:
    from odbAccess import openOdb
    from abaqusConstants import *

    odb = openOdb('frame_bw_dynamic.odb', readOnly=True)

    # Get the last step
    step_name = odb.steps.keys()[-1]
    step = odb.steps[step_name]
    print(f"Step: {step_name}")
    print(f"Number of frames: {len(step.frames)}")

    # Extract history output
    history_regions = step.historyRegions
    print(f"History regions: {history_regions.keys()}")

    # Get roof center displacement (node 11, U1)
    # Get iso top displacements (nodes 1, 2, 3, U1)
    # Try to get from history output first
    history_data = {}
    for region_name, region in history_regions.items():
        for var_name, var_data in region.historyOutputs.items():
            key = f"{region_name}_{var_name}"
            data = [(t, v) for t, v in var_data.data]
            history_data[key] = data
            print(f"  {key}: {len(data)} points, range={min(v for _,v in data):.4f} to {max(v for _,v in data):.4f}")

    # Save to JSON
    with open('abaqus_history.json', 'w') as f:
        json.dump({k: [(float(t), float(v)) for t,v in vals] for k, vals in history_data.items()}, f, indent=2)

    # Also extract field output for all nodes at key frames
    last_frame = step.frames[-1]
    print(f"\nLast frame time: {last_frame.frameValue}")

    field_outputs = last_frame.fieldOutputs
    print(f"Field outputs: {field_outputs.keys()}")

    if 'U' in field_outputs:
        u_field = field_outputs['U']
        node_data = {}
        for val in u_field.values:
            node_data[int(val.nodeLabel)] = {
                'U1': float(val.data[0]),
                'U2': float(val.data[1]) if len(val.data) > 1 else 0.0,
            }

        with open('abaqus_final_displacements.json', 'w') as f:
            json.dump(node_data, f, indent=2)

        # Print key nodes
        for node_id in [1, 2, 3, 11]:
            if node_id in node_data:
                print(f"  Node {node_id}: U1={node_data[node_id]['U1']:.6f}, U2={node_data[node_id]['U2']:.6e}")

    # Extract time history from field output (all frames)
    print("\nExtracting full time histories...")
    n_frames = len(step.frames)
    times = []
    roof_ux = []
    iso1_ux = []
    iso2_ux = []
    iso3_ux = []

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

        if (i+1) % 5000 == 0:
            print(f"  Frame {i+1}/{n_frames}: t={t:.3f}s")

    # Save time history
    output = {
        'software': 'ABAQUS',
        'max_time_s': max(times),
        'n_frames': n_frames,
        'times': times,
        'roof_ux_m': roof_ux,
        'iso1_ux_m': iso1_ux,
        'iso2_ux_m': iso2_ux,
        'iso3_ux_m': iso3_ux,
        'peak_roof_displacement_m': max(abs(r) for r in roof_ux),
        'peak_isolation_displacement_m': max(
            max(abs(v) for v in iso1_ux),
            max(abs(v) for v in iso2_ux),
            max(abs(v) for v in iso3_ux)
        ),
    }

    with open('abaqus_time_history.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nPeak roof displacement (0-{max(times):.1f}s): {output['peak_roof_displacement_m']:.6f} m")
    print(f"Peak isolation displacement (0-{max(times):.1f}s): {output['peak_isolation_displacement_m']:.6f} m")
    print(f"Results saved to abaqus_time_history.json and abaqus_history.json")

    odb.close()

except ImportError as e:
    print(f"ERROR: Cannot import ABAQUS modules. Must run with 'abaqus python'. {e}")
    print("Creating diagnostic output instead.")

    # Check STA file
    sta_file = 'frame_bw_dynamic.sta'
    if os.path.exists(sta_file):
        with open(sta_file) as f:
            lines = f.readlines()
        last_line = lines[-3] if len(lines) >= 3 else lines[-1]
        print(f"STA last line: {last_line.strip()}")

    sys.exit(1)
