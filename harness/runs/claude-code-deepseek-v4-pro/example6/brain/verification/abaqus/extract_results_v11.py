"""Extract v11 ODB results."""
import os, json

odb_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(odb_dir)

try:
    from odbAccess import openOdb

    odb = openOdb('frame_bw_dynamic_v11.odb', readOnly=True)
    step_name = odb.steps.keys()[-1]
    step = odb.steps[step_name]
    n_frames = len(step.frames)
    print(f"Frames: {n_frames}")

    times, roof_ux, iso1_ux, iso2_ux, iso3_ux = [], [], [], [], []
    for frame in step.frames:
        t = frame.frameValue
        u_field = frame.fieldOutputs['U']
        u_dict = {}
        for val in u_field.values:
            u_dict[int(val.nodeLabel)] = float(val.data[0])
        times.append(t)
        roof_ux.append(u_dict.get(11, 0.0))
        iso1_ux.append(u_dict.get(1, 0.0))
        iso2_ux.append(u_dict.get(2, 0.0))
        iso3_ux.append(u_dict.get(3, 0.0))

    peak_roof = max(abs(r) for r in roof_ux)
    peak_iso = max(max(abs(v) for v in iso1_ux), max(abs(v) for v in iso2_ux), max(abs(v) for v in iso3_ux))

    print(f"Time range: 0 to {max(times):.2f}s")
    print(f"Peak roof displacement: {peak_roof:.6f} m")
    print(f"Peak iso displacement: {peak_iso:.6f} m")

    # Print at selected times
    for i, t in enumerate(times):
        if t >= 1.0 and t <= 15.0 and abs(round(t) - t) < 0.001:
            print(f"  t={t:.1f}s: roof={roof_ux[i]:.6f}, iso1={iso1_ux[i]:.6f}")

    output = {
        'software': 'ABAQUS v11',
        'max_time_s': max(times),
        'n_frames': n_frames,
        'times': times,
        'roof_ux_m': roof_ux,
        'iso1_ux_m': iso1_ux,
        'iso2_ux_m': iso2_ux,
        'iso3_ux_m': iso3_ux,
        'peak_roof_displacement_m': peak_roof,
        'peak_isolation_displacement_m': peak_iso,
    }
    with open('abaqus_time_history_v11.json', 'w') as f:
        json.dump(output, f, indent=2)

    odb.close()
except ImportError as e:
    print(f"ERROR: {e}")
