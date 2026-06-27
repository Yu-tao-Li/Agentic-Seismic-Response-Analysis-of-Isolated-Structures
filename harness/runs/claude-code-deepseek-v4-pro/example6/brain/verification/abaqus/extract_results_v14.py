import os, json, sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))
try:
    from odbAccess import openOdb
    odb = openOdb('frame_bw_dynamic_v14.odb', readOnly=True)
    step = odb.steps[odb.steps.keys()[-1]]
    frames = step.frames
    print(f"Frames: {len(frames)}")
    times = [f.frameValue for f in frames]
    roofs = []
    isos = [[], [], []]
    for f in frames:
        u = {}
        for v in f.fieldOutputs['U'].values:
            u[int(v.nodeLabel)] = float(v.data[0])
        roofs.append(u.get(11, 0))
        isos[0].append(u.get(1, 0))
        isos[1].append(u.get(2, 0))
        isos[2].append(u.get(3, 0))
    all_iso = isos[0] + isos[1] + isos[2]
    pk_roof = max(max(roofs), -min(roofs))
    pk_iso = max(max(all_iso), -min(all_iso))
    print(f"Time range: 0 to {max(times):.2f}s")
    print(f"Peak roof: {pk_roof:.6f} m")
    print(f"Peak iso:  {pk_iso:.6f} m")
    for i, t in enumerate(times):
        if t >= 1.0 and t <= 16.0 and abs(round(t)-t) < 0.001:
            print(f"  t={t:.0f}s: roof={roofs[i]:.6f}, iso1={isos[0][i]:.6f}")
    odb.close()
except ImportError as e:
    print(f"ERROR: {e}")
