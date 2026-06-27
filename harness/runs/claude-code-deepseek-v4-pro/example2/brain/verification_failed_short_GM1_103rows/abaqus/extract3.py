from odbAccess import *
from abaqusConstants import *

odb = openOdb(path='test_conn4.odb')
step = odb.steps['Step-1']
frames = step.frames
print(f'Frames: {len(frames)}')
for i, frame in enumerate(frames):
    t = frame.frameValue
    for v in frame.fieldOutputs['U'].values:
        if v.nodeLabel == 2:
            print(f'Frame {i} t={t:.4f}: U1={v.data[0]:.8e}')

# Check history output
for region in step.historyRegions.values():
    print(f'History region: {region.name}')
    for hist in region.historyOutputs.values():
        print(f'  {hist.name}: {len(hist.data)} points')
        for t, v in hist.data:
            print(f'    t={t:.4f}: {v:.8e}')
odb.close()
