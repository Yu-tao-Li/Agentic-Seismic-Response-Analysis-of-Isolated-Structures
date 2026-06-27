from odbAccess import *
from abaqusConstants import *
odb = openOdb(path='test_xaxis.odb')
step = odb.steps['Step-1']
frames = step.frames
print(f'Frames: {len(frames)}')
for i, frame in enumerate(frames):
    t = frame.frameValue
    for v in frame.fieldOutputs['U'].values:
        if v.nodeLabel == 2:
            print(f'Frame {i} t={t:.4f}: U1={v.data[0]:.8e}')
odb.close()
