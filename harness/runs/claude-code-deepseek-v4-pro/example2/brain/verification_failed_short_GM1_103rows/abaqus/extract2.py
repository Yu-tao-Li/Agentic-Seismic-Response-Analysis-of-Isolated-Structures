from odbAccess import *
from abaqusConstants import *
import sys

odb = openOdb(path='test_conn3.odb')
step = odb.steps['Step-1']
regions = step.historyRegions
print('History regions:', regions.keys())

# Check field output frames
frames = step.frames
print(f'Number of frames in step: {len(frames)}')
for i, frame in enumerate(frames):
    print(f'Frame {i}: t={frame.frameValue:.6f}')

# Try field output
for i, frame in enumerate(frames):
    t = frame.frameValue
    try:
        uField = frame.fieldOutputs['U']
        nodeLabels = []
        for v in uField.values:
            nodeLabels.append(v.nodeLabel)
            if v.nodeLabel == 2:
                print(f'Frame {i} t={t:.4f}: Node 2 U1={v.data[0]:.8e}')
    except:
        print(f'Frame {i} t={t:.4f}: no U field')

odb.close()
