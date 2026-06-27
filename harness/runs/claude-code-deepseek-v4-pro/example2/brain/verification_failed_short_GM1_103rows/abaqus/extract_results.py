from odbAccess import *
from abaqusConstants import *
import sys

odb = openOdb(path='test_conn3.odb')
step = odb.steps['Step-1']
lastFrame = step.frames[-1]
roofNode = odb.rootAssembly.instances['PART-1-1'].nodes
# Get displacement at node 2 (index 1 in 0-based)
dispField = lastFrame.fieldOutputs['U']
for node in roofNode:
    disp = dispField.getSubset(region=node).values[0]
    if node.label == 2:
        print(f'Node {node.label}: U1={disp.data[0]:.6e}, U2={disp.data[1]:.6e}, U3={disp.data[2]:.6e}')

# Get all frames
print(f'Number of frames: {len(step.frames)}')
for i, frame in enumerate(step.frames):
    t = frame.frameValue
    dispField = frame.fieldOutputs['U']
    for node in roofNode:
        if node.label == 2:
            disp = dispField.getSubset(region=node).values[0]
            print(f'Frame {i}: t={t:.4f}s, U1={disp.data[0]:.8e}')
            break

odb.close()
