"""Explore ODB structure to find correct API for frequency extraction."""
from abaqus import session
import json

odb = session.openOdb(name='three_story_modal.odb')
step = odb.steps['FREQ']

output = []
output.append(f"Step name: {step.name}")
output.append(f"Number of frames: {len(step.frames)}")

for i, frame in enumerate(step.frames):
    output.append(f"\nFrame {i}:")
    output.append(f"  frameId: {frame.frameId}")
    output.append(f"  description: {frame.description}")
    output.append(f"  domain: {frame.domain}")

    # Check available attributes
    attrs = [a for a in dir(frame) if not a.startswith('_')]
    output.append(f"  attributes: {attrs}")

    # Check field outputs
    field_names = [f.name for f in frame.fieldOutputs.values()]
    output.append(f"  fieldOutputs: {field_names}")

    if i > 3:
        break

with open('explore_odb_output.txt', 'w') as f:
    f.write('\n'.join(output))
