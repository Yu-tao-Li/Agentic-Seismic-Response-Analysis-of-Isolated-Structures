"""Debug ODB structure."""
import sys
sys.path.insert(0, r'<abaqus_install>\EstProducts\2025\win_b64\code\python\lib')
from odbAccess import openOdb

odb = openOdb('isolated_building.odb')
print("Steps:", list(odb.steps.keys()))
for name, step in odb.steps.items():
    print(f"Step '{name}': {len(step.frames)} frames, time={step.totalTime}")
    # Check history regions
    print(f"  History regions: {list(step.historyRegions.keys())}")
    for hr_name, hr in step.historyRegions.items():
        print(f"    {hr_name}: outputs={list(hr.historyOutputs.keys())}")
        for ho_name, ho in hr.historyOutputs.items():
            print(f"      {ho_name}: {len(ho.data)} data points")

odb.close()
