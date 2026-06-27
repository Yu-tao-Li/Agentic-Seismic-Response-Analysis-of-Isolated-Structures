from odbAccess import openOdb

odb = openOdb("example6_boucwen.odb", readOnly=True)
print("steps", list(odb.steps.keys()))
for step_name, step in odb.steps.items():
    print("step", step_name, "frames", len(step.frames), "history_regions", len(step.historyRegions))
    if step.frames:
        print("frame0_fields", sorted(step.frames[0].fieldOutputs.keys()))
        print("last_fields", sorted(step.frames[-1].fieldOutputs.keys()))
    print("history_keys", list(step.historyRegions.keys())[:50])
odb.close()
