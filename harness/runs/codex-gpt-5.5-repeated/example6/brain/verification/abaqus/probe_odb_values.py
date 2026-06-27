from odbAccess import openOdb

odb = openOdb("example6_boucwen.odb", readOnly=True)
asm = odb.rootAssembly
print("instances", list(asm.instances.keys()))
print("node_sets", list(asm.nodeSets.keys())[:50])
print("element_sets", list(asm.elementSets.keys())[:50])
step = odb.steps["DYNAMIC_BOUCWEN"]
for idx in [0, 1, -1]:
    frame = step.frames[idx]
    print("frame", idx, "value", frame.frameValue)
    for field_name in ["U", "RF", "CF"]:
        field = frame.fieldOutputs[field_name]
        vals = []
        for v in field.values:
            if v.nodeLabel in [11, 12, 13, 101, 102, 103, 112, 122, 132]:
                vals.append((v.nodeLabel, tuple(v.data)))
        print(field_name, vals[:20])
odb.close()
