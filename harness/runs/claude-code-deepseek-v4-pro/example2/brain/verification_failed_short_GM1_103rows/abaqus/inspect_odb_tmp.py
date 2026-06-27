from odbAccess import openOdb
import json, os
odb = openOdb(path='three_story_shear.odb')
print('steps', list(odb.steps.keys()))
for sname, step in odb.steps.items():
    print('STEP', sname, 'frames', len(step.frames), 'historyRegions', len(step.historyRegions))
    print('history keys', list(step.historyRegions.keys())[:20])
    for rname, region in list(step.historyRegions.items())[:20]:
        print('REGION', rname, 'outputs', list(region.historyOutputs.keys()))
        for oname, hist in region.historyOutputs.items():
            print(' ', oname, 'n', len(hist.data), 'firstlast', (hist.data[0] if hist.data else None), (hist.data[-1] if hist.data else None))
    # field outputs
    if step.frames:
        f = step.frames[-1]
        print('last frame', f.frameValue, 'field keys', list(f.fieldOutputs.keys()))
        for key in ['U','A','V']:
            if key in f.fieldOutputs:
                vals=[]
                for v in f.fieldOutputs[key].values:
                    vals.append((getattr(v,'nodeLabel',None), v.data))
                print(key, vals[:20])
odb.close()
