"""Small ODB diagnostic for Example 5 story yielding."""
from odbAccess import openOdb
import numpy as np

odb = openOdb('example5_model.odb')
step = odb.steps['DYNAMIC']
inst = odb.rootAssembly.instances['PART-1-1']

for elset_name in ['SPRING1', 'SPRING2', 'SPRING3', 'SPRING4']:
    region = inst.elementSets[elset_name]
    stress = []
    peeq = []
    for frame in step.frames:
        s_data = frame.fieldOutputs['S'].getSubset(region=region).values[0].data
        pe_data = frame.fieldOutputs['PEEQ'].getSubset(region=region).values[0].data
        stress.append(s_data[0] if hasattr(s_data, '__len__') else s_data)
        peeq.append(pe_data[0] if hasattr(pe_data, '__len__') else pe_data)
    stress = np.asarray(stress)
    peeq = np.asarray(peeq)
    print(elset_name, 'stress_peak', np.max(np.abs(stress)), 'peeq_max', np.max(peeq))

odb.close()
