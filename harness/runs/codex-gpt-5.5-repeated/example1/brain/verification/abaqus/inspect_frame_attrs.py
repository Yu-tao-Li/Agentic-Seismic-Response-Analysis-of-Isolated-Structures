from odbAccess import openOdb

odb = openOdb("three_story_shear_spring_mass.odb", readOnly=True)
try:
    frame = odb.steps["FREQUENCY"].frames[1]
    for name in dir(frame):
        if name.startswith("_"):
            continue
        try:
            value = getattr(frame, name)
        except Exception as exc:
            value = "<error %s>" % exc
        if name in ("mode", "frequency", "description", "frameValue", "incrementNumber") or isinstance(
            value, (int, float, str)
        ):
            print(name, repr(value))
finally:
    odb.close()
