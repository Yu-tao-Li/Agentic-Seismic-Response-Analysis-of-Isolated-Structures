from odbAccess import openOdb

odb = openOdb("three_story_shear_spring_mass.odb", readOnly=True)
try:
    step = odb.steps["FREQUENCY"]
    for index, frame in enumerate(step.frames):
        print("FRAME", index)
        print("  mode:", getattr(frame, "mode", None))
        print("  frequency repr:", repr(getattr(frame, "frequency", None)))
        print("  description:", getattr(frame, "description", ""))
finally:
    odb.close()
