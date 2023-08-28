import json
import os.path as path

def generate_test(jsonfile, program, ref, ef, fout):
    basename = None
    with open(jsonfile, "r") as f:
        data = json.loads(f.read())
    basename = data['name']

    if program == "init":
        tname = "{:s}_init".format(basename)
        tout = "{:s}.out".format(tname)
        print("""add_test(NAME {:s} 
                          COMMAND bin/${{TESTLANG}}_test_SI_init 
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s} 
                          Testing/{:s})""".format(tname, jsonfile, tout),
              file=fout)
        print("""add_test(NAME {:s}_comp
                          COMMAND ${{CMAKE_COMMAND}} -E compare_files
                          Testing/{:s}
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s})""".format(tname, tout, ref),
              file=fout)
    elif program == "energy":
        tname = "{:s}_energy".format(basename)
        if ef.lower() != "none":
            tname += '_{:s}'.format(path.basename(ef)[:-4])
            ef_str = 'tests/'+ef
        else:
            ef_str = ''
        tout = "{:s}.out".format(tname)

        print("""add_test(NAME {:s} 
                          COMMAND bin/${{TESTLANG}}_test_SI_potential 
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                          Testing/{:s} {:s})""".format(tname, jsonfile, tout, ef_str),
              file=fout)
        print("""add_test(NAME {:s}_comp 
                          COMMAND python3 ${{CMAKE_SOURCE_DIR}}/tests/compare_potential.py
                          Testing/{:s}
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s} 
                          {:6.5g} {:6.5g})""".format(tname, tout, ref, 1e-5, 1e-6 ),
              file=fout)
    elif program == "ipd":
        tname = "{:s}_ipd".format(basename)
        if ef.lower() != "none":
            tname += '_{:s}'.format(path.basename(ef)[:-4])
            ef_str = 'tests/'+ef
        else:
            ef_str = ''
        tout = "{:s}.out".format(tname)

        print("""add_test(NAME {:s} 
                          COMMAND bin/${{TESTLANG}}_test_SI_potential 
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                          Testing/{:s} {:s})""".format(tname, jsonfile, tout, ef_str),
              file=fout)
        print("""add_test(NAME {:s}_comp 
                          COMMAND python3 ${{CMAKE_SOURCE_DIR}}/tests/compare_ipd.py
                          Testing/{:s}
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s} 
                          {:6.5g} {:6.5g})""".format(tname, tout, ref, 1e-5, 1e-6 ),
              file=fout)
    #else if program == "ipd":
    #    ...
    #else if program == "grad":
    #    ...
    else:
        print("message(FATAL_ERROR, \"Automatically generated test {:s} cannot be understood\")".format(program), file=fout)

with open("test_list") as f, \
     open("TestsCmake.txt", "w+") as fout:
    for l in f:
        if not l.startswith("#") and l.strip():
            jsonf, program, ref, ef = l.split()
            generate_test(jsonf, program, ref, ef, fout)
