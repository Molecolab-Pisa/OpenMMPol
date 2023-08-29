import json
import os.path as path

def generate_test(jsonfile, program, ref, ef, fout, atol, rtol):
    atol_ene = 1e-6
    rtol_ene = 1e-6
    atol_grad = 1e-4
    rtol_grad = 1e-3
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
        if atol is None:
            atol = atol_ene
        if rtol is None:
            rtol = rtol_ene

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
                          {:6.5g} {:6.5g})""".format(tname, tout, ref, rtol, atol),
              file=fout)
    elif program == "ipd":
        if atol is None:
            atol = 1e-5
        if rtol is None:
            rtol = 1e-6
        
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
                          {:6.5g} {:6.5g})""".format(tname, tout, ref, rtol, atol),
              file=fout)
    elif program == "grad-num":
        if atol is None:
            atol = atol_grad
        if rtol is None:
            rtol = rtol_grad
        
        tname = "{:s}_geomgrad".format(basename)
        tname_num = tname+'_num'
        tout_num = "{:s}.out".format(tname_num)
        tname_ana = tname+'_ana'
        tout_ana = "{:s}.out".format(tname_ana)
        doref = False

        if ref.lower() != "none":
            doref = True

        print("""add_test(NAME {:s} 
                          COMMAND bin/${{TESTLANG}}_test_SI_geomgrad_num
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                          Testing/{:s})""".format(tname_num, jsonfile, tout_num),
              file=fout)
        print("""add_test(NAME {:s} 
                          COMMAND bin/${{TESTLANG}}_test_SI_geomgrad
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                          Testing/{:s})""".format(tname_ana, jsonfile, tout_ana),
              file=fout)
        print("""add_test(NAME {:s}_comp_num_ana 
                          COMMAND python3 ${{CMAKE_SOURCE_DIR}}/tests/compare_geomgrad.py
                          Testing/{:s}
                          Testing/{:s} 
                          {:6.5g} {:6.5g})""".format(tname, tout_num, tout_ana, rtol, atol),
              file=fout)
        if doref:
            print("""add_test(NAME {:s}_comp_ana_ref 
                            COMMAND python3 ${{CMAKE_SOURCE_DIR}}/tests/compare_geomgrad.py
                            Testing/{:s}
                            ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                            {:6.5g} {:6.5g})""".format(tname, tout_ana, ref, rtol, atol),
                file=fout)
            print("""add_test(NAME {:s}_comp_num_ref 
                            COMMAND python3 ${{CMAKE_SOURCE_DIR}}/tests/compare_geomgrad.py
                            Testing/{:s}
                            ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                            {:6.5g} {:6.5g})""".format(tname, tout_num, ref, rtol, atol),
                file=fout)
    elif program == "grad":
        if atol is None:
            atol = atol_grad
        if rtol is None:
            rtol = rtol_grad
        
        tname = "{:s}_geomgrad".format(basename)
        tname_ana = tname+'_ana'
        tout_ana = "{:s}.out".format(tname_ana)
        doref = False

        print("""add_test(NAME {:s} 
                          COMMAND bin/${{TESTLANG}}_test_SI_geomgrad
                          ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                          Testing/{:s})""".format(tname_ana, jsonfile, tout_ana),
              file=fout)
        print("""add_test(NAME {:s}_comp_ana_ref 
                        COMMAND python3 ${{CMAKE_SOURCE_DIR}}/tests/compare_geomgrad.py
                        Testing/{:s}
                        ${{CMAKE_SOURCE_DIR}}/tests/{:s}
                        {:6.5g} {:6.5g})""".format(tname, tout_ana, ref, rtol, atol),
            file=fout)
    else:
        print("message(FATAL_ERROR, \"Automatically generated test {:s} cannot be understood\")".format(program), file=fout)

with open("test_list") as f, \
     open("TestsCmake.txt", "w+") as fout:
    for l in f:
        if not l.startswith("#") and l.strip():
            tok = l.split()
            jsonf, program, ref, ef = tok[:4]
            if len(tok) >= 5:
                atol = float(tok[4])
            else:
                atol = None
            if len(tok) >= 6:
                rtol = float(tok[6])
            else:
                rtol = None
            generate_test(jsonf, program, ref, ef, fout, atol, rtol)
