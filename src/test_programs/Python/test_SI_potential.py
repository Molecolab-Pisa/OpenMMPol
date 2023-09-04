import pyopenmmpol as ommp
import sys
import json

if len(sys.argv) != 4 and len(sys.argv) != 3:
    print("""Syntax expected"""
          """    $ test_init_xyz.exe <JSON FILE> <OUTPUT FILE> [<EF_FILE>]""")
    exit(1)

my_system, my_qmh = ommp.smartinput(sys.argv[1])
ommp.set_outputfile(sys.argv[2])

if len(sys.argv) == 4:
    ext_field_file = sys.argv[3]
else:
    ext_field_file = None

if ext_field_file is not None:
    ef = np.loadtxt(ext_field_file)
    my_system.set_external_field(ef, solver)

fake_qm = None
if my_qmh:
    if my_system.use_linkatoms:
        # Initialize a fake qm
        with open(sys.argv[1], 'r') as jsonFile:
            values = json.load(jsonFile)
        try:
            prm_file = values['qm']['prm_file']['path']
        except:
            prm_file = ""

        fake_qm = ommp.OMMPSystem(my_qmh, prm_file)
        fake_qm.turn_pol_off(range(fake_qm.mm_atoms)) # Turn all pol. off.

em =   my_system.get_fixedelec_energy()
ep =   my_system.get_polelec_energy()
ev =   my_system.get_vdw_energy()
if(my_qmh):
    evqmmm = ommp.qm_helper_vdw_energy(my_qmh, my_system)
else:
    evqmmm = 0.
eb =   my_system.get_bond_energy()
ea =   my_system.get_angle_energy()
eba =  my_system.get_strbnd_energy()
eub =  my_system.get_urey_energy()
eopb = my_system.get_opb_energy()
ept =  my_system.get_pitors_energy()
et =   my_system.get_torsion_energy()
ett =  my_system.get_tortor_energy()
eat =  my_system.get_angtor_energy()
ebt =  my_system.get_strtor_energy()
eit =  my_system.get_imptorsion_energy()
etotmm = my_system.get_full_energy()
if(fake_qm):
    eqm = fake_qm.get_full_energy()
else:
    eqm = 0.0
etot = etotmm + eqm + evqmmm;

eaa = 0.0
eopd = 0.0
eid = 0.0
er = 0.0
edsp = 0.0
ec = 0.0
ecd = 0.0
ed = 0.0
ect = 0.0
erxf = 0.0
es = 0.0
elf = 0.0
eg = 0.0
ex = 0.0

ommp.message("EM      {:20.12e}".format(em), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EP      {:20.12e}".format(ep), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EV      {:20.12e}".format(ev), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EVQMMM  {:20.12e}".format(evqmmm), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EB      {:20.12e}".format(eb), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EA      {:20.12e}".format(ea), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EBA     {:20.12e}".format(eba), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EUB     {:20.12e}".format(eub), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EOPB    {:20.12e}".format(eopb), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EPT     {:20.12e}".format(ept), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ET      {:20.12e}".format(et), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ETT     {:20.12e}".format(ett), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EAA     {:20.12e}".format(eaa), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EOPD    {:20.12e}".format(eopd), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EID     {:20.12e}".format(eid), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EIT     {:20.12e}".format(eit), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EBT     {:20.12e}".format(ebt), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EAT     {:20.12e}".format(eat), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ER      {:20.12e}".format(er), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EDSP    {:20.12e}".format(edsp), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EC      {:20.12e}".format(ec), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ECD     {:20.12e}".format(ecd), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ED      {:20.12e}".format(ed), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ECT     {:20.12e}".format(ect), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ERXF    {:20.12e}".format(erxf), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ES      {:20.12e}".format(es), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ELF     {:20.12e}".format(elf), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EG      {:20.12e}".format(eg), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EX      {:20.12e}".format(ex), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ETOTMM  {:20.12e}".format(etotmm), ommp.verbosity['none'], "TEST-ENE")
ommp.message("EQM     {:20.12e}".format(eqm), ommp.verbosity['none'], "TEST-ENE")
ommp.message("ETOT    {:20.12e}".format(etot), ommp.verbosity['none'], "TEST-ENE")

for i in range(my_system.n_ipd):
    for j in range(my_system.pol_atoms):
        ommp.message("%20.12e %20.12e %20.12e" % tuple(my_system.ipd[i,j,:]), ommp.verbosity['none'], "TEST-IPD")

if my_qmh:
    del my_qmh
del my_system
