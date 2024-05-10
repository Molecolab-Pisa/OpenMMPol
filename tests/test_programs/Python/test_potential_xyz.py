import numpy as np
import scipy.constants as pc
import pyopenmmpol as ommp
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--verbose", "-v",
                    help="Increase verbosity level to debug",
                    action="store_true")
parser.add_argument("--xyz",
                    help=".xyz file to be used as input",
                    required=True)
parser.add_argument("--prm",
                    help=".prm forcefield file to be used as input",
                    required=True)
parser.add_argument("--electric-field", "-f",
                    help="""File containing a formatted matrix (3xN)
                         with N number of induced dipoles,
                         containing the electric field at
                         polarizable sites""")
parser.add_argument("--out-file", "-o",
                    help="""Formatted file where to print the
                            results""", default=None)
parser.add_argument("--solver", "-s",
                    help="""Solver method for the linear
                            system""",
                    default='default',
                    choices=ommp.available_solvers.keys())

args = parser.parse_args()

inxyz = args.xyz
inprm = args.prm
outfile = args.out_file
solver = args.solver
if(args.verbose):
    ommp.set_verbose(3)
else:
    ommp.set_verbose(1)
ext_field_file = args.electric_field

ms = ommp.OMMPSystem(inxyz, inprm)
au2kcalmol = pc.physical_constants['Hartree energy'][0]/(1000.0 * pc.calorie / pc.N_A )

if ext_field_file is not None:
    ef = np.loadtxt(ext_field_file)
    ms.set_external_field(ef, solver)

em =   ms.get_fixedelec_energy() * au2kcalmol
ep =   ms.get_polelec_energy() * au2kcalmol
ev =   ms.get_vdw_energy() * au2kcalmol
eb =   ms.get_bond_energy() * au2kcalmol
ea =   ms.get_angle_energy() * au2kcalmol
eba =  ms.get_strbnd_energy() * au2kcalmol
eub =  ms.get_urey_energy() * au2kcalmol
eopb = ms.get_opb_energy() * au2kcalmol
ept =  ms.get_pitors_energy() * au2kcalmol
et =   ms.get_torsion_energy() * au2kcalmol
ett =  ms.get_tortor_energy() * au2kcalmol
eat =  ms.get_angtor_energy() * au2kcalmol
ebt =  ms.get_strtor_energy() * au2kcalmol
eit =   ms.get_imptorsion_energy() * au2kcalmol
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

print("EM      {:20.12e}".format(em))
print("EP      {:20.12e}".format(ep))
print("EV      {:20.12e}".format(ev))
print("EB      {:20.12e}".format(eb))
print("EA      {:20.12e}".format(ea))
print("EBA     {:20.12e}".format(eba))
print("EUB     {:20.12e}".format(eub))
print("EOPB    {:20.12e}".format(eopb))
print("EPT     {:20.12e}".format(ept))
print("ET      {:20.12e}".format(et))
print("ETT     {:20.12e}".format(ett))

print("EAA     {:20.12e}".format(eaa))
print("EOPD    {:20.12e}".format(eopd))
print("EID     {:20.12e}".format(eid))
print("EIT     {:20.12e}".format(eit))
print("EBT     {:20.12e}".format(ebt))
print("EAT     {:20.12e}".format(eat))
print("ER      {:20.12e}".format(er))
print("EDSP    {:20.12e}".format(edsp))
print("EC      {:20.12e}".format(ec))
print("ECD     {:20.12e}".format(ecd))
print("ED      {:20.12e}".format(ed))
print("ECT     {:20.12e}".format(ect))
print("ERXF    {:20.12e}".format(erxf))
print("ES      {:20.12e}".format(es))
print("ELF     {:20.12e}".format(elf))
print("EG      {:20.12e}".format(eg))
print("EX      {:20.12e}".format(ex))

del ms

