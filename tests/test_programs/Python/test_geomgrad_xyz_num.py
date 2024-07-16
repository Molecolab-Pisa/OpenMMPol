import numpy as np
import scipy.constants as pc
import pyopenmmpol as ommp
import argparse
from random import sample as randsample

def num_ana_compare(sys, gradf, name, delt):
    ana = gradf(sys)
    num = gradf(sys, numerical=True)

    print("DELTA NUM - ANA {:s}".format(name))
    Mdelta = 0.0
    for i in range(ana.shape[0]):
        delta = num[i] - ana[i]
        if max(abs(delta)) > Mdelta:
            Mdelta = max(abs(delta))

        print("{:5d} (A) {:+12.8g} {:+12.8g} {:+12.8g}".format(i+1,
                                                               ana[i,0],
                                                               ana[i,1],
                                                               ana[i,2]))
        print("      (N) {:+12.8g} {:+12.8g} {:+12.8g}".format(num[i,0],
                                                               num[i,1],
                                                               num[i,2]))
        print("      (D) {:+12.8g} {:+12.8g} {:+12.8g}".format(delta[0],
                                                               delta[1],
                                                               delta[2]))
        print()

    if Mdelta > delt:
        print("WARNING {:s} delta ({:.3g}) > max_delta ({:.3g})".format(name, Mdelta, delt))
        return 1
    else:
        return 0

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
parser.add_argument("--out-file", "-o",
                    help="""Formatted file where to print the
                            results""", default=None)
parser.add_argument("--frac-frozen",
                    help="""Fraction of atoms that should be
                            randomly frozen""", default="0.0")
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
ffr = float(args.frac_frozen)
if(args.verbose):
    ommp.set_verbose(3)
else:
    ommp.set_verbose(0)

ms = ommp.OMMPSystem(inxyz, inprm)
frozen = randsample(range(ms.mm_atoms), int(ms.mm_atoms*ffr))
print("Atomi attivi: ", ms.mm_atoms - int(ms.mm_atoms*ffr))
ms.set_frozen_atoms(frozen)

rc = 0
rc += num_ana_compare(ms, ommp.OMMPSystem.fixedelec_geomgrad, "fixedelec", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.vdw_geomgrad, "vdw", 1e-8)
#rc += num_ana_compare(ms, ommp.OMMPSystem.polelec_geomgrad, "polelec", 1e-8)
rc += num_ana_compare(ms, ommp.OMMPSystem.bond_geomgrad, "bond", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.angle_geomgrad, "angle", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.urey_geomgrad, "urey", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.torsion_geomgrad, "torsion", 1e-8)
rc += num_ana_compare(ms, ommp.OMMPSystem.imptorsion_geomgrad, "imptorsion", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.strbnd_geomgrad, "strbnd", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.angtor_geomgrad, "angtor", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.opb_geomgrad, "opb", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.strtor_geomgrad, "strtor", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.tortor_geomgrad, "tortor", 1e-9)
rc += num_ana_compare(ms, ommp.OMMPSystem.pitors_geomgrad, "pitors", 1e-9)

del ms

exit(rc)

