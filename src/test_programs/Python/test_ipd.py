import numpy as np
import pyopenmmpol as ommp
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--verbose", "-v",
                    help="Increase verbosity level to debug",
                    action="store_true")
parser.add_argument("--mmpol", "-i",
                    help=".mmpol file to be used as input",
                    required=True)
parser.add_argument("--electric-field", "-f",
                    help="""File containing a formatted matrix
                         (3xN) with N number of induced dipoles,
                         containing the electric field at
                         polarizable sites""")
parser.add_argument("--out-file", "-o",
                    help="""Formatted file where to print the
                            results""", default=None)
#parser.add_argument("--solver", "-s",
#                    help="""Solver method for the linear
#                            system""",
#                    default='default',
#                    choices=ommp.available_solvers.keys())

args = parser.parse_args()

ext_field_file = args.electric_field
infile = args.mmpol
outfile = args.out_file
#solver = args.solver

if(args.verbose):
    ommp.set_verbose(3)
else:
    ommp.set_verbose(1)

my_system = ommp.OMMPSystem(infile)

if ext_field_file is not None:
    ef = np.loadtxt(ext_field_file)
else:
    ef = np.zeros((my_system.pol_atoms, 3))

my_system.set_external_field(ef)
ipd = my_system.ipd

if outfile is None:
    for k in range(ipd.shape[0]):
        for i in range(ipd.shape[1]):
            print('{:20.12f}{:20.12f}{:20.12f}'.format(ipd[k,i,0],
                                                       ipd[k,i,1],
                                                       ipd[k,i,2]))
else:
    with open(outfile, 'w') as f:
        for k in range(ipd.shape[0]):
            for i in range(ipd.shape[1]):
                print('{:20.12f}{:20.12f}{:20.12f}'.format(ipd[k,i,0],
                                                           ipd[k,i,1],
                                                           ipd[k,i,2]),
                      file=f)

