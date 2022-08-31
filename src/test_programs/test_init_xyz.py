import numpy as np
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
if(args.verbose):
    ommp.set_verbose(3)
else:
    ommp.set_verbose(1)

ommp.init_xyz(inprm, inxyz)

if outfile is None:
    ommp.print_summary()
else:
    ommp.print_summary_to_file(outfile)

ommp.terminate()

