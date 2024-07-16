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
parser.add_argument("--out-file", "-o",
                    help="""Formatted file where to print the
                            results""", default=None)

args = parser.parse_args()

inxyz = args.xyz
inprm = args.prm
outfile = args.out_file
if(args.verbose):
    ommp.set_verbose(3)
else:
    ommp.set_verbose(1)

my_system = ommp.OMMPSystem(inxyz, inprm)

if outfile is None:
    my_system.print_summary()
else:
    my_system.print_summary_to_file(outfile)

del my_system
