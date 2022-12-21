import pyopenmmpol as ommp
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--verbose", "-v",
                    help="Increase verbosity level to debug",
                    action="store_true")
parser.add_argument("--mmpol", "-i",
                    help=".mmpol file to be used as input",
                    required=True)
parser.add_argument("--out-file", "-o",
                    help="""Formatted file where to print the
                            results""", default=None)

args = parser.parse_args()

my_system = ommp.OMMPSystem(args.mmpol)

if outfile is not None:
    my_system.print_summary(args.out_file)
else:
    my_system.print_summary()

del my_system
