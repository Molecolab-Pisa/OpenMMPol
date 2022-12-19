import pyopenmmpol as ommp
import sys

if len(sys.argv) not in [1, 2]:
    print("Syntax expected")
    print("    $ test_init.exe <INPUT FILE> [<OUTPUT FILE>]")
    exit(0)
elif len(sys.argv) == 1:
    outfile = None
else:
    outfile = sys.argv[1]

my_system = ommp.OMMPSystem(sys.argv[0])

if outfile is not None:
    my_system.print_summary(outfile)
else:
    my_system.print_summary()

del my_system
