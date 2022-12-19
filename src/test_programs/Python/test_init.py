import pyopenmmpol as ommp
import sys

print(sys.argv)
if len(sys.argv) not in [2, 3]:
    print("Syntax expected")
    print("    $ test_init.exe <INPUT FILE> [<OUTPUT FILE>]")
    exit(0)
elif len(sys.argv) == 2:
    outfile = None
else:
    outfile = sys.argv[2]

my_system = ommp.OMMPSystem(sys.argv[1])

if outfile is not None:
    my_system.print_summary(outfile)
else:
    my_system.print_summary()

del my_system
