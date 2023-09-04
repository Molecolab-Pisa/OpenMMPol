import pyopenmmpol as ommp
import sys
import json

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print("""Syntax expected"""
          """    $ test_init_xyz.exe <JSON FILE> [<OUTPUT FILE>]""")
    exit(1)

my_system, my_qmh = ommp.smartinput(sys.argv[1])
ommp.se~.

t
if len(sys.argv) == 3:
    my_system.print_summary(sys.argv[2])
else:
    my_system.print_summary()

if my_qmh:
    del my_qmh
del my_system
