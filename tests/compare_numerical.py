import numpy as np
import sys

A = np.loadtxt(sys.argv[1])
B = np.loadtxt(sys.argv[2])
if len(sys.argv) > 3:
    rtol = float(sys.argv[3])
else:
    rtol = 1e-05

if len(sys.argv) > 4:
    atol = float(sys.argv[4])
else:
    atol = 1e-08

print(rtol, atol)

if(np.allclose(A, B, rtol=rtol, atol=atol)):
    quit(0)
else:
    quit(1)
