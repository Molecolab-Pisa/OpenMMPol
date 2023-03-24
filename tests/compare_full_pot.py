import numpy as np
import sys

def get_term(fname, term):
    a = None
    with open(fname, 'r') as f:
        for l in f:
            if l.startswith("{:s} ".format(term)):
                a = float(l.split()[-1])
    return a

a = get_term(sys.argv[1], sys.argv[3])
b = get_term(sys.argv[2], sys.argv[3])

if len(sys.argv) > 5:
    rtol = float(sys.argv[4])
else:
    rtol = 1e-05

if len(sys.argv) > 6:
    atol = float(sys.argv[5])
else:
    atol = 1e-08

if abs(b) > 0:
    print(a, b, abs(a-b)/abs(b), abs(a-b))
else:
    print(a, b, 0.00, abs(a-b))

if np.isclose(a, b, rtol, atol):
    exit(0)
else:
    exit(1)
