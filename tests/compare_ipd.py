import numpy as np
import sys

def get_ipd(fname):
    a = []
    pres = ""
    with open(fname, 'r') as f:
        for l in f:
            if l.startswith('[OMMP]') and pres == '':
                # This file is the output of an OMMP program
                pres = "[OMMP]  [TEST-IPD] "
            if not pres or l.startswith("{:s}".format(pres)):
                if not pres:
                    dat = l
                else:
                    dat = l[len(pres):]
                a += [[float(x) for x in dat.split()]]
    return np.array(a)

A = get_ipd(sys.argv[1])
B = get_ipd(sys.argv[2])

if len(sys.argv) > 2:
    rtol = float(sys.argv[3])
else:
    rtol = 1e-05

if len(sys.argv) > 3:
    atol = float(sys.argv[4])
else:
    atol = 1e-08

#print(rtol, atol)

if(np.allclose(A, B, rtol=rtol, atol=atol)):
    quit(0)
else:
    quit(1)
