import numpy as np
import sys

def get_avail_terms(fname):
    tnames = []
    OMMP_out = False
    with open(fname, 'r') as f:
        for l in f:
            if not l.strip():
                break
            if l.startswith('[OMMP]') and not OMMP_out:
                # This file is the output of an OMMP program
                OMMP_out = True
            if OMMP_out and l.startswith("[OMMP]  [TEST-ENE] "):
                tnames += [l.split()[2]]
            elif not OMMP_out:
                tnames += [l.split()[0]]
    return tnames

def get_term(fname, term):
    a = None
    pres = ""
    with open(fname, 'r') as f:
        for l in f:
            if l.startswith('[OMMP]') and pres == '':
                # This file is the output of an OMMP program
                pres = "[OMMP]  [TEST-ENE] "
            if l.startswith("{:s}{:s} ".format(pres, term)):
                a = float(l.split()[-1])
    return a

tnames_a = get_avail_terms(sys.argv[1])
tnames_b = get_avail_terms(sys.argv[2])
if len(tnames_a) == 0 or len(tnames_b) == 0:
    print("No potential terms to compare!")
    exit(1)

if len(sys.argv) > 3:
    rtol = float(sys.argv[3])
else:
    rtol = 1e-05

if len(sys.argv) > 4:
    atol = float(sys.argv[4])
else:
    atol = 1e-08

outval = 0
fmt = "{:<10s} {:20.10g} {:20.10g} {:12.6g}/{:12.6g} {:12.6g}/{:12.6g}"
for tname in list(set(tnames_a) & set(tnames_b)):
    a = get_term(sys.argv[1], tname)
    b = get_term(sys.argv[2], tname)
    isok = np.isclose(a, b, rtol, atol)

    if not isok:
        print("*** ", end='')
        outval += 1
    else:
        print("    ", end='')
    
    if abs(b) > 0:
        print(fmt.format(tname, a, b, abs(a-b)/abs(b), rtol, abs(a-b), atol))
    else:
        print(fmt.format(tname, a, b, 0.00, rtol, abs(a-b), atol))

exit(outval)
