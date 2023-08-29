import numpy as np
import sys

def get_avail_terms(fname):
    tnames = []
    OMMP_out = False
    with open(fname, 'r') as f:
        for l in f:
            if l.startswith('[OMMP]') and not OMMP_out:
                # This file is the output of an OMMP program
                OMMP_out = True
            if OMMP_out and l.startswith("[OMMP]  [TEST-GRD] Grad "):
                tnames += [l.split()[3]]
            elif not OMMP_out and l.startswith("Grad "):
                tnames += [l.split()[1]]
    return tnames

def get_term(fname, term):
    a = None
    pres = ""
    with open(fname, 'r') as f:
        for l in f:
            if l.startswith('[OMMP]') and pres == '':
                pres = "[OMMP]  [TEST-GRD] "
            if l.startswith("{:s}Grad {:s}".format(pres, term)) and l[len("{:s}Grad {:s}".format(pres, term))] in ' \n':
                a = []
                l = next(f)
                offset = len(pres.split())
                try:
                    float(l.split()[offset])
                except:
                    offset += 1

                while l.strip():
                    a += [[float(x) for x in l.split()[offset:]]]
                    l = next(f)
                    if l.strip() == pres.strip():
                        return np.array(a)

tnames_a = get_avail_terms(sys.argv[1])
tnames_b = get_avail_terms(sys.argv[2])
if len(tnames_a) == 0 or len(tnames_b) == 0:
    print("No gradients terms to compare!")
    exit(1)

if len(sys.argv) >= 4:
    rtol = float(sys.argv[3])
else:
    rtol = 1e-05

if len(sys.argv) >= 5:
    atol = float(sys.argv[4])
else:
    atol = 1e-08

outval = 0
outerr = []

for tname in list(set(tnames_a) & set(tnames_b)):
    a = get_term(sys.argv[1], tname)
    b = get_term(sys.argv[2], tname)
    
    isok = np.isclose(a, b, rtol, atol).all()
    delta = np.abs(a - b)
    deltar = np.abs((a-b)/b)
    m = np.max(delta)
    mr = 100*np.nanmax(deltar)
    if np.isnan(mr):
        mr = 0.0

    if not isok:
        print("*** ", end='')
        outval += 1
    else:
        print("    ", end='')

    print("GRAD {:s}".format(tname))
    for i in range(a.shape[0]):
        print("[{:5d}] {:12.6f} {:12.6f} {:12.6f}".format(i, a[i,0], a[i,1], a[i,2]))
        print("        {:12.6f} {:12.6f} {:12.6f}".format(b[i,0], b[i,1], b[i,2]))
        print("        {:12.6f} {:12.6f} {:12.6f}".format(delta[i,0], delta[i,1], delta[i,2]))
    out = "Delta Max {:5s} {:6d}/{:1d} {:8.2g}({:8.2g}) {:8.2g}%({:8.2g}%) {}".format(tname, np.argmax(delta)//3, np.argmax(delta)%3, m, atol, mr, rtol*100, isok)

    if not isok:
        outerr += ['[{:8s}] {:s}'.format(tname, out)]

for i in outerr:
    print(i)

exit(outval)

