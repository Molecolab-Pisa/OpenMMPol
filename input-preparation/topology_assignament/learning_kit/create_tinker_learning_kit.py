#  Standard Amino Acids:  GLY, ALA, VAL, LEU, ILE, SER, THR, CYS, CYX, PRO,
#  PHE, TYR, TRP, HIS, ASP, ASN, GLU, GLN, MET, LYS, ARG, ORN, AIB
#
#  Alternative Protonation States:  CYD, TYD, HID, HIE, HIP, ASH, GLH, LYD
#
#  N-Terminal Cap Residues:  H2N=Deprotonated, FOR=Formyl, ACE=Acetyl,
#                            PCA=Pyroglutamic Acid
#  C-Terminal Cap Residues:  COH=Protonated, NH2=Amide, NME=N-MethylAmide


aminoacids = ['GLY',
              'ALA',
              'VAL',
              'LEU',
              'ILE',
              'SER',
              'THR',
              'CYS',
              'CYX',
              'PRO',
              'PHE',
              'TYR',
              'TRP',
              'HIS',
              'ASP',
              'ASN',
              'GLU',
              'GLN',
              'MET',
              'LYS',
              'ARG']

alt_aminoacids = ['CYD', 'TYD', 'HID', 'HIE', 'HIP', 'ASH', 'GLH', 'LYD']

c_term = ['H2N', 'FOR', 'ACE', 'PCA']

n_term = ['COH', 'NH2', 'NME']

ff_file = 'amoebabio18.prm'

for a in aminoacids + alt_aminoacids:
    with open("internal_{:s}.txt".format(a), 'w') as f:
        print("internal_{:s}".format(a), file=f)
        print("Internal {:s}".format(a), file=f)
        print("{:s}".format(ff_file), file=f)
        if a == 'CYX':
            print('CYX 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4', file=f)
            print('CYX 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3', file=f)
            print('CYX 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2', file=f)
            print('CYX 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1', file=f)
        else:
            print(a, file=f) # Default C-Terminal (COO-)
            print(a, file=f) # Internal
            print(a, file=f) # Default N-terminal (NH3+)
        print("\nN", file=f)

for ct in c_term:
    for a in aminoacids + alt_aminoacids:
        with open("ct_{:s}_{:s}.txt".format(ct, a), 'w') as f:
            print("ct_{:s}_{:s}".format(ct, a), file=f)
            print("C-terminal ({:s}) {:s}".format(ct, a), file=f)
            print("{:s}".format(ff_file), file=f)
            print("{:s}".format(ct), file=f)
            print(a, file=f)
            print("GLY", file=f)
            print("\nN", file=f)

for nt in n_term:
    for a in aminoacids + alt_aminoacids:
        with open("nt_{:s}_{:s}.txt".format(nt, a), 'w') as f:
            print("nt_{:s}_{:s}".format(nt, a), file=f)
            print("N-terminal ({:s}) {:s}".format(nt, a), file=f)
            print("{:s}".format(ff_file), file=f)
            print("GLY", file=f)
            print(a, file=f)
            print(nt, file=f)
            print("\nN", file=f)
