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
              #'CYX'
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

c_tem = ['H2N', 'FOR', 'ACE', 'PCA']

n_term = ['COH', 'NH2', 'NME']

ff_file = 'amoebabio18.prm'

for a in aminoacids:
    with open("internal_{:s}.txt".format(a), 'w') as f:
        print("internal_{:s}".format(a), file=f)
        print("Internal {:s}".format(a), file=f)
        print("{:s}".format(ff_file), file=f)
        print("GLY", file=f)
        print(a, file=f)
        print("GLY", file=f)
        print("\nN", file=f)
