from pysmiles import read_smiles, write_smiles
from topas import *
import sys

db = AssignamentDB()
with open('internal.txt', 'r') as f:
    for l in f:
        l = l.split('#')[0]
        if not l.strip() :
            continue
        tok = l.split()
        name = tok[0]
        smiles_str = tok[1]
        if tok[2] == '-':
            env_atm = []
        else:
            env_atm = [int(i) for i in tok[2].split(',')]
        threelc = tok[3]
        restype = tok[4]
        
        f = Fragment(name,
                        smiles_str,
                        env_atm,
                        default_resname=threelc,
                        restype=restype)
        db.add_fragment(f)


with open('c_terminal_backbones.txt', 'r') as f:
    for l in f:
        l = l.split('#')[0]
        if not l.strip() :
            continue
        tok = l.split()
        name = tok[0]
        name_pre = name.split('_')[0]
        bb_smiles_str = tok[1]
        if tok[2] == '-':
            bb_env_atm = []
        else:
            bb_env_atm = [int(i) for i in tok[2].split(',')]
        threelc = tok[3]
        restype = tok[4]
        with open(name_pre+'_c_terminal_res.txt', 'w') as fout:
            for f in db:
                if f.is_aminoacid_residue():
                    smiles_str, env_atm = f.replace_aminoacid_backbone(bb_smiles_str, bb_env_atm)
                    print('{:<20s}{:<50s}{:>20}{:>5s}\t{:s}'.format(name_pre+'_'+f.name,
                                                                smiles_str, 
                                                                ','.join([str(a) for a in env_atm]),
                                                                threelc,
                                                                restype), file=fout)

with open('n_terminal_backbones.txt', 'r') as f:
    for l in f:
        l = l.split('#')[0]
        if not l.strip() :
            continue
        tok = l.split()
        name = tok[0]
        name_pre = name.split('_')[0]
        bb_smiles_str = tok[1]
        if tok[2] == '-':
            bb_env_atm = []
        else:
            bb_env_atm = [int(i) for i in tok[2].split(',')]
        threelc = tok[3]
        restype = tok[4]
        with open(name_pre+'_n_terminal_res.txt', 'w') as fout:
            for f in db:
                if f.is_aminoacid_residue():
                    smiles_str, env_atm = f.replace_aminoacid_backbone(bb_smiles_str, bb_env_atm)
                    print('{:<20s}{:<50s}{:>20}{:>5s}\t{:s}'.format(name_pre+'_'+f.name,
                                                                smiles_str, 
                                                                ','.join([str(a) for a in env_atm]),
                                                                threelc,
                                                                restype), file=fout)

