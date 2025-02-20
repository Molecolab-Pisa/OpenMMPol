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
        
        if tok[3] == '-':
            bridge_atm = []
        else:
            bridge_atm = [int(i) for i in tok[3].split(',')]

        threelc = tok[4]
        restype = tok[5]
        
        f = Fragment(name,
                        smiles_str,
                        env_atm,
                        bridge_atm,
                        default_resname=threelc,
                        restype=restype)
        db.add_fragment(f)


for bbfile, out_base in [['c_terminal_backbones.txt', '_c_terminal_res.txt'],
                         ['n_terminal_backbones.txt', '_n_terminal_res.txt'],
                         ['double_cap_backbones.txt', '_double_cap_res.txt']]:
    with open(bbfile, 'r') as f:
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
            if tok[3] == '-':
                bb_bridge_atm = []
            else:
                bb_bridge_atm = [int(i) for i in tok[3].split(',')]
            threelc = tok[4]
            restype = tok[5]
            print(l)
            with open(name_pre+out_base, 'w') as fout:
                for f in db:
                    if f.is_aminoacid_residue():
                        smiles_str, env_atm, bridge_atm = f.replace_aminoacid_backbone(bb_smiles_str, bb_env_atm, bb_bridge_atm)
                        
                        env_str = ','.join([str(a) for a in env_atm])
                        if not env_str:
                            env_str = '-'
                        
                        bridge_str = ','.join([str(a) for a in bridge_atm])
                        if not bridge_str:
                            bridge_str = '-'
                        print('{:<20s}{:<50s}{:>20}{:>20s}{:>5s}\t{:s}'.format(name_pre+'_'+f.name,
                                                                    smiles_str, 
                                                                    env_str,
                                                                    bridge_str,
                                                                    f.resname,
                                                                    restype), file=fout)

