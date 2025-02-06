from pysmiles import read_smiles, write_smiles
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import json
from topas import *
import sys

smiles_aa = {'phe': 'C1=CC=C(C=C1)CC(C(=O)O)N',
             'ala': 'CC(C(=O)O)N',
             'gly': 'C(C(=O)O)N',
             'arg': 'NC(=N)NCCCC(N)C(=O)O',
             'arg+': 'NC(=[NH2])NCCCC(N)C(=O)O',
             'hid': 'C1=C(NC=N1)CC(C(=O)O)N',
             'hip': 'C1=C(NC=([NH]1))CC(C(=O)O)N',
             'hie': 'C1=C([N]C=([NH]1))CC(C(=O)O)N',
             'lys': 'NCCCCC(N)C(=O)O',
             'lys+': '[NH3]CCCCC(N)C(=O)O',
             'asp': 'C([CH](C(=O)O)N)C(=O)[OH]',
             'asp-': 'C([CH](C(=O)O)N)C(=O)[O]',
             'gln': 'NC(=O)CCC(N)C(=O)O',
             'ser': 'OCC(N)C(=O)O',
             'thr': 'CC(O)C(N)C(=O)O',
             'asn': 'C([CH](C(=O)O)N)C(=O)N',
             'glu-': '[O]C(=O)CCC(N)C(=O)O',
             'glu': '[OH]C(=O)CCC(N)C(=O)O',
             'cys': 'SCC(N)C(=O)O',
             'val': 'CC(C)C(N)C(=O)O',
             'ile': 'CC[CH](C)[CH](C(=O)O)N',
             'leu': 'CC(C)C[CH](C(=O)O)N',
             'met': 'CSCCC(C(=O)O)N',
             'tyr': 'Oc1ccc(CC(N)C(=O)O)cc1',
             'trp': 'C(N)(C(=O)O)CC1c2ccccc2NC=1'}

#json_out = {}
#json_out['res_assignament'] = [f.asdict() for f in db]
#for i in json_out['res_assignament']:
#    i['atomtypes'] = []
#
#with open('clean_db.json', 'w') as f:
#    print(json.dumps(json_out, indent=4), file=f)
#exit(1)

 #except:
high_p = ['asp', 'glu', 'lys+', 'hip', 'arg+']
db = []

#Internal residues
for resname in smiles_aa:
    smiles_str, env_atm = aa_to_resid(smiles_aa[resname])
    p = 0
    if resname in high_p:
        p = 1
    db += [Fragment(resname, 
                    smiles_str, 
                    env_atm, 
                    priority=p,
                    default_oln=three2one(resname[:3]),
                    default_resname=resname[:3].upper(),
                    restype='protein_res')]
db += [Fragment('pro', 
                'C1CCN([C](=O))C1C(=O)[N]', 
                [4, 5, 9], 
                priority=0,
                default_oln='P',
                default_resname='PRO',
                restype='protein_res')]

# N-terminal residues
## Formyl H-C(=O)-NH-|
for resname in smiles_aa:
    smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
                                    resid_bb='[CH]([C](=O)[N])(NH([C](=O)[H]))',
                                    resid_ca=0,
                                    resid_env=[3])
    p = 2
    if resname in high_p:
        p = 3
    db += [Fragment(resname+'_nt_nhcoh', 
                    smiles_str, 
                    env_atm, 
                    priority=p,
                    default_oln=three2one(resname[:3]),
                    default_resname=resname[:3].upper(),
                    restype='protein_nhcohterminal')]
## Protonated H3N-|
for resname in smiles_aa:
    smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
                                    resid_bb='[CH]([C](=O)[N])[NH3]',
                                    resid_ca=0,
                                    resid_env=[3])
    p = 2
    if resname in high_p:
        p = 3
    db += [Fragment(resname+'_nt_nh3', 
                    smiles_str, 
                    env_atm, 
                    priority=p,
                    default_oln=three2one(resname[:3]),
                    default_resname=resname[:3].upper(),
                    restype='protein_nh3terminal')]

# C-terminal residues
## NMA H3C-NH-C(=O)-|
for resname in smiles_aa:
    smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
                                    resid_bb='[CH]([C](=O)[NH][CH3])[NH]([C](=O))',
                                    resid_ca=0,
                                    resid_env=[6,7])
    p = 2
    if resname in high_p:
        p = 3
    db += [Fragment(resname+'_ct_conme', 
                    smiles_str, 
                    env_atm, 
                    priority=p,
                    default_oln=three2one(resname[:3]),
                    default_resname=resname[:3].upper(),
                    restype='protein_conmeterminal')]

db += [Fragment('pro_ct_conme', 
                'C1CCN([C](=O))C1C(=O)NC', 
                [4, 5], 
                priority=4,
                default_oln='P',
                default_resname='PRO',
                restype='protein_conmeterminal')]

# Deprotonated O2C-|
m = read_smiles('[CH]([C](=O)[O])[NH]([C](=O))', explicit_hydrogen=True)
tag_mol_graph(m)
for resname in smiles_aa:
    smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
                                    resid_bb='[CH]([C](=O)[O])[NH]([C](=O))',
                                    resid_ca=0,
                                    resid_env=[5,6])
    p = 2
    if resname in high_p:
        p = 3

    db += [Fragment(resname+'_ct_coo', 
                    smiles_str, 
                    env_atm, 
                    priority=p,
                    default_oln=three2one(resname[:3]),
                    default_resname=resname[:3].upper(),
                    restype='protein_cooterminal')]
db += [Fragment('pro', 'C1CCN([C](=O))C1C(=O)[N]', [4, 5, 9], priority=0)]

# Water
db += [Fragment('wat', 
                'O',
                env_atm=[],
                priority=0,
                default_oln='O',
                default_resname='WAT',
                restype='solvent')]
# Ions
db += [Fragment('Zn', 
                '[Zn]',
                env_atm=[],
                priority=0,
                default_oln='Zn',
                default_resname='Zn2',
                restype='solvent')]

json_out = {}
json_out['res_assignament'] = [f.asdict() for f in db]
with open(sys.argv[1], 'w') as f:
    print(json.dumps(json_out, indent=4), file=f)
