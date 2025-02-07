from pysmiles import read_smiles, write_smiles
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import json

import sys
from topas import *

logger = logging.getLogger('pysmiles')
logger.setLevel(level=logging.ERROR)

def topology_assign(input_structure, json_database):
    print("Loading db from file...")
    with open(json_database, 'r') as f:
        json_in = json.load(f)
    db = [dict_to_frag(f) for f in json_in['res_assignament']]

    print("Loading input structure and converting to graph...")
    uni = mda.Universe(input_structure)
    prot = selection_to_graph(uni.atoms)
    tag_mol_graph(prot)
    nat = len(prot.nodes)

    print("Sub graph recognition...")
    done_flags = []
    resnames = ['UNK'] * nat
    resid = np.zeros(nat)
    atmid = np.zeros(nat)
    for ifrag, frag in sorted(enumerate(db), key=lambda a: a[1].priority, reverse=True):
        resname = frag.name
        print("Searching for residue {:s}".format(resname))
        res = read_smiles(frag.smiles_str, explicit_hydrogen=True, strict=False)
        aa_matches = get_subgraph_matches(prot, res)
        if len(aa_matches) > 0:
            print("Found {:d} {:s} resids".format(len(aa_matches), resname))

        for m in aa_matches:
            candidate_resid = max(resid) + 1

            assignable_atoms = [i for i in m if m[i] not in frag.env_atm]
            print(m)
            print([a in done_flags for a in assignable_atoms])
            #if any([a in done_flags for a in assignable_atoms]):
            #    print("The structure is already assigned!")
            #    continue
            #assigning = None
            for i in m:
                if m[i] not in frag.env_atm:
                    if i in done_flags:
                        #if assigning == True:
                        #    print(m[i])
                        #    tag_mol_graph(res)
                        #    print_mol(res)
                        #    raise ValueError("Trying to partially assign a residue")
                        #else:
                        #    assigning = False
                        pass
                    else:
                        #if assigning == False:
                        #    print(m[i])
                        #    tag_mol_graph(res)
                        #    print_mol(res)
                        #    raise ValueError("Trying to partially assign a residue")
                        #else:
                        #    assigning = True

                        done_flags += [i]
                        resnames[i] = resname
                        resid[i] = candidate_resid
                        atmid[i] = m[i]

    lastresid = -1
    for i in range(nat):
        if resid[i] != lastresid:
            print(three2one(resnames[i]), end='')
            lastresid = resid[i]

    print()
    print("Assigned {:d} / {:d}".format(len(done_flags), nat))

topology_assign(sys.argv[1], sys.argv[2])
