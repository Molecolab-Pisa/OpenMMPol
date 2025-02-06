from pysmiles import read_smiles, write_smiles
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import json

import sys


class Fragment():
    def __init__(self,
                 name,
                 smiles_str,
                 env_atm=[],
                 atomnames=[],
                 atomtypes=[],
                 default_resname='UNK',
                 restype='unknown',
                 default_oln='-',
                 priority=0):
        self.name = name
        self.smiles_str = smiles_str
        self.env_atm = env_atm
        self.natoms = len(read_smiles(self.smiles_str, explicit_hydrogen=True).nodes)
        if default_resname is not None:
            self.resname = default_resname
        else:
            self.resname = name
        self.oln = default_oln
        self.priority = priority
        self.atomtypes = atomtypes
        self.atomnames = atomnames
        self.restype = restype

    @property
    def assigned_atm(self):
        return [i for i in range(self.natoms) if i not in self.env_atm]

    def set_atomtypes(self, atomtypes):
        if self.natoms == len(atomtypes):
            self.atomtypes = atomtypes.copy()
        else:
            raise ValueError("Length of the atomtypes array should be {:d}".format(natoms))

    def asdict(self):
        d = {}
        d['name'] = self.name
        d['smiles_str'] = self.smiles_str
        d['env_atm'] = list(self.env_atm)
        d['natoms'] = self.natoms
        d['resname'] = self.resname
        d['priority'] = self.priority
        d['atomtypes'] = [int(i) for i in self.atomtypes]
        d['atomnames'] = list(self.atomnames)
        d['restype'] = self.restype
        d['oln'] = self.oln
        return d

def dict_to_frag(d):
    return Fragment(d['name'],
                    d['smiles_str'],
                    env_atm=d['env_atm'],
                    default_resname=d['resname'],
                    priority=d['priority'],
                    atomtypes=d['atomtypes'],
                    atomnames=d['atomnames'],
                    restype=d['restype'],
                    default_oln=d['oln'])
def name_to_c(n):
    if n == 'C':
        return '#000000'
    elif n == 'H':
        return '#dddddd'
    elif n == 'O':
        return '#ff0000'
    elif n == 'N':
        return '#0000ff'
    else:
        return '#00ff00'

def tag_mol_graph(g):
    nx.set_node_attributes(g, -1, 'id')

    for i in range(len(g.nodes)):
        g.nodes[i]['id'] = i

def get_subgraph_matches(big_g, sub_g, exclude_list=None):
    def compare_atoms(a, b):
        if exclude_list is not None:
            if a['id'] in exclude_list:
                return False

        if a['element'] == b['element']:
            return True
        else:
            return False

    def compare_bonds(a, b):
        # All bonds are equal
        return True

    sub_iter = nx.isomorphism.GraphMatcher(big_g,
                                           sub_g,
                                           node_match=compare_atoms,
                                           edge_match=compare_bonds).subgraph_isomorphisms_iter()
    return remove_chemical_equivalence([a for a in sub_iter], big_g)

def remove_chemical_equivalence(subs, big_g):
    # Two subgraphs are chemically equivalent if
    #   a) all the nodes in the graph are the same
    cheq = {}
    for ii, s in enumerate(subs):
        done = False
        for i in cheq:
            if ii in cheq[i]:
                done = True
        if done:
            continue

        cheq[ii] = []
        s_nodes = set([i for i in s])
        for j, r in enumerate(subs):
            r_nodes = set([i for i in r])
            if len(s_nodes.intersection(r_nodes)) == len(s_nodes) and j != ii:
                cheq[ii] += [j]

    return [subs[i] for i in cheq]


def get_subgraph_match(big_g, sub_g, exclude_list=None):
    ms = get_subgraph_matches(big_g, sub_g, exclude_list)
    if len(ms) == 1:
        return ms[0]
    elif len(ms) == 0:
        return None
    else:
        raise ValueError("More than one match was found!")

def print_mol(m):
    nx.draw(m,
            node_color = [name_to_c(a[1]) for a in m.nodes(data='element')],
            labels = {n: m.nodes[n]['id'] for n in m}, font_color='#00ff00'
            )
    plt.show()

def aa_to_resid(in_aa_string, 
                resid_bb='[CH]([C](=O)[N])[NH]([C](=O))',
                resid_ca=0,
                resid_env=[3, 5, 6]):
    neutral_bb = '[CH](C(=O)[OH])[NH2]'
    m_neut_bb = read_smiles(neutral_bb,
                            explicit_hydrogen=True)
    neutral_ca = 0
    #tag_mol_graph(m_neut_bb)
    #print_mol(m_neut_bb)
    m_resid_bb = read_smiles(resid_bb,
                             explicit_hydrogen=True)
    
    m_aa = read_smiles(in_aa_string,
                       explicit_hydrogen=True)
    tag_mol_graph(m_aa)

    # Find the atom linked to CA and remove the whole backbone
    aa_ca = None
    try:
        aa2neutral = get_subgraph_match(m_aa, m_neut_bb)
    except ValueError:
        # This should be glycine!
        assert get_subgraph_match(m_aa, read_smiles('[CH2](C(=O)[OH])[NH2]')) is not None
        aa2neutral = get_subgraph_matches(m_aa, m_neut_bb)[0]
    if aa2neutral is None:
        print("Prolineeeee!!!")
        pro_m_neut_bb = read_smiles('[CH](C(=O)[OH])[NH]', explicit_hydrogen=True)
        aa2neutral = get_subgraph_match(m_aa, pro_m_neut_bb)
        print_mol(m_aa)
        tag_mol_graph(pro_m_neut_bb)
        print_mol(pro_m_neut_bb)


    toremove = []
    for i_aa in aa2neutral:
        i_neutral = aa2neutral[i_aa]
        toremove += [i_aa]

        if i_neutral == neutral_ca:
            for a, b in m_aa.edges(i_aa):
                try:
                    assert a == i_aa
                    c = b
                except AssertionError:
                    c = a

                if aa2neutral.get(c, None) is None:
                    if aa_ca is not None:
                        raise ValueError("Problem in locating sidechain")
                    aa_ca = c
            #print("Side chain joint point is ", c)
    m_aa.remove_nodes_from(toremove)
    m_resid_aa = nx.disjoint_union(m_resid_bb, m_aa)
    tag_mol_graph(m_resid_aa)

    resid2bb = get_subgraph_match(m_resid_aa, m_resid_bb)

    # Exclude the backbone from the next search...
    resid2sidechain = get_subgraph_match(m_resid_aa, m_aa, exclude_list=[i for i in resid2bb])
    residaa_s = None
    residaa_ca = None
    for i in resid2sidechain:
        if resid2sidechain[i] == aa_ca:
            if residaa_s is not None:
                raise ValueError("Problem in locating sidechain in final mol")
            residaa_s = i
    for i in resid2bb:
        if resid2bb[i] == resid_ca:
            if residaa_ca is not None:
                raise ValueError("Problem in locating CA in final mol")
            residaa_ca = i
    m_resid_aa.add_edge(residaa_s, residaa_ca)
    #print_mol(m_resid_aa)
    smiles_str = write_smiles(m_resid_aa)
    reload_map = get_subgraph_match(m_resid_aa, read_smiles(smiles_str))

    return smiles_str, [reload_map[i] for i in resid_env]

def selection_to_graph(sele):
    # Fix missing informations
    u = sele.universe
    if not hasattr(u.atoms, 'elements'):
        u.add_TopologyAttr('element',
                           [guess_atom_element(at.name) for at in u.atoms])

    if not hasattr(u.atoms, 'bonds') or sum([len(at.bonds) for at in sele]) == 0:
        # Create a VdW table that only relies on element information
        if hasattr(u.atoms, 'types'):
            vdw_table = {}
            for at in u.atoms:
                if at.type not in vdw_table:
                    vdw_table[at.type] = mda.topology.tables.vdwradii[at.element.upper()]
        else:
            vdw_table = None

        if not hasattr(u.atoms, 'bonds'):
            u.atoms.guess_bonds(vdwradii=vdw_table)
        else:
            sele.atoms.guess_bonds(vdwradii=vdw_table)

    G = nx.Graph(attr='element')
    natoms = 0
    for i, at in enumerate(sele):
        e = at.element
        if len(e) > 1:
            e = e[0] + e[1:].lower()
        G.add_node(at.index, element=e)
        natoms += 1
    nbond = 0
    for at in sele:
        for b in at.bonds:
            G.add_edge(at.index, b.partner(at).index)
            nbond += 1
    return G

def three2one(strin):
    three_one = {'phe': 'F',
                 'ala': 'A',
                 'gly': 'G',
                 'arg': 'R',
                 'hi':  'H',
                 'lys': 'K',
                 'asp': 'D',
                 'gln': 'Q',
                 'ser': 'S',
                 'thr': 'T',
                 'asn': 'N',
                 'glu': 'E',
                 'cys': 'C',
                 'pro': 'P',
                 'val': 'V',
                 'ile': 'I',
                 'leu': 'L',
                 'met': 'M',
                 'tyr': 'Y',
                 'trp': 'W',
                 'wat': 'O',
                 'Mg': 'Mg',
                 'Zn': 'Zn'}
    for t in three_one:
        if strin.startswith(t):
            return three_one[t]
    return '-'

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

print("Loading db from file...")
with open('clean_db.json', 'r') as f:
    json_in = json.load(f)
db = [dict_to_frag(f) for f in json_in['res_assignament']]

json_out = {}
json_out['res_assignament'] = [f.asdict() for f in db]
for i in json_out['res_assignament']:
    i['atomtypes'] = []

with open('clean_db.json', 'w') as f:
    print(json.dumps(json_out, indent=4), file=f)
exit(1)

## #except:
## high_p = ['asp', 'glu', 'lys+', 'hip', 'arg+']
## #    db = []
## #
## #    #Internal residues
## #    for resname in smiles_aa:
## #        smiles_str, env_atm = aa_to_resid(smiles_aa[resname])
## #        p = 0
## #        if resname in high_p:
## #            p = 1
## #        db += [Fragment(resname, 
## #                        smiles_str, 
## #                        env_atm, 
## #                        priority=p,
## #                        default_oln=three2one(resname[:3]),
## #                        default_resname=resname[:3].upper(),
## #                        restype='protein_res')]
## #    db += [Fragment('pro', 
## #                    'C1CCN([C](=O))C1C(=O)[N]', 
## #                    [4, 5, 9], 
## #                    priority=0,
## #                    default_oln='P',
## #                    default_resname='PRO',
## #                    restype='protein_res')]
## #
##     #N-terminal residues (NH3+)
## m = read_smiles('[CH]([C](=O)[N])(NH(C(=O)H))', explicit_hydrogen=True)
## tag_mol_graph(m)
## print_mol(m)
## for resname in smiles_aa:
##     smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
##                                     resid_bb='[CH]([C](=O)[N])(NH([C](=O)[H]))',
##                                     resid_ca=0,
##                                     resid_env=[3])
##     p = 2
##     if resname in high_p:
##         p = 3
##     db += [Fragment(resname+'_nt_nhcoh', 
##                     smiles_str, 
##                     env_atm, 
##                     priority=p,
##                     default_oln=three2one(resname[:3]),
##                     default_resname=resname[:3].upper(),
##                     restype='protein_nhcohterminal')]
##     #db += [Fragment('pro', 'C1CCN([C](=O))C1C(=O)[N]', [4, 5, 9], priority=0)]
## #    #N-terminal residues (NH3+)
## #    for resname in smiles_aa:
## #        smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
## #                                        resid_bb='[CH]([C](=O)[N])[NH3]',
## #                                        resid_ca=0,
## #                                        resid_env=[3])
## #        p = 2
## #        if resname in high_p:
## #            p = 3
## #        db += [Fragment(resname+'_nt_nh3', 
## #                        smiles_str, 
## #                        env_atm, 
## #                        priority=p,
## #                        default_oln=three2one(resname[:3]),
## #                        default_resname=resname[:3].upper(),
## #                        restype='protein_nh3terminal')]
## #    #db += [Fragment('pro', 'C1CCN([C](=O))C1C(=O)[N]', [4, 5, 9], priority=0)]
## #
## #C_terminal residues (COO)
## for resname in smiles_aa:
##     smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
##                                     resid_bb='[CH]([C](=O)[NH][CH3])[NH]([C](=O))',
##                                     resid_ca=0,
##                                     resid_env=[6,7])
##     #ma = read_smiles(smiles_str)
##     #tag_mol_graph(ma)
##     #print_mol(ma)
##     p = 5
##     if resname in high_p:
##         p = 4
##     db += [Fragment(resname+'_ct_conme', 
##                     smiles_str, 
##                     env_atm, 
##                     priority=p,
##                     default_oln=three2one(resname[:3]),
##                     default_resname=resname[:3].upper(),
##                     restype='protein_conmeterminal')]
## 
## ma = read_smiles('C1CCN([C](=O))C1C(=O)NC')
## tag_mol_graph(ma)
## print_mol(ma)
## db += [Fragment('pro_ct_conme', 
##                 'C1CCN([C](=O))C1C(=O)NC', 
##                 [4, 5], 
##                 priority=4,
##                 default_oln='P',
##                 default_resname='PRO',
##                 restype='protein_conmeterminal')]
##     #C_terminal residues (CONHMe)
##     #m = read_smiles('[CH]([C](=O)[O])[NH]([C](=O))', explicit_hydrogen=True)
##     #tag_mol_graph(m)
##     #for resname in smiles_aa:
##     #    smiles_str, env_atm = aa_to_resid(smiles_aa[resname],
##     #                                    resid_bb='[CH]([C](=O)[O])[NH]([C](=O))',
##     #                                    resid_ca=0,
##     #                                    resid_env=[5,6])
##     #    p = 2
##     #    if resname in high_p:
##     #        p = 3
##     #    db += [Fragment(resname+'_ct_coo', 
##     #                    smiles_str, 
##     #                    env_atm, 
##     #                    priority=p,
##     #                    default_oln=three2one(resname[:3]),
##     #                    default_resname=resname[:3].upper(),
##     #                    restype='protein_cooterminal')]
##     #db += [Fragment('pro', 'C1CCN([C](=O))C1C(=O)[N]', [4, 5, 9], priority=0)]
## #
## #    #Water
## #    db += [Fragment('wat', 
## #                    'O',
## #                    env_atm=[],
## #                    priority=0,
## #                    default_oln='O',
## #                    default_resname='WAT',
## #                    restype='solvent')]
## db += [Fragment('Zn', 
##                 '[Zn]',
##                 env_atm=[],
##                 priority=0,
##                 default_oln='Zn',
##                 default_resname='Zn2',
##                 restype='solvent')]
## 
## json_out = {}
## json_out['res_assignament'] = [f.asdict() for f in db]
## with open('clean_db.json', 'w') as f:
##     print(json.dumps(json_out, indent=4), file=f)

uni = mda.Universe(sys.argv[1])
prot = selection_to_graph(uni.atoms)
tag_mol_graph(prot)
nat = len(prot.nodes)
done_flags = []
resnames = ['UNK'] * nat
resid = np.zeros(nat)
atmid = np.zeros(nat)
for ifrag, frag in sorted(enumerate(db), key=lambda a: a[1].priority, reverse=True):
    resname = frag.name
    res = read_smiles(frag.smiles_str, explicit_hydrogen=True)
    aa_matches = get_subgraph_matches(prot, res)
    print("Found {:d} {:s} resids".format(len(aa_matches), resname))
    #print(frag.env_atm)
    #if len(aa_matches) > 0:
    #    tag_mol_graph(res)
    #    print_mol(res)
    #attype = np.zeros(frag.natoms, dtype=np.int64)
    #if len(frag.atomtypes) == 0 or 0 in [frag.atomtypes[i] for i in frag.assigned_atm]:
    #    missing_prm = True
    #else:
    #    missing_prm = False

    for m in aa_matches:
        candidate_resid = max(resid) + 1
        for i in m:
            if m[i] not in frag.env_atm:
                if i in done_flags:
                    pass
                else:
                    done_flags += [i]
                    resnames[i] = resname
                    resid[i] = candidate_resid
                    atmid[i] = m[i]
                    # if missing_prm:
                    #     attype[m[i]] = int(uni.atoms[i].type)
                    # else:
                    #     try:
                    #         assert frag.atomtypes[m[i]] == int(uni.atoms[i].type)
                    #     except:
                    #         print("Atom {}: expected {} found {}".format(i, frag.atomtypes[m[i]], int(uni.atoms[i].type)))
        #if missing_prm and not 0 in attype[frag.assigned_atm]:
        #    print("Added parameters for {}".format(frag.name))
        #    db[ifrag].set_atomtypes(attype)
    print(len(done_flags))



lastresid = -1
for i in range(nat):
    if resid[i] != lastresid:
        print(three2one(resnames[i]), end='')
        lastresid = resid[i]

print()
print("Assigned {:d} / {:d}".format(len(done_flags), nat))
