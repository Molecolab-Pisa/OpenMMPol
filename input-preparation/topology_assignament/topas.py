from pysmiles import read_smiles, write_smiles
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import json

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
