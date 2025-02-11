from pysmiles import read_smiles, write_smiles
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import json
import logging

logger = logging.getLogger('pysmiles')
logger.setLevel(level=logging.ERROR)

class PrmAssignament():
    def __init__(self, uni):
        self.u = uni
        # For each atom the residue it belongs to
        self.assigned_resid = np.zeros(self.natoms, dtype=np.int64) - 1
        # For each residue, the fragment id in the database
        self.assigned_fragid = np.zeros(self.natoms, dtype=np.int64) - 1
        # For each atom, the atomid within the assigned residue in the db
        self.assigned_atomid = np.zeros(self.natoms, dtype=np.int64) - 1
        # Molecule graph
        self.mol_graph = selection_to_graph(self.u.atoms)
        # Tag the mol graph
        tag_mol_graph(self.mol_graph)
        # Check for correctness
        assert len(self.mol_graph.nodes) == self.natoms
        # If needed, generate sub-graphs for each molecule
        self.split_into_molecules()
        # Database used for the assignament
        self.db = None

    @property
    def natoms(self):
        """Return the number of atoms in the system"""
        return len(self.u.atoms)

    def is_assigned(self, iatm: int):
        """Check if an atom with index iatm is assigned or not"""
        return self.assigned_fragid[iatm] >= 0

    @property
    def tot_assigned(self):
        """Count the total number of assigned atoms"""
        return np.count_nonzero(self.assigned_fragid >= 0)

    def get_next_free_resid(self):
        """Get the lowest residue-id that is not used"""
        return max(self.assigned_resid) + 1

    def split_into_molecules(self):
        """Check if more than one molecule is present in the
        system and generate the subgraphs for each molecule.
        Used to speedup the assignament procedure"""
        cc_gen = nx.connected_components(self.mol_graph)
        self.mol_sub_graphs = []
        for comp in cc_gen:
            self.mol_sub_graphs += [nx.subgraph(self.mol_graph, comp)]

    @property
    def nmol(self):
        """Number of molecules in the system"""
        if self.mol_sub_graphs is None:
            self.split_into_molecules()
        return len(self.mol_sub_graphs)

    def set_db(self, db):
        self.db = db

    # def get_ol_sequence(self):
    #     try:
    #         assert self.tot_assigned == self.natoms
    #     except AssertionError:
    #         print("The structure should be completely assigned before asking for sequence")

    def get_assigned_atomtypes(self):
        at = []
        for i in range(self.natoms):
            if self.is_assigned(i):
                try:
                    frag = self.db[self.assigned_fragid[i]]
                except KeyError:
                    print("Can't find fragment with fragment id {:d}".format(self.assigned_fragid[i]))
                    return None
                try:
                    at += [frag.atomtypes[self.assigned_atomid[i]]]
                except KeyError:
                    print("Atom {:d} assigned as atom {:d} of residue {:d} (...) hasn't an atomtype".format(i, 
                                                                                                            self.assigned_atomid[i],
                                                                                                            frag.fragment_id))
                    return None
            else:
                return None
        return at

    def topology_assign(self):
        """Perform the assignament using the current db"""
        try:
            assert self.db is not None
        except AssertionError:
            raise ValueError("A database should be set before starting the assignament")
        
        
        if self.mol_sub_graphs is None:
            mgs = [self.mol_graph]
        else:
            mgs = self.mol_sub_graphs

        for mg in mgs:
            print("Assigning molecule with {:d} nodes and {:d} arcs".format(len(mg.nodes), len(mg.edges)))

            for frag in self.db:
                
                # Skip if the number of atoms in fragment is larger than the
                # number of atoms in current molecule
                if frag.natoms > len(mg.nodes):
                    continue

                print("Searching for residue {:s} (priority {:d})".format(frag.name,
                                                                        frag.priority))
                aa_matches = get_subgraph_matches(mg, frag.frag_graph)
                if len(aa_matches) > 0:
                    print("Found {:d} {:s} resids".format(len(aa_matches), frag.name))

                for m in aa_matches:
                    candidate_resid = self.get_next_free_resid()

                    assignable_atoms = [i for i in m if m[i] not in frag.env_atm]
                    print(assignable_atoms)

                    if any([self.is_assigned(a) for a in assignable_atoms]):
                        print("The structure is already assigned!")
                        res_to_unassign = []
                        at_to_unassign = []
                        new_assigned_at = []
                        for at in assignable_atoms:
                            if self.is_assigned(at):
                                print("Atom {:03d} is in resid {:03d}".format(at, self.assigned_resid[at]))
                                res_to_unassign += [self.assigned_resid[at]]
                                at_to_unassign += [at]
                            else:
                                print("Atom {:03d} is not in a residue!")
                                new_assigned_at += [at]
                        #print(set(res_to_unassign))
                        #print(len(at_to_unassign))
                        #print(new_assigned_at)
                        if len(set(res_to_unassign)) > 1:
                            raise NotImplementedError("Multiple residues have to be unassigned!")
                        if len(new_assigned_at) == 0:
                            print("No improvement here.")
                            continue
                        else:
                            #print(np.count_nonzero(resid == res_to_unassign[0]))
                            #print(len(at_to_unassign))
                            #print(len(new_assigned_at))
                            if np.count_nonzero(resid == res_to_unassign[0]) == len(at_to_unassign):
                                print("Ok reassigning is needed here!")
                                exit()
                            else:
                                raise NotImplementedError("Partial reassignement is not expected!")
                        continue

                    print("Assigning {:d} atoms to residue {:d}".format(len(m), candidate_resid))
                    for i in m:
                        if m[i] not in frag.env_atm:
                            if self.is_assigned(i):
                                #if assigning == True:
                                #    print(m[i])
                                #    tag_mol_graph(res)
                                #    print_mol(res)
                                #    raise ValueError("Trying to partially assign a residue")
                                #else:
                                #    assigning = False
                                raise ValueError("Attempting unhandled ressignament")
                            else:
                                #if assigning == False:
                                #    print(m[i])
                                #    tag_mol_graph(res)
                                #    print_mol(res)
                                #    raise ValueError("Trying to partially assign a residue")
                                #else:
                                #    assigning = True

                                self.assigned_resid[i] = candidate_resid
                                self.assigned_fragid[i] = frag.fragment_id
                                self.assigned_atomid[i] = m[i]

    def learn_from(self, use_name=True, use_type=True):
        if use_name:
            raise NotImplementedError("Learning atomnames is still not implemented")
        if use_type:
            self.learn_type()

    def learn_type(self):
        if not hasattr(self.u.atoms[0], 'type'):
            raise ValueError("Cannot learn, sine your input file has no atom types")

        for i in range(self.natoms):
            if self.is_assigned(i):
                print("Correct atom type for {:03d} is {}".format(i, self.u.atoms[i].type))
                # Find the correct fragment
                frag = self.db[self.assigned_fragid[i]]
                print("Atom is assigned to {:s}-{:d}".format(frag.name, self.assigned_atomid[i]))
                #print(frag.atomtypes)
                if frag.has_atomtypes(self.assigned_atomid[i]):
                    try:
                        assert self.u.atoms[i].type == frag.atomtypes[self.assigned_atomid[i]]
                        print("Already in DB with consistent atomtype")
                    except AssertionError:
                        raise ValueError("""Inconsistent information between learning input 
                        and database for atom {:d} (atom {:d} of residue {:s})""".format(i, self.assigned_atomid[i], frag.name))
                else:
                    frag.set_atomtype(self.assigned_atomid[i], self.u.atoms[i].type)
                    print("Inserting in db")

            else:
                print("Atom {:03d} is unassigned".format(i))

    def save_tinker_xyz(self, fname, header_str=None):
        if header_str is None:
            header = 'Generated by TopologyAssignament by Mattia Bondanza'
        else:
            header = header_str.replace('\n', ' ')
            
        assigned_atomtypes = self.get_assigned_atomtypes()
        if assigned_atomtypes is None:
            print("Cannot continue without atom types")
            return

        with open(fname, 'w') as f:
            print("{:d} {:s}".format(self.natoms, header), file=f)
            for i in range(self.natoms):
                at = self.u.atoms[i]
                print("{:>5d} {:>5s} {:10.6f} {:10.6f} {:10.6f} {:>5d} ".format(at.id, 
                                                                                at.element,
                                                                                at.position[0],
                                                                                at.position[1],
                                                                                at.position[2],
                                                                                int(assigned_atomtypes[i])), 
                      end='', file=f)
                for b in at.bonds:
                    if b.atoms[0] == at:
                        at2 = b.atoms[1]
                    else:
                        at2 = b.atoms[0]
                    print('{:>5d} '.format(at2.id), end='', file=f)
                print(file=f)

    def assignament_log(self):
        s = ''
        for i, mol in enumerate(self.mol_sub_graphs):
            s += 'Molecule {:d} ({:d} atoms): '.format(i, len(mol.nodes))
            resids = []
            s1 = ''
            nuna = 0
            for atid in mol.nodes:
                if self.assigned_resid[atid] < 0:
                    if not s1.endswith('UNASSIGNED'):
                        s1 += '-UNASSIGNED'
                    nuna += 1
                else:
                    if self.assigned_resid[atid] in resids:
                        if self.assigned_resid[atid] == resids[-1]:
                            pass
                        else:
                            s1 += '?(uncontiguous residues found)'
                    else:
                        if s1.endswith('UNASSIGNED'):
                            s1 = s1.replace('UNASSIGNED', '?({:d} unassigned atoms)'.format(nuna))
                            nuna = 0
                        if len(s1) > 0:
                            s1 += '-'
                        resids += [self.assigned_resid[atid]]
                        s1 += self.db[self.assigned_fragid[atid]].resname
            if s1.endswith('UNASSIGNED'):
                s1 = s1.replace('UNASSIGNED', '?({:d} unassigned atoms)'.format(nuna))
            s += s1 + '\n'
        return s


class AssignamentDB():
    def __init__(self, dbfile=None):
        self.dbfile = dbfile

        if self.dbfile is not None:
            with open(self.dbfile, 'r') as f:
                json_in = json.load(f)
            db = [Fragment.from_dict(f) for f in json_in['res_assignament']]
            self.db = sorted(db, key=lambda frag: frag.priority, reverse=True)
            self.id_to_index = {}
            for ifrag, frag in enumerate(self.db):
                frag.fragment_id = ifrag
                self.id_to_index[ifrag] = ifrag
        else:
            self.db = []
            self.id_to_index = {}

    def __iter__(self):
        return iter(sorted(self.db, key=lambda frag: frag.priority, reverse=True))

    def __getitem__(self, k):
        if type(k) == int or np.can_cast(type(k), int):
            try:
                return self.db[self.id_to_index[k]]
            except KeyError:
                raise StopIteration
        else:
            raise KeyError

    def save_as_json(self, fout):
        json_out = {}
        json_out['res_assignament'] = [f.asdict() for f in self.db]
        with open(fout, 'w') as f:
            print(json.dumps(json_out, indent=4), file=f)

    def get_min_free_id(self):
        for i, k in enumerate(sorted(self.id_to_index.keys())):
            if i != k:
                return i
        return len(self.id_to_index)

    def add_fragment(self, frag):
        # TODO add check on new added fragments
        self.db += [frag]
        fragment_id = self.get_min_free_id()
        self.id_to_index[fragment_id] = len(self.db) - 1
        self.db[-1].fragment_id = fragment_id

    def remove_fragment(self, frag):
        fragment_id = frag.fragment_id
        try:
            index = self.id_to_index[fragment_id]
        except IndexError:
            print("Trying to remove an unexistent fragment")
            return
        self.db.pop(index)
        del self.id_to_index[fragment_id]
        for fid in self.id_to_index:
            if self.id_to_index[fid] > index:
                self.id_to_index[fid] -= 1

class Fragment():
    def __init__(self,
                 name,
                 smiles_str,
                 env_atm=[],
                 atomnames={},
                 atomtypes=None,
                 default_resname='UNK',
                 restype='unknown',
                 default_oln='-',
                 priority=None):

        self.name = name
        self.smiles_str = smiles_str
        if '/' in self.smiles_str or '\\' in self.smiles_str:
            print("Stereochemical information will be discarded!")
            self.smiles_str = self.smiles_str.replace('/', '-')
            self.smiles_str = self.smiles_str.replace('\\', '-')
        self.frag_graph = read_smiles(self.smiles_str,
                                      explicit_hydrogen=True,
                                      strict=False)
        tag_mol_graph(self.frag_graph)
        self.fragment_id = None
        self.env_atm = env_atm
        self.natoms = len(self.frag_graph.nodes)

        if default_resname is not None:
            self.resname = default_resname
        else:
            self.resname = name
        self.oln = default_oln
        if priority is None:
            self.priority = (self.natoms - len(env_atm)) * 10 + len(env_atm)
        else:
            self.priority = priority
        if atomtypes is None:
            self.atomtypes = {}
        else:
            self.atomtypes = atomtypes
        self.atomnames = atomnames
        self.restype = restype

    @property
    def assignable_atm(self):
        return [i for i in range(self.natoms) if i not in self.env_atm]

    def has_atomtypes(self, iatm=None):
        if iatm is None:
            return len(self.atomtypes) == self.assignable_atm
        else:
            if iatm in self.atomtypes:
                return self.atomtypes[iatm] is not None
            else:
                return False

    def set_atomtype(self, iatm, atomtype):
        if iatm in self.atomtypes and self.atomtypes[iatm] is not None:
            raise ValueError("Cannot reassign atomtypes")
        else:
            self.atomtypes[iatm] = atomtype
    
    def reset_atomtype(self, iatm):
        self.atomtypes[iatm] = None

    def asdict(self):
        d = {}
        d['name'] = self.name
        d['smiles_str'] = self.smiles_str
        d['env_atm'] = list(self.env_atm)
        d['natoms'] = self.natoms
        d['resname'] = self.resname
        d['priority'] = self.priority
        d['atomtypes'] = {int(a): self.atomtypes[a] for a in self.atomtypes}
        d['atomnames'] = {int(a): self.atomnames[a] for a in self.atomnames}
        d['restype'] = self.restype
        d['oln'] = self.oln
        return d

    @classmethod
    def from_dict(cls, d):
        dic_atomtypes = d['atomtypes']
        dic_atomnames = d['atomnames']

        atomtypes = {int(a): dic_atomtypes[a] for a in dic_atomtypes}
        atomnames = {int(a): dic_atomnames[a] for a in dic_atomnames}
        return cls(d['name'],
                   d['smiles_str'],
                   env_atm=d['env_atm'],
                   default_resname=d['resname'],
                   priority=d['priority'],
                   atomtypes=atomtypes,
                   atomnames=atomnames,
                   restype=d['restype'],
                   default_oln=d['oln'])

    def draw(self):
        def name_to_color(n):
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
        nx.draw(self.frag_graph,
                pos = nx.spring_layout(self.frag_graph, iterations=100*self.natoms),
                node_color = [name_to_color(a[1]) for a in self.frag_graph.nodes(data='element')],
                labels = {n: self.frag_graph.nodes[n]['id'] for n in self.frag_graph}, font_color='#00ff00'
                )
        plt.show()



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

def aa_to_resid(in_aa_string,
                resid_bb='[CH]([C](=O)[N])[NH]([C](=O))',
                resid_ca=0,
                resid_env=[3, 5, 6]):
    neutral_bb = '[CH](C(=O)[OH])[NH2]'
    m_neut_bb = read_smiles(neutral_bb,
                            explicit_hydrogen=True,
                            strict=False)
    neutral_ca = 0
    #tag_mol_graph(m_neut_bb)
    #print_mol(m_neut_bb)
    m_resid_bb = read_smiles(resid_bb,
                             explicit_hydrogen=True,
                             strict=False)

    m_aa = read_smiles(in_aa_string,
                       explicit_hydrogen=True,
                       strict=False)
    tag_mol_graph(m_aa)

    # Find the atom linked to CA and remove the whole backbone
    aa_ca = None
    try:
        aa2neutral = get_subgraph_match(m_aa, m_neut_bb)
    except ValueError:
        # This should be glycine!
        assert get_subgraph_match(m_aa, read_smiles('[CH2](C(=O)[OH])[NH2]', strict=False)) is not None
        aa2neutral = get_subgraph_matches(m_aa, m_neut_bb)[0]
    if aa2neutral is None:
        pro_m_neut_bb = read_smiles('[CH](C(=O)[OH])[NH]', explicit_hydrogen=True, strict=False)
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
    reload_map = get_subgraph_match(m_resid_aa, read_smiles(smiles_str, strict=False))

    return smiles_str, [reload_map[i] for i in resid_env]

def selection_to_graph(sele):
    # Fix missing informations
    u = sele.universe
    if not hasattr(u.atoms, 'elements'):
        u.add_TopologyAttr('element',
                           [guess_atom_element(at.name) for at in u.atoms])

    if not hasattr(u.atoms, 'bonds') or sum([len(at.bonds) for at in sele]) == 0:
        # Create a VdW table that only relies on element information
        u.atoms.guess_bonds()
        # if hasattr(u.atoms, 'types'):
        #     vdw_table = {}
        #     for at in u.atoms:
        #         if at.type not in vdw_table:
        #             vdw_table[at.type] = mda.topology.tables.vdwradii[at.element.upper()]
        # else:
        #     vdw_table = None

        # if not hasattr(u.atoms, 'bonds'):
        #     u.atoms.guess_bonds(vdwradii=vdw_table)
        # else:
        #     sele.atoms.guess_bonds(vdwradii=vdw_table)

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
