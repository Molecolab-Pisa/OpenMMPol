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
    """Class for performing assignament between a molecular structure and a database."""
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
        """Set the database to be used in the assignament"""
        self.db = db

    def get_assigned_atomtypes(self):
        """For each atom in the molecular structure, returns the assigned atomtypes"""
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
                aa_matches = get_subgraph_matches(mg, 
                                                  frag.frag_graph,
                                                  environment_nodes_sub=frag.env_atm,
                                                  bridge_atm=frag.bridge_atm,
                                                  full_chemical_equivalence=False)
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
                                print("Atom {:03d} is not in a residue!".format(at))
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
                            #print(np.count_nonzero(self.assigned_resid == res_to_unassign[0]))
                            #print(len(at_to_unassign))
                            #print(len(new_assigned_at))
                            if np.count_nonzero(self.assigned_resid == res_to_unassign[0]) == len(at_to_unassign):
                                print("Ok reassigning is needed here!")
                                raise NotImplementedError
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
        """Use the assignement and the atom types/names present in the
        structure to improve the database and verify it against the new
        set of data""" 
        if use_name:
            raise NotImplementedError("Learning atomnames is still not implemented")
        if use_type:
            self.learn_type()

    def learn_type(self):
        """Use the assignement and the atom types present in the
        structure to improve the database and verify it against the new
        set of data""" 
        
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
        "Save the present structure with assigned atomtypes as a Tinker xyz file."
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
        """Write a log for the structure identification"""
        s = ''
        for i, mol in enumerate(self.mol_sub_graphs):
            s += 'Molecule {:d} ({:d} atoms): '.format(i, len(mol.nodes))
            resids = []
            s1 = ''
            nuna = 0
            for atid in sorted(mol.nodes):
                if self.assigned_resid[atid] < 0:
                    if not s1.endswith('UNASSIGNED'):
                        if len(s1) > 0:
                            s1 += '-'
                        s1 += 'UNASSIGNED'
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
    """Database class, currently is very similar to an array of Fragment object sorted by priority"""
    def __init__(self, dbfile=None):
        """Database init, it can be initialized as an empty db without arguments or from a json file"""
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
        """Save database as json file"""
        json_out = {}
        json_out['res_assignament'] = [f.asdict() for f in self.db]
        with open(fout, 'w') as f:
            print(json.dumps(json_out, indent=4), file=f)

    def get_min_free_id(self):
        """Return the lowest id for inserting a new residue"""
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
    """A Fragment is litterally a piece of molecule (or even a whole molecule)
    that is searched during the assignament process. It is represented by a 
    smiles string.
    Each fragment could have "environment atoms" (env_atm) which
    should be there in the topology for getting a match, but are not assigned.
    Each fragment could have "bridge atoms" (bridge_atm) which are allowed
    to be connected with atoms that are not matched (neither as env_atm nor as
    regular atom.
    Examples: 
        Water
            smiles     "O"
            env_atm    []
            bridge_atm []
    
        Glycine (residue inside a protein chain)
            smiles     "O=C([N])CN[C]=O"
            (id)       "0 1  2  34 5  6"
            env_atm    [2,5,6] (to enforce that is connected with amide bonds)
            bridge_atm [2, 5] (the bridge to other residues, all the other
                               atoms are connected only to atoms of the matched 
                               substructure)"""
    def __init__(self,
                 name,
                 smiles_str,
                 env_atm=[],
                 bridge_atm=[],
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
        self.bridge_atm = bridge_atm
        self.natoms = len(self.frag_graph.nodes)

        if default_resname is not None:
            self.resname = default_resname
        else:
            self.resname = name
        self.oln = default_oln
        if priority is None:
            self.priority = self.natoms#10 * (self.natoms-len(self.env_atm)) + len(self.env_atm)
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
            return len(self.atomtypes) == len(self.assignable_atm)
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
        d['bridge_atm'] = list(self.bridge_atm)
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
                   bridge_atm=d['bridge_atm'],
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

    def is_aminoacid_residue(self):
        minimal_backbone = '[CH]([C])([C](=O))[N]'
        bb_graph = read_smiles(minimal_backbone,
                               explicit_hydrogen=True,
                               strict=False)
        tag_mol_graph(bb_graph)

        match = get_subgraph_matches(self.frag_graph, bb_graph)

        if len(match) > 1:
            # Threonine and serine
            thr_ser_minimal_backbone = '[CH]([C](OH))([C](=O))[N]'
            thr_ser_bb_graph = read_smiles(thr_ser_minimal_backbone,
                                           explicit_hydrogen=True,
                                           strict=False)
            tag_mol_graph(thr_ser_bb_graph)
            match = get_subgraph_match(self.frag_graph, thr_ser_bb_graph)
            is_thr_or_ser = match is not None
            return is_thr_or_ser

        elif len(match) == 1:
            return True
        else:
            # Glycine
            return self.is_glycine_residue()
    
    def is_glycine_residue(self):
        glycine_minimal_backbone = '[CH2]([C](=O))[N]'
        glycine_bb_graph = read_smiles(glycine_minimal_backbone,
                                        explicit_hydrogen=True,
                                        strict=False)
        tag_mol_graph(glycine_bb_graph)
        match = get_subgraph_match(self.frag_graph, glycine_bb_graph)
        is_glycine = match is not None
        return is_glycine

    def get_CA_index(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None

        minimal_backbone = '[CH]([C](=O))[N]'
        CA_idx = 0

        bb_graph = read_smiles(minimal_backbone,
                               explicit_hydrogen=True,
                               strict=False)
        tag_mol_graph(bb_graph)

        match = get_subgraph_matches(self.frag_graph, bb_graph)

        for i in match[0]:
            if match[0][i] == CA_idx:
                return i
        return None
    
    def get_C_index(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None

        minimal_backbone = '[CH]([C](=O))[N]'
        C_idx = 1

        bb_graph = read_smiles(minimal_backbone,
                               explicit_hydrogen=True,
                               strict=False)
        tag_mol_graph(bb_graph)

        match = get_subgraph_matches(self.frag_graph, bb_graph)

        for i in match[0]:
            if match[0][i] == C_idx:
                return i
        return None
    
    def get_O_index(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None

        minimal_backbone = '[CH]([C](=O))[N]'
        O_idx = 2

        bb_graph = read_smiles(minimal_backbone,
                               explicit_hydrogen=True,
                               strict=False)
        tag_mol_graph(bb_graph)

        match = get_subgraph_matches(self.frag_graph, bb_graph)

        for i in match[0]:
            if match[0][i] == O_idx:
                return i
        return None
    
    def get_N_index(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None

        minimal_backbone = '[CH]([C](=O))[N]'
        N_idx = 3

        bb_graph = read_smiles(minimal_backbone,
                               explicit_hydrogen=True,
                               strict=False)
        tag_mol_graph(bb_graph)

        match = get_subgraph_matches(self.frag_graph, bb_graph)

        for m in match:
            for i in m:
                if m[i] == N_idx and i not in self.env_atm:
                    return i
        return None
    
    def get_CB_index(self):
        try:
            assert self.is_aminoacid_residue() and not self.is_glycine_residue()
        except AssertionError:
            return None

        minimal_backbone = '[CH]([C])([C](=O))[N]'
        CB_idx = 1

        bb_graph = read_smiles(minimal_backbone,
                               explicit_hydrogen=True,
                               strict=False)
        tag_mol_graph(bb_graph)

        match = get_subgraph_matches(self.frag_graph, bb_graph)
        if len(match) > 1:
            minimal_oh_backbone = '[CH]([C]O)([C](=O))[N]'
            CB_idx = 1
            
            bb_graph = read_smiles(minimal_oh_backbone,
                                   explicit_hydrogen=True,
                                   strict=False)
            tag_mol_graph(bb_graph)

            match = get_subgraph_matches(self.frag_graph, bb_graph)

        for i in match[0]:
            if match[0][i] == CB_idx:
                return i
        return None
    
    def get_proline_CD_index(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None
        
        i_N = self.get_N_index()
        i_CA = self.get_CA_index()
        
        graph = self.frag_graph.copy()
        # 2. Try to disconnect the N terminal side
        to_disconnect = []
        for ed in graph.edges(i_N):
            if ed[0] == i_N:
                i_conn = ed[1]
            else:
                i_conn = ed[0]
            if i_conn == i_CA:
                continue
            try:
                test_graph = graph.copy()
                test_graph.remove_edge(i_N, i_conn)
                assert len(list(nx.connected_components(test_graph))) != 1
            except AssertionError:
                return i_conn
        return None

    def get_sidechain(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None

        if self.is_glycine_residue():
            gly_side = read_smiles('[H]', 
                                   explicit_hydrogen=True, 
                                   strict=False)
            tag_mol_graph(gly_side)
            return gly_side

        else:
            ca_id = self.get_CA_index()
            cb_id = self.get_CB_index()

            # graph = read_smiles(self.smiles_str,
            #                     explicit_hydrogen=True,
            #                     strict=False)
            # tag_mol_graph(graph)
            graph = self.frag_graph.copy()
            graph.remove_edge(ca_id, cb_id)
            
            side = nx.node_connected_component(graph, cb_id)
            side = nx.subgraph(self.frag_graph, side)
            if len(side) == len(self.frag_graph):
                # Proline should be treated differently
                raise NotImplementedError("No idea of how to define sidechain for a proline")
            return side
    
    def get_backbone(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None

        minimal_backbone = '[C]([C](=O))[N]'
        ox_id = 2

        bb_graph = read_smiles(minimal_backbone,
                                explicit_hydrogen=True,
                                strict=False)
        tag_mol_graph(bb_graph)
        matches = get_subgraph_matches(self.frag_graph, bb_graph)
        spurious_matches = []
        for ii, i in enumerate(matches):
            for n in i:
                if i[n] == ox_id:
                    if len(self.frag_graph.edges(n)) > 1:
                        # This is not a carbonyl
                        spurious_matches += [ii]
        for i in spurious_matches:
            matches.pop(i)

        if len(matches) != 1:
            return None
        else:
            return nx.subgraph(self.frag_graph, matches[0].keys())

    def get_base_struct(self):
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None
        
        i_O = self.get_O_index()
        i_C = self.get_C_index()
        i_N = self.get_N_index()
        i_CA = self.get_CA_index()
        # 1. Disconnect what is connected to CO
        graph = self.frag_graph.copy()
        to_disconnect = []
        for ed in graph.edges(i_C):
            if ed[0] == i_C:
                i_conn = ed[1]
            else:
                i_conn = ed[0]
            if i_conn != i_CA and i_conn != i_O:
                to_disconnect += [i_conn]
        for i_conn in to_disconnect:
            graph.remove_edge(i_C, i_conn)
        
        # 2. Disconnect the N terminal side
        to_disconnect = []
        for ed in graph.edges(i_N):
            if ed[0] == i_N:
                i_conn = ed[1]
            else:
                i_conn = ed[0]
            if i_conn == i_CA:
                continue
            try:
                test_graph = graph.copy()
                test_graph.remove_edge(i_N, i_conn)
                assert len(list(nx.connected_components(test_graph))) != 1
            except AssertionError:
                # This is the proline loop carbon
                continue
            to_disconnect += [i_conn]

        for i in to_disconnect:
            graph.remove_edge(i_N, i)
        
        base = nx.subgraph(self.frag_graph, nx.node_connected_component(graph, i_CA))
        return base
        
    def replace_aminoacid_backbone(self, new_backbone, new_backbone_env_atm, new_backbone_bridge_atm):
        # 0. Check input
        try:
            assert self.is_aminoacid_residue()
        except AssertionError:
            return None
        
        bbfrag = Fragment("new_bb",
                        new_backbone,
                        env_atm = new_backbone_env_atm,
                        bridge_atm = new_backbone_bridge_atm)
        
        try:
            assert bbfrag.is_glycine_residue()
        except AssertionError:
            print("New backbone should be provide as a glycine residue")
            return None

        # 1.1 If the aminoacid is glycine, just skip the whole thing
        #     and return the new_backbone
        if self.is_glycine_residue():
            return new_backbone, new_backbone_env_atm, new_backbone_bridge_atm

        # 1. Get structural information about the current residue
        target_base = self.get_base_struct()
        j_CA = self.get_CA_index()
        j_CB = self.get_CB_index()
        j_N = self.get_N_index()
        j_CD = self.get_proline_CD_index()
        is_proline = j_CD is not None

        i_CA = bbfrag.get_CA_index()
        i_N = bbfrag.get_N_index()

        target_base.nodes[j_CA]['id'] = 'old_CA'
        target_base.nodes[j_CB]['id'] = 'old_CB'
        bbfrag.frag_graph.nodes[i_CA]['id'] = 'new_CA'
        if is_proline:
            target_base.nodes[j_N]['id'] = 'old_N'
            target_base.nodes[j_CD]['id'] = 'old_CD'
            bbfrag.frag_graph.nodes[i_N]['id'] = 'new_N'


        for i in bbfrag.env_atm:
            if type(bbfrag.frag_graph.nodes[i]['id']) == str:
                bbfrag.frag_graph.nodes[i]['id'] += 'env'
            else:
                bbfrag.frag_graph.nodes[i]['id'] = 'env'

        for i in self.env_atm:
            if i in target_base.nodes:
                if type(target_base.nodes[i]['id']) == str:
                    target_base.nodes[i]['id'] += 'env'
                else:
                    target_base.nodes[i]['id'] = 'env'

        for i in bbfrag.bridge_atm:
            if type(bbfrag.frag_graph.nodes[i]['id']) == str:
                bbfrag.frag_graph.nodes[i]['id'] += 'bridge'
            else:
                bbfrag.frag_graph.nodes[i]['id'] = 'bridge'

        for i in self.bridge_atm:
            if i in target_base.nodes:
                if type(target_base.nodes[i]['id']) == str:
                    target_base.nodes[i]['id'] += 'bridge'
                else:
                    target_base.nodes[i]['id'] = 'bridge'
        
        print([bbfrag.frag_graph.nodes[n]['id'] for n in bbfrag.frag_graph])
        print([target_base.nodes[n]['id'] for n in target_base])

        new_frag = nx.disjoint_union(target_base, bbfrag.frag_graph)
        i_CA = [n for n in new_frag if new_frag.nodes[n]['id'] == 'new_CA'][0]
        j_CA = [n for n in new_frag if new_frag.nodes[n]['id'] == 'old_CA'][0]
        j_CB = [n for n in new_frag if new_frag.nodes[n]['id'] == 'old_CB'][0]
        if is_proline:
            print([new_frag.nodes[n]['id'] for n in new_frag])
            j_CD = [n for n in new_frag if new_frag.nodes[n]['id'] == 'old_CD'][0]
            j_N = [n for n in new_frag if new_frag.nodes[n]['id'] == 'old_N'][0]
            i_N = [n for n in new_frag if new_frag.nodes[n]['id'] == 'new_N'][0]
        
        new_frag.remove_edge(j_CA, j_CB)
        new_frag.add_edge(i_CA, j_CB)

        # Remove one H on C-alpha
        for ed in new_frag.edges(i_CA):
            if ed[0] == i_CA:
                i_H = ed[1]
            else:
                i_H = ed[0]

            if new_frag.nodes[i_H]['element'] == 'H':
                break
        new_frag.remove_edge(i_CA, i_H)

        if is_proline:
            new_frag.remove_edge(j_N, j_CD)
            new_frag.add_edge(i_N, j_CD)
            for ed in new_frag.edges(i_N):
                if ed[0] == i_N:
                    i_H = ed[1]
                else:
                    i_H = ed[0]

                if new_frag.nodes[i_H]['element'] == 'H':
                    break
            new_frag.remove_edge(i_N, i_H)

        out_frag = nx.node_connected_component(new_frag, i_CA)
        out_graph = nx.subgraph(new_frag, out_frag)
        smiles_str = write_smiles(out_graph)
        sm_graph = read_smiles(smiles_str, 
                            strict=False,
                            explicit_hydrogen=True)
        tag_mol_graph(sm_graph)
        eq = get_subgraph_match(out_graph, sm_graph)
        env_atm = [eq[i] for i in out_graph if 'env' in str(out_graph.nodes[i]['id'])]
        bridge_atm = [eq[i] for i in out_graph if 'bridge' in str(out_graph.nodes[i]['id'])]

        return smiles_str, env_atm, bridge_atm

def compare_atoms(a, b, exclude_list=None, exclude_list_b=None):
    if exclude_list is not None and a['id'] in exclude_list:
        if exclude_list_b is not None and b['id'] not in exclude_list_b:
            return False

    if a['element'] == b['element']:
        return True
    else:
        return False

def compare_bonds(a, b):
    # All bonds are equal
    return True

def tag_mol_graph(g):
    nx.set_node_attributes(g, -1, 'id')

    for i in range(len(g.nodes)):
        g.nodes[i]['id'] = i

def get_subgraph_matches(big_g, sub_g,
                         exclude_list = None,
                         environment_nodes_sub = [],
                         bridge_atm = None,
                         full_chemical_equivalence = True):

    sub_iter = nx.isomorphism.GraphMatcher(big_g,
                                           sub_g,
                                           node_match=lambda a, b: compare_atoms(a, b, exclude_list, environment_nodes_sub),
                                           edge_match=compare_bonds).subgraph_isomorphisms_iter()
    
    # NOTE Matches are returned as a generator of dictionaries {N_big: n_small}
    
    _sub_iter = []
    exclude_list = []
    for match in sub_iter:
        if bridge_atm is not None:
            bridge_atm_exc = False

            for at in match:
                if match[at] in bridge_atm:
                    continue
                else:
                    for ed in big_g.edges(at):
                        bonded = ed[0]
                        if bonded == at:
                            bonded = ed[1]
                        if bonded not in match:
                            bridge_atm_exc = True
                            break
                if bridge_atm_exc:
                    break
            if bridge_atm_exc:
                break
        
        _sub_iter += [match]
        #print(">>", match)
        exclude_list += [n for n in list(match.keys()) if match[n] not in environment_nodes_sub]
        #print("exclude_list", exclude_list, "(", environment_nodes_sub, ")")
    sub_iter = _sub_iter
    
    return remove_chemical_equivalence([a for a in sub_iter], big_g, full_chemical_equivalence)

def remove_chemical_equivalence(subs, big_g, full_chemical_equivalence=True):
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
    if not full_chemical_equivalence or len(cheq) < 2:
        return [subs[i] for i in cheq]

    #   b) all the different matched atoms in the big_g are chemically equivalent
    cheq2 = {}
    chemical_equivalent_atoms = {}
    tocheck = list(cheq.keys())
    for i in tocheck:
        for j in tocheck[i:]: 
            equivalent = True
            for at_i in subs[i]:
                ref_at = subs[i][at_i]
                for at_j in subs[j]:
                    if subs[j][at_j] == ref_at:
                        break
                if at_i != at_j:
                    if not are_chemically_equivalent_atoms(big_g, at_i, at_j):
                        equivalent = False
                
                if equivalent:
                    if i in cheq2:
                        cheq2[i] += [j]
                    else:
                        cheq2[i] = [j]

    return [subs[i] for i in cheq2]

def are_chemically_equivalent_atoms(g, i, j):
    test_graph_1 = g.copy()
    test_graph_2 = g.copy()
    test_graph_1.remove_edges_from(g.edges(i))
    test_graph_2.remove_edges_from(g.edges(j))
    
    return nx.is_isomorphic(test_graph_1, test_graph_2, compare_atoms, compare_bonds)


def get_subgraph_match(big_g, sub_g, exclude_list=None):
    ms = get_subgraph_matches(big_g, sub_g, exclude_list)
    if len(ms) == 1:
        return ms[0]
    elif len(ms) == 0:
        return None
    else:
        print(ms)
        raise ValueError("More than one match was found!")


def selection_to_graph(sele):
    # Fix missing informations
    u = sele.universe
    if not hasattr(u.atoms, 'elements') or any(u.atoms.elements == ''):
        u.add_TopologyAttr('element',
                           [guess_atom_element(at.name) for at in u.atoms])

    if not hasattr(u.atoms, 'bonds') or sum([len(at.bonds) for at in sele]) == 0:
        # Create a VdW table that only relies on element information
        if hasattr(u.atoms, 'types'):
            vdw_table = {}
            for at in u.atoms:
                if at.type not in vdw_table:
                    vdw_table[at.type] = mda.guesser.tables.vdwradii[at.element.upper()]
        else:
            vdw_table = None

        u.atoms.guess_bonds(vdwradii=vdw_table)

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
