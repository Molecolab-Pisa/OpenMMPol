from pysmiles import read_smiles, write_smiles
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
import json
import logging

logger_pysmiles = logging.getLogger('pysmiles')
logger_pysmiles.setLevel(level=logging.ERROR)

logger = logging.getLogger(__name__)

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
                    logger.error("Can't find fragment with fragment id {:d}".format(self.assigned_fragid[i]))
                    return None
                try:
                    at += [frag.atomtypes[self.assigned_atomid[i]]]
                except KeyError:
                    logger.warning("Atom {:d} assigned as atom {:d} of residue {:d} (...) hasn't an atomtype".format(i, 
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
            logger.critical("A database should be set before starting the assignament")
            raise
        
        
        if self.mol_sub_graphs is None:
            mgs = [self.mol_graph]
        else:
            mgs = self.mol_sub_graphs

        for mg in mgs:

            logger.info("Assigning molecule with {:d} nodes and {:d} arcs".format(len(mg.nodes), len(mg.edges)))

            for frag in self.db:
                
                # Skip if the number of atoms in fragment is larger than the
                # number of atoms in current molecule
                if frag.natoms > len(mg.nodes):
                    continue

                logger.debug("Searching for residue {:s} (priority {:d})".format(frag.name,
                                                                                 frag.priority))
                aa_matches = get_subgraph_matches(mg, 
                                                  frag.frag_graph,
                                                  environment_nodes_sub=frag.env_atm,
                                                  bridge_atm=frag.bridge_atm,
                                                  full_chemical_equivalence=False)
                if len(aa_matches) > 0:
                    logger.debug("Found {:d} {:s} resids".format(len(aa_matches), frag.name))

                for m in aa_matches:
                    candidate_resid = self.get_next_free_resid()

                    assignable_atoms = [i for i in m if m[i] not in frag.env_atm]

                    if any([self.is_assigned(a) for a in assignable_atoms]):
                        logger.debug("The structure is already assigned!")
                        res_to_unassign = []
                        at_to_unassign = []
                        new_assigned_at = []
                        for at in assignable_atoms:
                            if self.is_assigned(at):
                                logger.debug("Atom {:03d} is in resid {:03d}".format(at, self.assigned_resid[at]))
                                res_to_unassign += [self.assigned_resid[at]]
                                at_to_unassign += [at]
                            else:
                                logger.debug("Atom {:03d} is not in a residue!".format(at))
                                new_assigned_at += [at]
                        
                        if len(set(res_to_unassign)) > 1:
                            logger.critical("Multiple residues have to be unassigned!")
                            raise NotImplementedError
                        if len(new_assigned_at) == 0:
                            logger.debug("No improvement here.")
                            continue
                        else:
                            if np.count_nonzero(self.assigned_resid == res_to_unassign[0]) == len(at_to_unassign):
                                logger.critical("Reassigning would be needed here! This probably means that the priority of your fragment is not correct")
                                raise NotImplementedError
                            else:
                                logger.critical("Partial reassignement is not expected!")
                                raise NotImplementedError
                        continue

                    logger.debug("Assigning {:d} atoms to residue {:d}".format(len(m), candidate_resid))
                    for i in m:
                        if m[i] not in frag.env_atm:
                            if self.is_assigned(i):
                                logger.critical("Attempting unhandled ressignament")
                                raise NotImplementedError
                            else:
                                self.assigned_resid[i] = candidate_resid
                                self.assigned_fragid[i] = frag.fragment_id
                                self.assigned_atomid[i] = m[i]

            logger.info("Assignament completed for molecule. Atoms assigned {:d}/{:d}".format(len(mg.nodes), sum(self.assigned_resid[mg.nodes] > 0)))

    def learn_from(self, use_name=True, use_type=True):
        """Use the assignement and the atom types/names present in the
        structure to improve the database and verify it against the new
        set of data""" 
        if use_type:
            self.learn_type()
        if use_name:
            self.learn_names()

    def learn_type(self):
        """Use the assignement and the atom types present in the
        structure to improve the database and verify it against the new
        set of data""" 
        
        if not hasattr(self.u.atoms[0], 'type'):
            logger.critical("Cannot learn, since your input file has no atom types")
            raise ValueError

        for i in range(self.natoms):
            if self.is_assigned(i):
                logger.debug("Correct atom type for atom {:03d} is {}".format(i, self.u.atoms[i].type))
                # Find the correct fragment
                frag = self.db[self.assigned_fragid[i]]
                logger.info("Atom is assigned to frag {:s} - id {:d}".format(frag.name, self.assigned_atomid[i]))
                if frag.has_atomtypes(self.assigned_atomid[i]):
                    try:
                        assert self.u.atoms[i].type == frag.atomtypes[self.assigned_atomid[i]]
                        logger.info("Already in DB with consistent atomtype")
                    except AssertionError:
                        logger.critical("""Inconsistent information between learning input 
                        and database for atom {:d} (atom {:d} of residue {:s})""".format(i, self.assigned_atomid[i], frag.name))
                        raise
                else:
                    frag.set_atomtype(self.assigned_atomid[i], self.u.atoms[i].type)
                    logger.info("Inserting in db")

            else:
                logger.info("Atom {:03d} is unassigned".format(i))
    
    def learn_names(self):
        """Use the assignement and the atom names present in the
        structure to improve the database and verify it against the new
        set of data""" 
        
        if not hasattr(self.u.atoms[0], 'name'):
            logger.critical("Cannot learn, since your input file has no atom names")
            raise ValueError

        for i in range(self.natoms):
            if self.is_assigned(i):
                logger.debug("Correct atom name for atom {:03d} is {}".format(i, self.u.atoms[i].name))
                # Find the correct fragment
                frag = self.db[self.assigned_fragid[i]]
                logger.info("Atom is assigned to frag {:s} - id {:d}".format(frag.name, self.assigned_atomid[i]))
                if frag.has_atomnames(self.assigned_atomid[i]):
                    try:
                        assert self.u.atoms[i].name == frag.atomnames[self.assigned_atomid[i]]
                        logger.info("Already in DB with consistent atomname")
                    except AssertionError:
                        logger.critical("""Inconsistent information between learning input 
                        and database for atom {:d} (atom {:d} of residue {:s} - name {:s})""".format(i, 
                                                                                                     self.assigned_atomid[i], 
                                                                                                     frag.name,
                                                                                                     frag.atomnames[self.assigned_atomid[i]]))
                        raise
                else:
                    frag.set_atomname(self.assigned_atomid[i], self.u.atoms[i].name)
                    logger.info("Inserting in db")

            else:
                logger.info("Atom {:03d} is unassigned".format(i))

    def save_tinker_xyz(self, fname, header_str=None):
        "Save the present structure with assigned atomtypes as a Tinker xyz file."
        if header_str is None:
            header = 'Generated by TopologyAssignament by Mattia Bondanza'
        else:
            header = header_str.replace('\n', ' ')
            
        assigned_atomtypes = self.get_assigned_atomtypes()
        if assigned_atomtypes is None:
            logger.error("Cannot continue without atom types")
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
        s = '\n'
        # First write completely assigned molecules
        mol_stat = {}
        mol_repo = {}
        for i, mol in enumerate(self.mol_sub_graphs):
            mol_stat[i] = np.sum(self.assigned_resid[mol.nodes]>=0) / len(mol.nodes)
            resids = []
            s0 = 'Molecule {:d} ({:d} atoms) {:.1f} % assigned: '.format(i, len(mol.nodes), mol_stat[i]*100)
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
            mol_repo[i] = s0 + s1 + '\n'
        
        sorted_mol_repo = sorted(mol_stat, key=lambda k: mol_stat[k], reverse=True)
        for i in sorted_mol_repo:
            s += mol_repo[i]
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

    def __contains__(self, item):
        if type(item) != Fragment:
            logger.critical("Only fragments can be contained in database")
            raise TypeError
        else:
            for fr in self:
                # check if current fragment is isomorphic to some other fragments
                if nx.faster_could_be_isomorphic(item.frag_graph, fr.frag_graph):
                    isom = nx.is_isomorphic(item.frag_graph, 
                                            fr.frag_graph, 
                                            node_match=compare_atoms, 
                                            edge_match=compare_bonds)
                    if isom:
                        return True
            return False

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
        if not frag in self:
            self.db += [frag]
            fragment_id = self.get_min_free_id()
            self.id_to_index[fragment_id] = len(self.db) - 1
            self.db[-1].fragment_id = fragment_id
        else:
            logger.warning("An isomorphic fragment is already present in db, skipping {}".format(frag.name))

    def remove_fragment(self, frag):
        fragment_id = frag.fragment_id
        try:
            index = self.id_to_index[fragment_id]
        except IndexError:
            logger.error("Trying to remove an unexistent fragment")
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
                 atomnames=None,
                 atomtypes=None,
                 default_resname='UNK',
                 restype='unknown',
                 default_oln='-',
                 priority=None):

        self.name = name
        self.smiles_str = smiles_str
        if '/' in self.smiles_str or '\\' in self.smiles_str:
            logger.warning("Stereochemical information will be discarded!")
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
            self.atomtypes = dict()
        else:
            self.atomtypes = atomtypes
        
        if atomnames is None:
            self.atomnames = dict()
        else:
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
    
    def has_atomnames(self, iatm=None):
        if iatm is None:
            return len(self.atomnames) == len(self.assignable_atm)
        else:
            if iatm in self.atomnames:
                return self.atomnames[iatm] is not None
            else:
                return False

    def set_atomname(self, iatm, atomname):
        if iatm in self.atomnames and self.atomnames[iatm] is not None:
            raise ValueError("Cannot reassign atomname")
        else:
            self.atomnames[iatm] = atomname
    
    def reset_atomname(self, iatm):
        self.atomnames[iatm] = None

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
            logger.error("New backbone should be provide as a glycine residue")
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
        
        #print([bbfrag.frag_graph.nodes[n]['id'] for n in bbfrag.frag_graph])
        #print([target_base.nodes[n]['id'] for n in target_base])

        new_frag = nx.disjoint_union(target_base, bbfrag.frag_graph)
        i_CA = [n for n in new_frag if new_frag.nodes[n]['id'] == 'new_CA'][0]
        j_CA = [n for n in new_frag if new_frag.nodes[n]['id'] == 'old_CA'][0]
        j_CB = [n for n in new_frag if new_frag.nodes[n]['id'] == 'old_CB'][0]
        if is_proline:
            #print([new_frag.nodes[n]['id'] for n in new_frag])
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
        logger.critical("More than one subgraph match was found when only one was expected")
        raise ValueError


def selection_to_graph(sele):
    u = sele.universe
    try:
        assert hasattr(u.atoms, 'elements')
        assert all(u.atoms.elements != '')
    except AssertionError:
        logger.critical("Cannot create a molecule graph without element informations in selection")
        raise

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

def guess_elements(u):
    """Interface to the standard function to guess element for each atom
    present in MDAnalysis"""
    
    from MDAnalysis.topology.guessers import guess_atom_element

    if not hasattr(u.atoms, 'elements') or any(u.atoms.elements == ''):
        logger.info("Guessing atom elements")
        u.add_TopologyAttr('element',
                           [guess_atom_element(at.name) for at in u.atoms])

def guess_bonds(u):
    """Modified algorithm to guess bonds from coordinates, it basically
    uses standard MDAnalysis functions adapted to implement the first two
    steps of the algorithm proposed in 10.1186/1758-2946-4-26"""

    # Covalent radii from 10.1186/1758-2946-4-26
    COVALENT_RADII = {
            'H':	0.23,
            'HE':	0.93,
            'LI':	0.68,
            'BE':	0.35,
            'B':	0.83,
            'C':	0.68,
            'N':	0.68,
            'O':	0.68,
            'F':	0.64,
            'NE':	1.12,
            'NA':	0.97,
            'MG':	1.10,
            'AL':	1.35,
            'SI':	1.20,
            'P':	1.05,
            'S':	1.02,
            'CL':	0.99,
            'AR':	1.57,
            'K':	1.33,
            'CA':	0.99,
            'SC':	1.44,
            'TI':	1.47,
            'V':	1.33,
            'CR':	1.35,
            'MN':	1.35,
            'FE':	1.34,
            'CO':	1.33,
            'NI':	1.50,
            'CU':	1.52,
            'ZN':	1.45,
            'GA':	1.22,
            'GE':	1.17,
            'AS':	1.21,
            'SE':	1.22,
            'BR':	1.21,
            'KR':	1.91,
            'RB':	1.47,
            'SR':	1.12,
            'Y':	1.78,
            'ZR':	1.56,
            'NB':	1.48,
            'MO':	1.47,
            'TC':	1.35,
            'RU':	1.40,
            'RH':	1.45,
            'PD':	1.50,
            'AG':	1.59,
            'CD':	1.69,
            'IN':	1.63,
            'SN':	1.46,
            'TE':	1.47,
            'I':	1.40,
            'XE':	1.98,
            'CS':	1.67,
            'BA':	1.34,
            'LA':	1.87,
            'CE':	1.83,
            'PR':	1.82,
            'ND':	1.81,
            'PM':	1.80,
            'SM':	1.80,
            'EU':	1.99,
            'GD':	1.79,
            'TB':	1.76,
            'DY':	1.75,
            'HO':	1.74,
            'ER':	1.73,
            'TM':	1.72}
    
    if not hasattr(u.atoms, 'bonds') or sum([len(at.bonds) for at in u.atoms]) == 0:
        logger.info("Guessing bonds")
        if hasattr(u.atoms, 'types'):
            vdw_table = {}
            for at in u.atoms:
                if at.type not in vdw_table:
                    vdw_table[at.type] = COVALENT_RADII[at.element.upper()]
        else:
            logger.warning("VdW radii table cannot be created as the universe has no atom-type information")
            vdw_table = None

        if vdw_table is not None:
            # To guess bonds we are basically using a modified version of eq. (1) from
            # 10.1186/1758-2946-4-26.
            #
            # A bond is present if 
            #     d_ij < ff * (Ri + Rj) + delta
            # To be compliant with MDAnalysis scheme delta is incorparated in Ri and Rj

            logger.info("Guessing bonds with algorithm by Zhang et al. (10.1186/1758-2946-4-26)")
            delta = .4
            ff = 1.0
            u.atoms.guess_bonds(vdwradii={el: vdw_table[el] + delta/(2.0 * ff) for el in vdw_table},
                                fudge_factor=ff)
        else:
            logger.warning("Guessing bonds with default MDAnalysis algorithm, this could have some issues")
            u.atoms.guess_bonds()


def check_for_overconnected_atoms(u, remove=True):
    
    """This function check if there are some atoms with a number of bonds 
    too high for the element, in case it can either remove the bonds or 
    rise an error"""

    oct = {'H': 1,
           'NA': 0,
           'K': 0,
           'CL': 0,
           'MG': 0,
           'ZN': 0,
           'C': 4,
           'N': 4,
           'O': 2,
           'P': 4,
           'S': 4,
           }
    
    for a in u.atoms:
        if a.element in oct:
            if len(a.bonds) > oct[a.element.upper()]:
                dists = np.zeros(len(a.bonds))
                for i, b in enumerate(a.bonds):
                    dists[i] = np.linalg.norm(b[0].position - b[1].position)
                if remove:
                    logger.warning("Atom {} ({}) has {} bonds (max {}). The {} longest bond(s) will be removed".format(a.index, 
                                                                                                                    a.element, 
                                                                                                                    len(a.bonds), 
                                                                                                                    oct[a.element.upper()],
                                                                                                                    len(a.bonds)-oct[a.element]))
                else:
                    logger.error("Atom {} ({}) has {} bonds (max {}). The {} longest bond(s) will be removed".format(a.index, 
                                                                                                                     a.element, 
                                                                                                                     len(a.bonds), 
                                                                                                                     oct[a.element.upper()],
                                                                                                                     len(a.bonds)-oct[a.element]))
                    raise ValueError

                bond_to_del = np.argsort(dists)[::-1][:len(a.bonds)-oct[a.element.upper()]]
                for i in bond_to_del:
                    b = a.bonds[i]
                    if remove:
                        logger.info("Removing bond between atom {:d} ({:s}) and {:d} ({:s}) with length of {:.2f}A".format(b[0].index,
                                                                                                                        b[0].element,
                                                                                                                        b[1].index,
                                                                                                                        b[1].element,
                                                                                                                        dists[i]))
                u.delete_bonds(a.bonds[bond_to_del])
