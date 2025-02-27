import numpy as np
import MDAnalysis as mda
import json
import argparse
from scipy.spatial.distance import cdist
from time import time
import hashlib

def match_atoms(uni1, uni2):
    """This function takes two universes and return a 
    dictionary {at1@uni1: at1@uni2, ... } containing 
    the atoms that match (in general coordinates and 
    element are checked)."""

    tstart = time()

    DIST_THR = 1e-2

    #Define which system is larger
    if len(uni1.atoms) >= len(uni2.atoms):
        uniL = uni1
        uniS = uni2
    else:
        uniL = uni2
        uniS = uni1
    
    distmatrix = cdist(uniS.atoms.positions, 
                       uniL.atoms.positions)
    matchlist = {}
    for i, at in enumerate(uniS.atoms):
        mindist = np.min(distmatrix[i,:])
        if mindist < DIST_THR:
            # Ok this atom has a match!
            match_at = np.argmin(distmatrix[i,:])
            # TODO CHECK
            if uniL == uni1:
                matchlist[match_at] = at
            else:
                matchlist[i] = uniL.atoms[match_at]
        else:
            print("No matching atom found for atom {}".format(at))
    
    print("Time elap in matchlist: {:g} s".format(time()-tstart))

    return matchlist

def input_load(fin_txyz, fin_pdb=None):
    print("Loading Tinker XYZ")
    try:
        universe = mda.Universe(fin_txyz, format='txyz')
    except FileNotFoundError:
        print("{:s} not found, cannot open input.".format(fin_txyz))
        return None
    except (ValueError, IndexError):
        print("""{:s} doesen't seems to be a correctly"""
              """formatted Tinker xyz.""".format(fin_txyz))
        return None

    print("Loading PDB")
    if fin_pdb is not None:
        try:
            universe_pdb = mda.Universe(fin_pdb, format='pdb')
        except FileNotFoundError:
            print("{:s} not found, cannot open input.".format(fin_pdb))
            return None
        except (ValueError, IndexError):
            print("""{:s} doesen't seems to be a correctly"""
                  """formatted PDB.""".format(fin_pdb))
            return None
    else:
        return universe

    if not hasattr(universe_pdb.atoms, 'elements'):
        universe_pdb.add_TopologyAttr('element', [mda.topology.guessers.guess_atom_element(at.name) for at in universe_pdb.atoms])
        
    # Now let's start dancing
    # Match the atoms of the two universes
    match_list = match_atoms(universe, universe_pdb)
    # Create a list containing residues index for each atom,
    # it is incredibly more efficient to do that here than
    # assigning later one by one.
    tstart = time()
    atom_resindex = np.zeros(len(universe.atoms))
    for at in universe.atoms:
        if at.index in match_list:
            atom_resindex[at.index] = match_list[at.index].residue.ix
    for at in universe.atoms:
        if at.index not in match_list:
            if len(at.bonds) == 1:
                atom_resindex[at.index] = atom_resindex[at.bonds[0].partner(at).index]
            else:
                print("""Atom {} does not match to any atom in PDB and """
                      """we do not understand why it is here (number of """
                      """bond different from 1).""".format(at))
                return None


    print("Atom resindex list: {:g} s".format(time()-tstart))
    tstart = time()
    # 1. create a new universe
    new_universe = mda.Universe.empty(n_atoms=len(universe.atoms),
                                      n_residues=len(universe_pdb.residues),
                                      n_segments=len(universe_pdb.segments),
                                      atom_resindex=atom_resindex,
                                      residue_segindex=universe_pdb.residues.segindices,
                                      trajectory=True)

    new_universe.add_TopologyAttr('segid', [s.segid for s in universe_pdb.segments])
    new_universe.add_TopologyAttr('resid', [r.resid for r in universe_pdb.residues])
    new_universe.add_TopologyAttr('resnum', [r.resnum for r in universe_pdb.residues])
    new_universe.add_TopologyAttr('resname', [r.resname for r in universe_pdb.residues])

    new_universe.add_TopologyAttr('element')
    new_universe.add_TopologyAttr('name')
    new_universe.add_TopologyAttr('type')
    new_universe.add_TopologyAttr('bonds', [(b[0].index, b[1].index) for b in universe.bonds])
    print("Creation of new universe: {:g} s".format(time()-tstart))
    
    tstart = time()
    # Let's start trivial
    for at in new_universe.atoms:
        if at.index in match_list:
            # From PDB
            at.element = match_list[at.index].element
            # Overwrite the tinker names, that here are just harmful
            at.name = match_list[at.index].name
            # From tinker XYZ
            at.type = universe.atoms[at.index].type
    print("Matched atoms: {:g} s".format(time()-tstart))

    tstart = time()
    # Then think to the unmatched atoms
    for at in new_universe.atoms:
        if at.index not in match_list:
            # Only go ahead if the atom has just one bond!
            # Already checked before
            guessed_element = mda.topology.guessers.guess_atom_element(universe.atoms[at.index].name)
            # Check that it looks like an hydrogen...
            if guessed_element == 'H':
                # Guess that it is an hydrogen atom
                # It was probably added during input preparation so it is
                # allowed to be here.
                at.element = 'H'
                at.name = 'H'
                at.type = universe.atoms[at.index].type
            elif guessed_element == 'O':
                # Ok that's a carbonyl oxygen let's see if it is a terminal one.
                partn = at.bonds[0].partner(at)
                nox = 0
                nca = 0
                for b in partn.bonds:
                    if b.partner(partn).element == 'O' \
                       or b.partner(partn) == at:
                        nox += 1
                    elif b.partner(partn).element == 'C':
                        nca += 1
                if nox == 2 and nca == 1:
                    # Ok it actually seems a terminal carbon the oxygen
                    # is allowed to stay there
                    at.element = 'O'
                    at.name = 'O'
                    at.type = universe.atoms[at.index].type
            else:
                print("""Atom {} does not match to any atom in PDB and """
                      """we do not understand why it is here (element """
                      """or connectivity unrecognized).""".format(at))
                return None

    new_universe.load_new(universe.atoms.positions[np.newaxis,:,:])

    print("Unmatched atoms: {:g} s".format(time()-tstart))
    return new_universe

def txyz_writer(atgrp, outfile):
    with open(outfile, 'w+') as f:
        print("{:d} Automatically generated by OMMP-IP".format(len(atgrp)), file=f)
        for i, at in enumerate(atgrp):
            print("{:8d} {:5s} {:12.6f} {:12.6f} {:12.6f} {:5d}".format(i+1,
                                                                        at.element,
                                                                        at.position[0],
                                                                        at.position[1],
                                                                        at.position[2],
                                                                        int(at.type)),
                  end='', file=f)
            for bnd in at.bonds:
                partn = bnd.partner(at)
                if partn in atgrp:
                    print("{:8d} ".format(np.argwhere(atgrp.indices == partn.index).flatten()[0]+1), end='', file=f)
            print('', file=f)

def file_md5sum(fname):
    with open(fname, 'r') as f:
        md5sum = hashlib.md5(f.read().encode()).hexdigest()
    return md5sum

def ommp_ip_main():
    parser = argparse.ArgumentParser(prog='ommp-ip',
                        description="""OpenMMPol Input Preparation program, 
                                    this script should assist the process 
                                    creating input for OpenMMPol library 
                                    starting from standard structure file.
                                    It is mainly oriented to protein 
                                    complexes currently.""",
                        epilog='Good luck with your calculation!')

    parser.add_argument('-i', '--input-file',
                        required=True,
                        metavar='<TINKER_XYZ_INPUT_FILE>',
                        dest='input_file',
                        help="""Tinker-xyz/arc file containing both QM and MM 
                                parts of the system""")
    parser.add_argument('-p', '--pdb-input-file',
                        required=False,
                        metavar='<PDB_INPUT_FILE>',
                        dest='pdb_input_file',
                        default = None,
                        help="""PDB input file used to generate the Tinker-xyz 
                                file, this is used to extract extra information
                                about atom names (standard) and residues that
                                are not present in Tinker-xyz, in order to
                                simplify selection strings""")
    parser.add_argument('-d', '--prm-file',
                        required=False,
                        metavar='<TINKER_PRM_INPUT_FILE>',
                        dest='prm_input_file',
                        default = None,
                        help="""Tinker parameter file, it should be the one used
                                to create the Tinker XYZ file in input.""")
    parser.add_argument('-q', '--qm-selection',
                        required=False,
                        metavar='<QM_SELECTION_STRING>',
                        dest='qm_sele',
                        default='not all',
                        help="""Selection string used to select atoms in the system
                                that will be in the QM part""")
    parser.add_argument('-f', '--frozen-selection',
                        required=False,
                        metavar='<FROZEN_SELECTION_STRING>',
                        dest='frozen_sele',
                        default='not all',
                        help="""Selection string used to select atoms in the system
                                that are frozen, that is fixed during a MD/optimization""")
    parser.add_argument('-o', '--output-basename',
                        required=False,
                        metavar='<OUTPUT_BASENAME>',
                        dest='output_basename',
                        default='out',
                        help="""All output file will be created using this 
                                as basename. Default is 'out'""")

    args = parser.parse_args()

    mmtxyz_path = args.output_basename+'_mm.arc'
    qmxyz_path = args.output_basename+'_qm.xyz'
    json_si_path = args.output_basename+'_si.json'
    mmpdb_path = args.output_basename+'_mm.pdb'
    frozenpdb_path = args.output_basename+'_frozen.pdb'

    qmmm_sys = input_load(args.input_file, args.pdb_input_file)
    #qmmm_sys.atoms.write(args.output_basename+'_full_system.pdb')dd
    if qmmm_sys is None:
        exit(1)

    qm_sys = qmmm_sys.select_atoms(args.qm_sele)
    mm_sys = qmmm_sys.atoms - qm_sys
    print("The whole system has {:d} atoms.".format(len(qmmm_sys.atoms)))
    print("{:d} atoms selected for QM part".format(len(qm_sys)))
    print("{:d} atoms selected for MM part".format(len(mm_sys)))

    # Frozen atoms
    frozen = qmmm_sys.select_atoms(args.frozen_sele)
    notfrozen = qmmm_sys.atoms - frozen
    mm_frozen = []
    qm_frozen = []
    if len(frozen) > 0:
        for at in frozen:
            if at in mm_sys:
                mm_frozen += [int(np.argwhere(at.index == mm_sys.indices).flatten()[0])+1]
            else:
                qm_frozen += [int(np.argwhere(at.index == qm_sys.indices).flatten()[0])+1]
    print("{:d} atoms are frozen.".format(len(frozen)))
    print("{:d} frozen are in QM part".format(len(qm_frozen)))
    print("{:d} frozen are in MM part".format(len(mm_frozen)))

    
    # Search for required link atoms
    linkatoms = []
    for at in qm_sys:
        for bnd in at.bonds:
            partn = bnd.partner(at)
            if partn not in qm_sys:
                print("Link atom required between atoms {} and {}".format(at, partn))
                qmat = at
                qmat_i = np.argwhere(qmat.index == qm_sys.indices).flatten()[0]
                mmat = partn
                mmat_i = np.argwhere(mmat.index == mm_sys.indices).flatten()[0]
                la_pos = qmat.position + (mmat.position - qmat.position) / np.linalg.norm(mmat.position - qmat.position) * 1.0
                la_element = 'H'
                la_name = 'H'
                # If possible use the same kind of another H bonded to the QM atom
                la_type = 0
                for qmbnd in qmat.bonds:
                    p = qmbnd.partner(qmat)
                    if p in qm_sys and p.element == 'H':
                        la_type = p.type
                if qmat.element == 'H':
                    print("Creating a linkatom between an hydrogen (QM) and an MM atom is a very ill choice.")
                    print("Atoms involved QM:{} MM:{}".format(qmat, mmat))
                    exit(1)
                if la_type == 0:
                    print("""Creating a link atom with a QM atom without any other H. """
                          """Link atom type cannot be guessed, it is left to 0, you should"""
                          """ likely chage it before running your calculations""")
                linkatoms += [{'mm': mmat_i+1, 
                               'qm': qmat_i+1, 
                               'la': len(qm_sys)+len(linkatoms)+1, 
                               'pos': la_pos,
                               'la_element': la_element,
                               'la_name': la_name,
                               'la_type': la_type,
                               'la_dist': None}]
    if len(linkatoms) > 0:
        # Create a universe just for link atoms
        la_universe = mda.Universe.empty(n_atoms=len(linkatoms),
                                        trajectory=True)
        la_universe.add_TopologyAttr('element', [la['la_element'] for la in linkatoms])
        la_universe.add_TopologyAttr('name', [la['la_name'] for la in linkatoms])
        la_universe.add_TopologyAttr('type', [la['la_type'] for la in linkatoms])
        la_universe.load_new(np.array([la['pos'] for la in linkatoms])[np.newaxis,:,:])
        # and merge it with the QM subsystem
        qm_and_la_sys = mda.Merge(qm_sys, la_universe.atoms)
    else:
        qm_and_la_sys = qm_sys

    # 1. Output files needed for the calculation

    txyz_writer(mm_sys, mmtxyz_path)
    qm_and_la_sys.atoms.write(qmxyz_path)

    # 2. Assemble JSON file
    json_si_data = {}
    json_si_data['name'] = ''
    json_si_data['description'] = 'Generated by OMMP-IP'
    json_si_data['version'] = '0.5.0'
    json_si_data['xyz_file'] = {}
    json_si_data['xyz_file']['path'] = mmtxyz_path
    json_si_data['xyz_file']['md5sum'] = file_md5sum(mmtxyz_path)
    json_si_data['prm_file'] = {}
    if args.prm_input_file is not None:
        json_si_data['prm_file']['path'] = args.prm_input_file
        json_si_data['prm_file']['md5sum'] = file_md5sum(args.prm_input_file)
    else:
        json_si_data['prm_file']['path'] = "<PUT PRM FILE HERE>"
        print("Since no prm file has been provide in input, please add it in the JSON file")
    json_si_data['verbosity'] = 'low'
    if len(mm_frozen) > 0:
        json_si_data['frozen_atoms'] = mm_frozen
    if len(qm_and_la_sys.atoms) > 0:
        json_si_data['qm'] = {}
        json_si_data['qm']['qm_atoms'] = [at.element for at in qm_and_la_sys.atoms]
        json_si_data['qm']['qm_coords'] = [[float(at.position[0]),
                                            float(at.position[1]),
                                            float(at.position[2])] for at in qm_and_la_sys.atoms]
        json_si_data['qm']['qm_atom_types'] = [int(at.type) for at in qm_and_la_sys.atoms]
        json_si_data['qm']['prm_file'] = {}
        if args.prm_input_file is not None:
            json_si_data['qm']['prm_file']['path'] = args.prm_input_file
            json_si_data['qm']['prm_file']['md5sum'] = file_md5sum(args.prm_input_file)
        else:
            json_si_data['qm']['prm_file']['path'] = "<PUT PRM FILE HERE>"
            print("Since no prm file has been provide in input, please add it in the JSON file")
        if len(qm_frozen) > 0:
            json_si_data['qm']['qm_frozen_atoms'] = qm_frozen
    if len(linkatoms) > 0:
        json_si_data['link_atoms'] = []
        for la in linkatoms:
            json_si_data['link_atoms'] += [{'MM_id': int(la['mm']),
                                            'QM_id': int(la['qm']),
                                            'LA_id': int(la['la'])}]
            if la['la_dist'] is not None:
                json_si_data['link_atoms'][-1]['bond_length'] = float(la['la_dist'])
            #if ...:
            #     json_si_data['link_atoms'][-1]['eel_remove'] =

    # Write on file
    with open(json_si_path, 'w+') as f:
        print(json.dumps(json_si_data, indent=4), file=f)

    # 3. Assemble utility files
    mm_sys.atoms.write(mmpdb_path)
    if len(frozen.atoms) > 0:
        frozen.atoms.write(frozenpdb_path)

if __name__ == '__main__':
    ommp_ip_main()



    
