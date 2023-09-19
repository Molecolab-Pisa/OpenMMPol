import numpy as np
import MDAnalysis as mda
import json
import argparse
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

def match_atoms(uni1, uni2):
    """This function takes two universes and return a 
    dictionary {at1@uni1: at1@uni2, ... } containing 
    the atoms that match (in general coordinates and 
    element are checked)."""

    DIST_THR = 1e-4

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
                matchlist[match_at] = i
            else:
                matchlist[i] = match_at
        else:
            print("No matching atom found for atom {}".format(at))
    
    return matchlist

def input_load(fin_txyz, fin_pdb=None):
    try:
        universe = mda.Universe(fin_txyz, format='txyz')
    except FileNotFoundError:
        print("{:s} not found, cannot open input.".format(fin_txyz))
        return None
    except (ValueError, IndexError):
        print("""{:s} doesen't seems to be a correctly"""
              """formatted Tinker xyz.""".format(fin_txyz))
        return None

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

    # Match the atoms of the two universes
    match_list = match_atoms(universe, universe_pdb)
    
    # Now let's start dancing
    # 1. create a new universe
    new_universe = mda.Universe.empty(n_atoms=len(universe.atoms),
                                      n_residues=len(universe_pdb.residues),
                                      n_segments=len(universe_pdb.segments),
                                      atom_resindex=[0] * len(universe.atoms),
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
    
    # Let's start trivial
    for at in new_universe.atoms:
        if at.index in match_list:
            # From PDB
            at.element = universe_pdb.atoms[match_list[at.index]].element
            # Overwrite the tinker names, that here are just harmful
            at.name = universe_pdb.atoms[match_list[at.index]].name
            match_residue = universe_pdb.atoms[match_list[at.index]].residue
            at.residue = new_universe.residues[match_residue.ix]
            # From tinker XYZ
            at.type = universe.atoms[at.index].type

    ## Then think to the unmatched atoms
    #for at in new_universe.atoms:
    #    if at.index not in match_list:
    #        # Only go ahead if the atom has just one bond!
    #        if len(at.bonds) != 1:
    #            print("""Atom {} does not match to any atom in PDB and """
    #                  """we do not understand why it is here (number of """
    #                  """bond different from 1).""")
    #            return None
    #            
    #        guessed_element = mda.topology.guessers.guess_atom_element(universe.atoms[at.index].name)
    #        # Check that it looks like an hydrogen...
    #        if guessed_element == 'H':
    #            # Guess that it is an hydrogen atom
    #            # It was probably added during input preparation so it is
    #            # allowed to be here.
    #            at.element = 'H'
    #            at.name = 'H'
    #            at.type = universe.atoms[at.index].type
    #            at.residue = at.bonds[0].partner(at).residue
    #        elif guessed_element == 'O':
    #            # Ok that's a carbonyl oxygen let's see if it is a terminal one.
    #            partn = at.bonds[0].partner(at)
    #            nox = 0
    #            nca = 0
    #            for b in partn.bonds:
    #                if b.partner(partn).element == 'O' \
    #                   or b.partner(partn) == at:
    #                    nox += 1
    #                elif b.partner(partn).element == 'C':
    #                    nca += 1
    #            if nox == 2 and nca == 1:
    #                # Ok it actually seems a terminal carbon the oxygen
    #                # is allowed to stay there
    #                at.element = 'O'
    #                at.name = 'O'
    #                at.type = universe.atoms[at.index].type
    #                at.residue = at.bonds[0].partner(at).residue
    #        else:
    #            print("""Atom {} does not match to any atom in PDB and """
    #                  """we do not understand why it is here (element """
    #                  """or connectivity unrecognized).""")
    #            return None
    return new_universe

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='OMMP-IP',
                        description="""OpenMMPol Input Preparation program, 
                                    this script should assist the process 
                                    creating input for OpenMMPol library 
                                    starting from standard structure file.
                                    It is mainly oriented to protein 
                                    complexes currently.""",
                        epilog='Good luck with your calculation!')

    parser.add_argument('-i', '--input-file',
                        required=True,
                        metavar='<INPUT_FILE>',
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
    parser.add_argument('-o', '--output-basename',
                        required=False,
                        metavar='<OUTPUT_BASENAME>',
                        dest='output_basename',
                        default='out',
                        help="""All output file will be created using this 
                                as basename. Default is 'out'""")

    args = parser.parse_args()
    qmmm_sys = input_load(args.input_file, args.pdb_input_file)
    
    if qmmm_sys is None:
        exit(1)

    print(qmmm_sys.select_atoms('all').resnums)
    print(qmmm_sys.select_atoms('backbone'))
