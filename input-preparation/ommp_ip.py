import numpy as np
import MDAnalysis as mda
import json
import argparse
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

def match_atoms(uni1, uni2):
    """This function takes two universes and return a list
    of atoms touple [(at1@uni1, at1@uni2), ...] containing 
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

    matchlist = []
    for i, at in enumerate(uniS.atoms):
        mindist = np.min(distmatrix[i,:])
        if mindist < DIST_THR:
            # Ok this atom has a match!
            match_at = np.argmin(distmatrix[i,:])
            if uniL == uni1:
                matchlist += [(uniL.atoms[match_at].id, at.id)]
            else:
                matchlist += [(at.id, uniL.atoms[match_at].id)]
        else:
            print("No matching atom found for atom {}".format(at))
    
    return matchlist

def input_load(fin_txyz, fin_pdb=None):
    try:
        universe = mda.Universe(fin_txyz, format='txyz')
    except FileNotFoundError:
        print("{:s} not found, cannot open input.".format(fin_txyz))
        universe = None
    except (ValueError, IndexError):
        print("""{:s} doesen't seems to be a correctly"""
              """formatted Tinker xyz.""".format(fin_txyz))
        universe = None

    if fin_pdb is not None:
        try:
            universe_pdb = mda.Universe(fin_pdb, format='pdb')
        except FileNotFoundError:
            print("{:s} not found, cannot open input.".format(fin_pdb))
            universe_pdb = None
        except (ValueError, IndexError):
            print("""{:s} doesen't seems to be a correctly"""
                  """formatted PDB.""".format(fin_pdb))
            universe_pdb = None

        if universe_pdb is None:
            # To signal an error in the calling function
            universe = None

    print(match_atoms(universe, universe_pdb))

    return universe

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
