import MDAnalysis as mda
import sys, os
from topas import *
import argparse
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

allowed_cmd = ['assign', 'database']

if len(sys.argv) < 2 or sys.argv[1] not in allowed_cmd:
    print("Use as first argument one of the following commands:")
    for c in allowed_cmd:
        print("\t{:s}".format(c))

elif sys.argv[1] == 'assign':
    parser = argparse.ArgumentParser(prog='assign',
                                     description="""Perform atom type assignament on a file containing a molecule structure""")
    parser.add_argument('-d', '--db',
                        required=True,
                        type=str,
                        help='Database to be used for the assignament',
                        metavar='<db_name.json>',
                        dest='db')
    parser.add_argument('-f', '--mol-file',
                        required=True,
                        type=str,
                        help='''Molecular structure to be assigned.
                                Any kind of structural file that can be read by MDAnalysis is accepted.''',
                        metavar='<mol_file.xyz>',
                        dest='molf')
    
    parser.add_argument('--out-txyz-file',
                        type=str,
                        help='''Output Tinker xyz file containing the assignament information.''',
                        metavar='<out_file.xyz>',
                        dest='out_txyz')
    args = parser.parse_args(args=sys.argv[2:])

    mya = PrmAssignament(mda.Universe(args.molf))
    print("Number of molecule in the system {:d}".format(mya.nmol))

    print("Loading db from file...")
    mydb = AssignamentDB(args.db)
    mya.set_db(mydb)

    mya.topology_assign()

    print("Assigned {:d} / {:d}".format(mya.tot_assigned, mya.natoms))

    if hasattr(args, 'out_txyz'):
        mya.save_tinker_xyz(args.out_txyz)

elif sys.argv[1] == 'database':
    parser = argparse.ArgumentParser(prog='database',
                                     description="""Perform operations on database""",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--db',
                        required=True,
                        type=str,
                        help='Database to create or edit',
                        metavar='<db_name.json>',
                        dest='db')
    
    parser.add_argument('--add-models',
                        type=str,
                        nargs='+',
                        help='''File(s) containing a list of model to be added to the database.
                        It should be a normal text file with the following space
                        separated fields:
                        fragment_name smile_string    enviroment_atoms 3-letter_code residue_type
                        water         O               -                WAT           solvent
                        gly_intern    O=C([N])CN[C]=O 2,5,6            GLY           protein_res''',
                        metavar='<models.txt>',
                        dest='models')
    
    parser.add_argument('--learn-from',
                        type=str,
                        nargs='+',
                        help='''File(s) containing structures with assigned atom types/atom names.
                        Each file undergoes a topological assignament using current database, then
                        if atom types or atom names are not present in the db, they are extracted
                        from the file and inserted in the db. If the informations are present in 
                        both the db and the structure file, they are compared and if inconsistencies
                        are found an error is returned. All files accepted by MDAnalysis are taken
                        as valid. Note that the default extension for Tinker xyz is .arc.''',
                        metavar='<struct.mol2>',
                        dest='learn_struct')



    args = parser.parse_args(args=sys.argv[2:])

    if os.path.exists(args.db):
        logger.info("Database exists, it will be loaded from file.")
        db_exists = True
        db = AssignamentDB(args.db)
    else:
        logger.info("Database does not exists, it will be created.")
        db_exists = False
        db = AssignamentDB()

    db_modified = False

    if hasattr(args, 'models'):
        for modf in args.models:
            with open(modf, 'r') as f:
                for l in f:
                    tok = l.split()
                    name = tok[0]
                    smiles_str = tok[1]
                    if tok[2] == '-':
                        env_atm = []
                    else:
                        env_atm = [int(i) for i in tok[2].split(',')]
                    threelc = tok[3]
                    restype = tok[4]
                    
                    db.add_fragment(Fragment(name,
                                            smiles_str,
                                            env_atm,
                                            default_resname=threelc,
                                            restype=restype))
                    db_modified = True
    
    if hasattr(args, 'learn_struct'):
        for molsf in args.learn_struct:
            logger.info("Learning from {:s}".format(molsf))
            mya = PrmAssignament(mda.Universe(molsf))
            mya.set_db(db)
            mya.topology_assign()
            mya.learn_from()
            db_modified = True

    if db_modified or not db_exists:
        db.save_as_json(args.db)

#mya.learn_from()
#mydb.save_as_json(sys.argv[2])

