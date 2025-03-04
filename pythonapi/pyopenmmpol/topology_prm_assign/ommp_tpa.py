import MDAnalysis as mda
import sys, os
from pyopenmmpol.topology_prm_assign import *
import argparse
import logging

def ommp_tpa_main():
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)

    allowed_cmd = ['assign', 'database']

    if len(sys.argv) < 2 or sys.argv[1] not in allowed_cmd:
        print("Use as first argument one of the following commands:")
        for c in allowed_cmd:
            print("\t{:s}".format(c))

    elif sys.argv[1] == 'assign':
        parser = argparse.ArgumentParser(prog='ommp-tpa assign',
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
                            default=None,
                            dest='out_txyz')
        parser.add_argument('--no-guess-elements',
                            action='store_true',
                            help='''Do not guess elements of atoms in input structure.''',
                            default=False,
                            dest='no_guess_elements')
        parser.add_argument('--no-guess-bonds',
                            action='store_true',
                            help='''Do not guess bonds of input structure.''',
                            default=False,
                            dest='no_guess_bonds')
        parser.add_argument('--no-overconnected-correction',
                            action='store_true',
                            help='''Do not correct bonds on overconnected atoms.''',
                            default=False,
                            dest='no_overconnected_corr')
        parser.add_argument('--no-overconnected-check',
                            action='store_true',
                            help='''Do not check bonds on overconnected atoms.''',
                            default=False,
                            dest='no_overconnected_check')

        args = parser.parse_args(args=sys.argv[2:])
        
        universe = mda.Universe(args.molf)
        
        if not args.no_guess_elements:
            guess_elements(universe)
        if not args.no_guess_bonds:
            guess_bonds(universe)
        if not args.no_overconnected_corr and not args.no_overconnected_check:
            check_for_overconnected_atoms(universe, True)
        elif not args.no_overconnected_check:
            check_for_overconnected_atoms(universe, False)

        mya = PrmAssignament(universe)
        logger.info("Number of molecule in the system {:d}".format(mya.nmol))

        logger.info("Loading db from file")
        mydb = AssignamentDB(args.db)
        mya.set_db(mydb)

        mya.topology_assign()

        logger.info("Assigned {:d} / {:d}".format(mya.tot_assigned, mya.natoms))
        logger.info(mya.assignament_log())

        if hasattr(args, 'out_txyz') and args.out_txyz is not None:
            mya.save_tinker_xyz(args.out_txyz)

    elif sys.argv[1] == 'database':
        parser = argparse.ArgumentParser(prog='ommp-tpa database',
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
                            fragment_name smile_string    enviroment_atoms bridge_atoms 3-letter_code residue_type
                            water         O               -                -             WAT           solvent
                            gly_intern    O=C([N])CN[C]=O 2,5,6            2,5           GLY           protein_res''',
                            metavar='<models.txt>',
                            dest='models')
        
        parser.add_argument('--plot-models',
                            action='store_true',
                            help='''Plot graph of models added to database.''',
                            default=False,
                            dest='plot_models')
        
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

        parser.add_argument('--learn-only',
                            choices=['name', 'type'],
                            help='''Take only atom name/type from the input structural file''',
                            dest='learn_only')
        
        parser.add_argument('--remove-fragment-without-types',
                            action='store_true',
                            help='''All parameter in database without atomtypes are just removed.''',
                            default=False,
                            dest='remove_notypes')

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

        if hasattr(args, 'models') and args.models is not None:
            for modf in args.models:
                with open(modf, 'r') as f:
                    for l in f:
                        if not l.strip():
                            continue
                        tok = l.split()
                        name = tok[0]
                        smiles_str = tok[1]
                        if tok[2] == '-':
                            env_atm = []
                        else:
                            env_atm = [int(i) for i in tok[2].split(',')]
                        
                        if tok[3] == '-':
                            bridge_atm = []
                        else:
                            bridge_atm = [int(i) for i in tok[3].split(',')]

                        threelc = tok[4]
                        restype = tok[5]
                        
                        f = Fragment(name,
                                    smiles_str,
                                    env_atm,
                                    bridge_atm,
                                    default_resname=threelc,
                                    restype=restype)
                        if args.plot_models:
                            f.draw()
                        db.add_fragment(f)
                        db_modified = True
        
        if hasattr(args, 'learn_struct') and args.learn_struct is not None:
            learn_name = True
            learn_type = True
            if hasattr(args, 'learn_only'):
                if args.learn_only == 'name':
                    learn_type = False
                elif args.learn_only == 'type':
                    learn_name = False

            for molsf in args.learn_struct:
                logger.info("Learning from {:s}".format(molsf))
                
                universe = mda.Universe(molsf)
                guess_elements(universe)
                guess_bonds(universe)
                check_for_overconnected_atoms(universe, True)
                
                mya = PrmAssignament(universe)
                mya.set_db(db)
                mya.topology_assign()
                mya.learn_from(use_name=learn_name, 
                            use_type=learn_type)
                db_modified = True

        if args.remove_notypes:
            logger.info("Removing fragments without atom-types")
            frag_to_remove = []
            for frag in db:
                if not frag.has_atomtypes():
                    frag_to_remove += [frag]
            for frag in frag_to_remove:
                db.remove_fragment(frag)
                db_modified = True


        if db_modified or not db_exists:
            logger.info("Writing database")
            db.save_as_json(args.db)

