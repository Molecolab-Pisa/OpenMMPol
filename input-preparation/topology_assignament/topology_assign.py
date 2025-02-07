import numpy as np
import MDAnalysis as mda
import json
import sys
from topas import *

mya = PrmAssignament(mda.Universe(sys.argv[1]))

print("Loading db from file...")
with open(sys.argv[2], 'r') as f:
    json_in = json.load(f)
db = [dict_to_frag(f) for f in json_in['res_assignament']]
db = sorted(db, key=lambda frag: frag.priority, reverse=True)
mya.set_db(db)

mya.topology_assign()

print("Assigned {:d} / {:d}".format(mya.tot_assigned, mya.natoms))

mya.learn_from()

