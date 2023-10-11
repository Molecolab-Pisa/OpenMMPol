import json
import hashlib
import os
import sys

# 1. Import json file
# 2. Edit json file
# 2a. Remove qm section
# 2b. Remove linkatom section
# 3. Create hdf5 file
# 4. Create a new json file for hdf5
# 4a. Starting from original json, remove xyz_file and prm_file
#     or mmpol_file
# 4b. Add to the new json hdf5_file (and its md5sum)
# 4c. Save new json on file

json_input = sys.argv[1]
base_out = sys.argv[2]
ommp_pp = sys.argv[3]

# 1. Import json file
with open(json_input, 'r') as f:
    jin = json.loads(f.read())

# 2. Edit json file
jtmp = jin.copy()
# 2a. Remove qm section
jtmp.pop('qm', None)
# 2b. Remove linkatom section
jtmp.pop('link_atoms', None)
# 2c. Save tempjson on file
json_tmp = base_out+'_tmp.json'
with open(json_tmp, 'w') as f:
    print(json.dumps(jtmp, indent=4), file=f)
# 3. Create hdf5 file
hdf5_out = base_out+'.hdf5'
if os.path.exists(hdf5_out):
    os.remove(hdf5_out)
os.system(ommp_pp+' '+json_tmp+' '+hdf5_out)
os.remove(json_tmp)
# 4. Create a new json file for hdf5
jout = jin.copy()
# 4a. Starting from original json, remove xyz_file and prm_file
#     or mmpol_file
jout.pop('xyz_file', None)
jout.pop('prm_file', None)
jout.pop('mmpol_file', None)
# 4b. Add to the new json hdf5_file (and its md5sum)
jout['hdf5_file'] = {}
jout['hdf5_file']['path'] =  hdf5_out
with open(hdf5_out, 'rb') as f:
    hdf5_out_md5 = hashlib.md5(f.read()).hexdigest()
jout['hdf5_file']['md5sum'] = hdf5_out_md5

# 4c. Save new json on file
json_out = base_out+'.json'
with open(json_out, 'w') as f:
    print(json.dumps(jout, indent=4), file=f)

