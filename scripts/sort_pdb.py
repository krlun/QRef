#!/Applications/phenix-1.20.1-4487/build/bin/cctbx.python

import sys, os
import mmtbx.model
import iotbx.pdb

def main(args):
    if (len(args) != 1):
        raise RuntimeError("Please specify one pdb file name.")
    model_file = args[0]
    sorted_pdb_out = os.path.splitext(model_file)[0] + "_sorted.pdb"
    pdb = iotbx.pdb.input(file_name=model_file)
    model = mmtbx.model.manager(model_input=pdb)

    # print('Writing file: ', sorted_pdb_out) # Only works in Python 3
    print('Writing file:  ' + sorted_pdb_out) # Python 2 (cctbx.python)
    model.get_hierarchy().write_pdb_file(file_name=sorted_pdb_out, crystal_symmetry=pdb.crystal_symmetry())


if (__name__ == "__main__"):
    main(sys.argv[1:])
