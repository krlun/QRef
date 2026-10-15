"""Support code for the scripts.

Everything the scripts share with the QRef interface lives in qref/common.py
and is re-exported here, so that there is one definition of it. The qref module
is found either because it is installed in the installation of Phenix being
used, or, when a script is run from the QRef directory, beside the script.
"""
from __future__ import absolute_import, division, print_function

import os
import re
import sys

try:
    import qref
except ImportError:
    root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    if not os.path.isfile(os.path.join(root, 'qref', '__init__.py')):
        raise SystemExit('Cannot find the qref module. Install QRef with '
                         'install.py, or run the script from the QRef '
                         'directory. Exiting...')
    sys.path.insert(0, root)

from qref.common import Region
from qref.common import read_settings
from qref.common import read_region


def read_junctfactor(junc_factor_file):
    # a line 'residue atom1 atom2 type' starts a bond, followed by lines
    # 'ltype distance' up to a blank line; distances are written as Fortran
    # doubles, e.g. 1.1070d0
    junctfactor = dict()
    distances = None
    with open(junc_factor_file, 'r') as handle:
        for line in handle:
            fields = re.split('[#!]', line)[0].split()
            if not line.strip():
                distances = None
            elif len(fields) == 4:
                resname, bond = fields[0], '-'.join(fields[1:3])
                distances = junctfactor.setdefault(resname, dict()).setdefault(bond, dict())
            elif distances is not None and fields:
                distances[int(fields[0])] = float(fields[1].split('d')[0])
    return junctfactor
