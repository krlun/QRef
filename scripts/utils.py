"""Support code for the scripts.

Everything the scripts share with the QRef interface lives in qref/common.py
and is re-exported here, so that there is one definition of it. The qref module
is found either because it is installed in the installation of Phenix being
used, or, when a script is run from the QRef directory, beside the script.
"""
from __future__ import absolute_import, division, print_function

import os
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
from qref.common import Settings
from qref.common import read_region
