"""Fixtures that set up the 7s4h_CuD example in a temporary directory, ready
for qref.run()."""
import json
import os
import pickle
import re
import shutil
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.join(ROOT, 'tests')
# ahead of the qref installed under Phenix, so the tests exercise the working copy
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

EXAMPLE = os.path.join(ROOT, 'examples', '7s4h_CuD')
MODEL = '7s4h_sorted.pdb'
CIFS = ('HXG.cif', 'PLC.cif')
STUB = os.path.join(HERE, 'stub_orca.py')
REFERENCE = os.path.join(HERE, 'reference.json')

# install.py is Python 3 only
collect_ignore = [] if sys.version_info[0] >= 3 else ['test_install.py']

# the second occurrence of a serial makes it a link atom
ONE = """4555 ! ASN C 227
4558-4565 ! ASN C 227
4597-4605 ! HIS C 231
4725-4733 ! HIS C 245
20943 ! CU C
25036-25038 ! HOH C 406
25047-25049 ! HOH C 415

4555 ! CA ASN C 227
4597 ! CB HIS C 231
4725 ! CB HIS C 245
"""

# the bonded partner of a link atom must be in the same region, or qref_prep.py
# exits: 4555 with 4558, 4597 with 4598, 4725 with 4726
TWO_A = """4555 ! ASN C 227
4558-4565 ! ASN C 227
20943 ! CU C
25036-25038 ! HOH C 406
25047-25049 ! HOH C 415

4555 ! CA ASN C 227
"""

TWO_B = """4597-4605 ! HIS C 231
4725-4733 ! HIS C 245

4597 ! CB HIS C 231
4725 ! CB HIS C 245
"""

# region, atoms, R row-wise, t: a quarter turn about z and 10 Angstrom along x,
# so run() has to rotate the gradients back
TRANSFORM = ['2', '4725-4733',
             '0', '-1', '0',
             '1', '0', '0',
             '0', '0', '1',
             '10.0', '0.0', '0.0']

# applied after TRANSFORM to the same atoms: a quarter turn about x and 5 Angstrom
# along y, which does not commute with TRANSFORM
SECOND_TRANSFORM = ['2', '4725-4733',
                    '1', '0', '0',
                    '0', '0', '-1',
                    '0', '1', '0',
                    '0.0', '5.0', '0.0']

# atoms put ahead of the example, so that each of its serials goes past 99999 and
# is written in hybrid-36
PADDING = 100000


def shifted(syst1, offset):
    """syst1 with each serial raised by offset, and the comments dropped."""
    return ''.join(re.sub(r'\d+', lambda match: str(int(match.group()) + offset),
                          re.split('[#!]', line)[0]) + '\n'
                   for line in syst1.splitlines())


def write_padded(source, destination, count):
    """source behind a chain of count water oxygens, well away from it, so that
    each of its atoms is count further along.  The padded model is put in the box
    source gets, which the padding lies outside of; a box of its own would change
    the last digits of the MM1 energies."""
    from iotbx.data_manager import DataManager
    from iotbx.pdb import hierarchy
    from iotbx.pdb.hybrid_36 import hy36encode
    manager = DataManager()
    manager.process_model_file(source)
    original = manager.get_model(filename=source)
    original.add_crystal_symmetry_if_necessary()
    root = hierarchy.root()
    model = hierarchy.model()
    root.append_model(model)
    chain = hierarchy.chain(id='ZZ')
    model.append_chain(chain)
    for i in range(count):
        residue_group = hierarchy.residue_group(resseq=hy36encode(4, i + 1))
        atom_group = hierarchy.atom_group(resname='HOH')
        atom = hierarchy.atom()
        atom.name = ' O  '
        atom.element = ' O'
        atom.xyz = (1000 + 3*(i % 50), 1000 + 3*(i//50 % 50), 1000 + 3*(i//2500))
        atom.occ = 1.0
        atom.b = 30.0
        atom_group.append_atom(atom)
        residue_group.append_atom_group(atom_group)
        chain.append_residue_group(residue_group)
    for original_chain in original.get_hierarchy().only_model().chains():
        model.append_chain(original_chain.detached_copy())
    root.atoms().reset_serial()
    root.write_pdb_file(file_name=destination,
                        crystal_symmetry=original.crystal_symmetry())


CASES = {
    'one region': {
        'syst1': {'syst1': ONE},
        'options': [],
    },
    'two regions with a transform and restraints': {
        'syst1': {'syst11': TWO_A, 'syst12': TWO_B},
        'options': (['-t'] + TRANSFORM
                    + ['-rd', '1', '20943', '4561', '2.1', '2500']
                    + ['-ra', '2', '4598', '4600', '4602', '109.5', '10']
                    + ['-r']),
    },
    'two transforms on the same atoms': {
        'syst1': {'syst11': TWO_A, 'syst12': TWO_B},
        'options': ['-t'] + TRANSFORM + ['-t'] + SECOND_TRANSFORM,
    },
    'one region past 99999': {
        'syst1': {'syst1': shifted(ONE, PADDING)},
        'options': ['-r'],
        'padding': PADDING,
    },
}

# the cases every test taking `prepared` runs on; a test asks for any other case
# by name, with indirect parametrization
DEFAULT_CASES = ('one region', 'two regions with a transform and restraints')


def pytest_addoption(parser):
    parser.addoption('--record', action='store_true',
                     help='write the reference values instead of checking them')


@pytest.fixture(scope='session')
def recording(request):
    return request.config.getoption('--record')


@pytest.fixture(scope='session')
def installation():
    """Name of the Phenix installation in use.  The MM1 gradients come from
    cctbx, so a reference value holds for one installation only."""
    version = os.environ.get('PHENIX_VERSION')
    if version:
        # phenix_env.sh sets the bare version, the build path carries the prefix
        return version if version.startswith('phenix') else 'phenix-' + version
    import libtbx.load_env
    path = abs(libtbx.env.build_path)
    while path not in (os.sep, ''):
        path, name = os.path.split(path)
        if name.startswith('phenix'):
            return name
    return 'unknown'


@pytest.fixture(params=DEFAULT_CASES, ids=DEFAULT_CASES)
def prepared(request, tmp_path, monkeypatch):
    """The example taken through qref_prep.py, once per region split.  Returns
    the case name and the full model."""
    name = request.param
    case = CASES[name]
    work = str(tmp_path)
    for copied in (MODEL, 'junctfactor') + CIFS:
        shutil.copy(os.path.join(EXAMPLE, copied), work)
    if 'padding' in case:
        write_padded(os.path.join(EXAMPLE, MODEL), os.path.join(work, MODEL),
                     case['padding'])
    regions = sorted(case['syst1'])
    for index, region in enumerate(regions, 1):
        with open(os.path.join(work, region), 'w') as handle:
            handle.write(case['syst1'][region])
        with open(os.path.join(work, 'qm_%d.inp' % index), 'w') as handle:
            handle.write('! TPSS D4 DEF2-SV(P)\n! ENGRAD\n'
                         '*pdbfile 0 1 qm_%d_h.pdb\n' % index)
    monkeypatch.chdir(work)

    import qref_prep
    monkeypatch.setattr(sys, 'argv',
                        ['qref_prep.py', MODEL, '-c'] + list(CIFS)
                        + ['-s'] + regions + case['options'])
    qref_prep.main()

    from iotbx.data_manager import DataManager
    manager = DataManager()
    for cif in CIFS:
        manager.process_restraint_file(cif)
    manager.process_model_file(MODEL)
    model = manager.get_model(filename=MODEL)
    model.add_crystal_symmetry_if_necessary()
    # run() reads settings.pickle, which the patched mmtbx/model/model.py writes
    # during refinement and qref_prep.main() deletes
    with open('settings.pickle', 'wb') as handle:
        pickle.dump(model.get_default_pdb_interpretation_params(), handle)

    with open('qref.dat') as handle:
        settings = json.load(handle)
    settings['orca_binary'] = STUB
    with open('qref.dat', 'w') as handle:
        json.dump(settings, handle, indent=4, sort_keys=True)
    return name, model
