"""Fixtures that set up the 7s4h_CuD example in a temporary directory, ready
for qref.run()."""
import json
import os
import pickle
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

CASES = {
    'one region': {
        'syst1': {'syst1': ONE},
        'options': [],
    },
    'two regions with a transform and restraints': {
        'syst1': {'syst11': TWO_A, 'syst12': TWO_B},
        'options': (['-t'] + TRANSFORM
                    + ['-rd', '1', '20943', '4561', '2.1', '2500']
                    + ['-ra', '2', '4598', '4600', '4602', '109.5', '10']),
    },
}


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
        return version
    import libtbx.load_env
    path = abs(libtbx.env.build_path)
    while path not in (os.sep, ''):
        path, name = os.path.split(path)
        if name.startswith('phenix'):
            return name
    return 'unknown'


@pytest.fixture(params=sorted(CASES), ids=sorted(CASES))
def prepared(request, tmp_path, monkeypatch):
    """The example taken through qref_prep.py, once per region split.  Returns
    the case name and the full model."""
    name = request.param
    case = CASES[name]
    work = str(tmp_path)
    for copied in (MODEL, 'junctfactor') + CIFS:
        shutil.copy(os.path.join(EXAMPLE, copied), work)
    regions = sorted(case['syst1'])
    for index, region in enumerate(regions, 1):
        with open(os.path.join(work, region), 'w') as handle:
            handle.write(case['syst1'][region])
        with open(os.path.join(work, f'qm_{index}.inp'), 'w') as handle:
            handle.write('! TPSS D4 DEF2-SV(P)\n! ENGRAD\n'
                         f'*pdbfile 0 1 qm_{index}_h.pdb\n')
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
