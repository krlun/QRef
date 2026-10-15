import json
import os

import pytest

from conftest import REFERENCE, ROOT
from qref import qref
from qref.common import Region
import qref_prep

# at a link atom the step reaches Orca scaled by the junction factor g and
# rounded to the three decimals of qm_i_h.pdb, so one near 0.001 is mostly
# rounding error; a much larger one picks up the curvature of the stub surface
STEP = 0.05
# 0.0005/(g*STEP) bounds the rounding error at 0.014, above this tolerance; it
# passes only because g*STEP lands on the grid, so recheck if g or STEP changes
TOLERANCE = 0.01


def contribution(model, sites_cart=None):
    """QRef's contribution alone, with the MM part going in as zero.  run()
    adds into the gradient array in place, so each call needs a fresh one."""
    sites = model.get_sites_cart() if sites_cart is None else sites_cart
    target, gradients = qref.run(sites_cart=sites,
                                 mm_residual_sum=0.0,
                                 mm_gradients=model.get_sites_cart() * 0.0)
    return target, gradients


def moved(serial, axis, sites_cart, step):
    sites = sites_cart.deep_copy()
    position = list(sites[serial - 1])
    position[axis] += step
    sites[serial - 1] = tuple(position)
    return sites


def test_unchanged(prepared, installation, recording):
    """Target and gradients against tests/reference.json, to the last digit;
    --record writes the entry for the installation in use."""
    name, model = prepared
    target, gradients = contribution(model)
    # 17 significant digits round-trip a double, so comparing the strings compares
    # the values exactly
    result = {'target': format(float(target), '.17g'),
              'gradients': dict((str(i), [format(c, '.17g') for c in g])
                                for i, g in enumerate(gradients)
                                if g != (0.0, 0.0, 0.0))}

    reference = {}
    if os.path.exists(REFERENCE):
        with open(REFERENCE) as handle:
            reference = json.load(handle)

    if recording:
        reference.setdefault(installation, {})[name] = result
        with open(REFERENCE, 'w') as handle:
            json.dump(reference, handle, indent=2, sort_keys=True)
        pytest.skip('recorded %s for %s' % (name, installation))

    expected = reference.get(installation, {}).get(name)
    if expected is None:
        pytest.skip('no reference for %s under %s; run with --record'
                    % (name, installation))
    assert result['target'] == expected['target']
    assert sorted(result['gradients']) == sorted(expected['gradients'])
    for index in sorted(expected['gradients'], key=int):
        assert result['gradients'][index] == expected['gradients'][index], \
            'serial %d' % (int(index) + 1)


@pytest.mark.parametrize('serial, axis', [(4555, 0), (4555, 1), (4558, 2),
                                          (4600, 0), (4602, 1), (20943, 2)])
def test_gradient_matches_finite_difference(prepared, serial, axis):
    """4555 is a link atom and 4558 its bonded partner, so the QM gradient there is
    split by g instead of passed straight through.  With two regions, 4600 and
    4602 carry the angle restraint and 20943 the distance restraint."""
    name, model = prepared
    sites_cart = model.get_sites_cart()
    _, gradients = contribution(model, sites_cart)
    analytic = gradients[serial - 1][axis]

    plus, _ = contribution(model, moved(serial, axis, sites_cart, +STEP))
    minus, _ = contribution(model, moved(serial, axis, sites_cart, -STEP))
    numeric = (plus - minus) / (2.0 * STEP)

    assert analytic == pytest.approx(numeric, rel=TOLERANCE), (
        'atom %d axis %d: gradient %.8f, finite difference %.8f'
        % (serial, axis, analytic, numeric))


def test_restart_file_follows_the_coordinates(prepared):
    """run() rewrites the restart file from sites_cart, to three decimals."""
    name, model = prepared
    with open('qref.dat') as handle:
        restart = json.load(handle)['restart']
    if restart is None:
        pytest.skip('%s has no restart file' % name)

    serial, axis = 4555, 0
    sites_cart = moved(serial, axis, model.get_sites_cart(), +0.5)
    contribution(model, sites_cart)

    with open(restart) as handle:
        for line in handle:
            if line[0:6].strip() in ('ATOM', 'HETATM') \
                    and int(line[6:11]) == serial:
                written = [float(line[30:38]), float(line[38:46]),
                           float(line[46:54])]
                break
        else:
            raise AssertionError('serial %d not in %s' % (serial, restart))
    assert written == pytest.approx(list(sites_cart[serial - 1]), abs=5e-4)


def test_log_gains_a_row_per_call(prepared):
    """run() appends one row per call to qref_i.log, numbered from 1."""
    name, model = prepared
    contribution(model)
    contribution(model)

    with open('qref_1.log') as handle:
        header, first, second = handle.read().splitlines()
    assert header.split() == ['iter', 'QM', 'energy', '(kJ/mol)', 'MM', 'energy',
                              '("kJ/mol")', 'MM1', 'energy', '("kJ/mol")',
                              'QM/MM', 'energy', '("kJ/mol")']
    for row, iteration in ((first, '1'), (second, '2')):
        fields = row.split()
        assert fields[0] == iteration
        assert len(fields) == 5
        for energy in fields[1:]:
            float(energy)


@pytest.mark.parametrize('failure', ['exit', 'noscf'])
def test_failed_orca_falls_back_on_the_previous_gradient(prepared, monkeypatch,
                                                         failure):
    """Logged as Failed and refined on with the gradient Orca left last time, so
    an unconverged step does not end the refinement."""
    name, model = prepared
    _, expected = contribution(model)

    monkeypatch.setenv('STUB_ORCA_FAIL', failure)
    _, gradients = contribution(model)

    assert list(gradients) == list(expected)
    with open('qref_1.log') as handle:
        assert handle.read().splitlines()[-1].split() == ['2', 'Failed']


def test_failed_orca_raises_with_no_gradient_to_fall_back_on(prepared,
                                                             monkeypatch):
    """Nothing to refine on when the first call is the one that fails."""
    name, model = prepared
    monkeypatch.setenv('STUB_ORCA_FAIL', 'exit')
    with pytest.raises(qref.OrcaFailed):
        contribution(model)


def test_lock_is_released(prepared, monkeypatch):
    """A lock left behind silently disables QRef on the next refinement started
    in the same directory."""
    name, model = prepared
    monkeypatch.setenv('STUB_ORCA_FAIL', 'exit')
    with pytest.raises(qref.OrcaFailed):
        contribution(model)
    assert not os.path.exists('qm.lock')

    monkeypatch.delenv('STUB_ORCA_FAIL')
    contribution(model)
    assert not os.path.exists('qm.lock')


# syst1 files in order, and how qref_prep.py exits on them (None: it carries on);
# the second occurrence of a serial makes it a link atom
OVERLAPS = {
    'disjoint': (['1-5\n5\n', '6-9\n'], None),
    'link atom in all three': (['1-5\n5\n', '5-9\n5\n', '5 10-12\n5\n'], None),
    'link atom, then QM atom': (['1-5\n5\n', '5-9\n'],
                                'syst12 shares atoms 5 with syst11.'),
    'QM atom, then link atom': (['1-5\n', '5-9\n5\n'],
                                'syst12 shares atoms 5 with syst11.'),
    'QM atoms in the first and third': (['1-5\n5\n', '6-9\n', '3-4 10\n'],
                                        'syst13 shares atoms 3, 4 with syst11.'),
}


@pytest.mark.parametrize('name', sorted(OVERLAPS))
def test_regions_share_only_link_atoms(name, tmp_path, monkeypatch):
    """E_MM counts the restraints among atoms that QM regions share once, but each
    region subtracts them, so only link atoms, which have none among them, may be
    shared."""
    texts, message = OVERLAPS[name]
    monkeypatch.chdir(tmp_path)
    files = ['syst1%d' % index for index in range(1, len(texts) + 1)]
    for syst1, text in zip(files, texts):
        with open(syst1, 'w') as handle:
            handle.write(text)
    if message is None:
        qref_prep.check_overlap(files)
        return
    with pytest.raises(SystemExit) as raised:
        qref_prep.check_overlap(files)
    assert str(raised.value).startswith(message)


def test_link_atoms_do_not_count_as_altlocs(tmp_path, capsys):
    """A region of altloc A cut at blank link atoms is one conformation."""
    # MET A 102 in 5n6y: CA to CE altloc A, cut at the blank N and C
    syst1 = str(tmp_path / 'syst1')
    with open(syst1, 'w') as handle:
        handle.write('793 794 796-800\n793 794\n')
    region = Region(1, syst1)
    with open(os.path.join(ROOT, 'examples', '5n6y_CO3', '5n6y_sorted.pdb')) as handle:
        lines = [line for line in handle if line.startswith('CRYST1')
                 or (line[0:6].strip() in ('ATOM', 'HETATM')
                     and int(line[6:11]) in region.serials)]
    path = str(tmp_path / 'met102_a.pdb')
    with open(path, 'w') as handle:
        handle.writelines(lines)

    from iotbx.data_manager import DataManager
    manager = DataManager()
    manager.process_model_file(path)
    model = manager.get_model(filename=path)
    # as qref_prep.main() does before check_altlocs()
    region.restore_serials(model)

    qref_prep.check_altlocs(model, region.link_atoms)
    assert capsys.readouterr().out == 'Altloc:  "and altloc A"\n'
    qref_prep.check_altlocs(model, link_atoms=set())
    assert capsys.readouterr().out.startswith('Warning, multiple altlocs')
