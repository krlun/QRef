import json
import os

import pytest

from conftest import REFERENCE
from qref import qref

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
    gradients, target = qref.run(sites_cart=sites,
                                 mm_gradients=model.get_sites_cart() * 0.0,
                                 mm_residual_sum=0.0)
    return gradients, target


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
    gradients, target = contribution(model)
    # 17 significant digits round-trip a double, so comparing the strings compares
    # the values exactly
    result = {'target': f'{float(target):.17g}',
              'gradients': dict((str(i), [f'{c:.17g}' for c in g])
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
        pytest.skip(f'recorded {name} for {installation}')

    expected = reference.get(installation, {}).get(name)
    if expected is None:
        pytest.skip(f'no reference for {name} under {installation}; '
                    f'run with --record')
    assert result['target'] == expected['target']
    assert sorted(result['gradients']) == sorted(expected['gradients'])
    for serial in sorted(expected['gradients'], key=int):
        assert result['gradients'][serial] == expected['gradients'][serial], serial


@pytest.mark.parametrize('serial, axis', [(4555, 0), (4555, 1), (4558, 2)])
def test_gradient_matches_finite_difference(prepared, serial, axis):
    """4555 is a link atom and 4558 its bonded partner, so the QM gradient there is
    split by g instead of passed straight through."""
    name, model = prepared
    sites_cart = model.get_sites_cart()
    gradients, _ = contribution(model, sites_cart)
    analytic = gradients[serial - 1][axis]

    _, plus = contribution(model, moved(serial, axis, sites_cart, +STEP))
    _, minus = contribution(model, moved(serial, axis, sites_cart, -STEP))
    numeric = (plus - minus) / (2.0 * STEP)

    assert analytic == pytest.approx(numeric, rel=TOLERANCE), (
        f'atom {serial} axis {axis}: gradient {analytic:.8f}, '
        f'finite difference {numeric:.8f}')


def test_restart_file_follows_the_coordinates(prepared):
    """run() rewrites the restart file from sites_cart, to three decimals."""
    name, model = prepared
    with open('qref.dat') as handle:
        restart = json.load(handle)['restart']
    if restart is None:
        pytest.skip(f'{name} has no restart file')

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
            raise AssertionError(f'serial {serial} not in {restart}')
    assert written == pytest.approx(list(sites_cart[serial - 1]), abs=5e-4)
