"""Tests for QRefManager and digester(), on the one-region case of the
7s4h_CuD example."""
import os

import pytest

from cctbx.geometry_restraints.manager import manager as standard_manager
from scitbx.array_family import flex

from qref import manager
from qref import qref

ONE_REGION = pytest.mark.parametrize('prepared', ['one region'], indirect=True)

WEIGHT_CALCULATION = '''
def energies(geometry, sites_cart):
    return geometry.energies_sites(sites_cart=sites_cart, compute_gradients=True)
'''


class MacroCycle(object):
    """A caller in macro-cycle macro_cycle."""

    def __init__(self, macro_cycle):
        self.macro_cycle = macro_cycle

    def energies(self, geometry, sites_cart):
        return geometry.energies_sites(sites_cart=sites_cart, compute_gradients=True)


@pytest.fixture
def geometry(prepared, monkeypatch):
    """The model, its geometry restraints manager and the parameters it was made
    with, processed in a directory without qref.dat."""
    name, model = prepared
    work = os.getcwd()
    os.mkdir('standard')
    monkeypatch.chdir('standard')
    params = model.get_default_pdb_interpretation_params()
    model.process(pdb_interpretation_params=params, make_restraints=True)
    monkeypatch.chdir(work)
    geometry = model.get_restraints_manager().geometry
    assert type(geometry) is standard_manager
    return model, geometry, params


@pytest.fixture
def no_qref(monkeypatch):
    """Make qref.run() fail the test."""
    def run(**kwargs):
        raise AssertionError('qref.run() was called')
    monkeypatch.setattr(qref, 'run', run)


def assert_same(result, expected):
    assert result.residual_sum == expected.residual_sum
    assert result.target == expected.target
    if expected.gradients is None:
        assert result.gradients is None
    else:
        assert list(result.gradients) == list(expected.gradients)


@ONE_REGION
def test_digester_needs_qref_dat(geometry):
    model, standard, params = geometry
    os.remove('qref.dat')
    assert manager.digester(model, standard, params) is standard


@ONE_REGION
def test_digester_needs_the_whole_model(geometry):
    model, standard, params = geometry
    part = model.select(model.selection('not water'))
    assert manager.digester(part, standard, params) is standard


@ONE_REGION
def test_energies_with_gradients_add_qref(geometry):
    """The energies of a QRefManager are those of qref.run() on the standard
    energies."""
    model, standard, params = geometry
    qref_manager = manager.digester(model, standard, params)
    assert isinstance(qref_manager, manager.QRefManager)
    sites_cart = model.get_sites_cart()

    result = qref_manager.energies_sites(sites_cart=sites_cart, compute_gradients=True)
    mm = standard.energies_sites(sites_cart=sites_cart, compute_gradients=True)
    target, gradients = qref.run(sites_cart=sites_cart, mm_residual_sum=mm.residual_sum,
                                 mm_gradients=mm.gradients, pdb_interpretation_params=params)
    assert target != mm.residual_sum
    assert result.residual_sum == result.target == target
    assert list(result.gradients) == list(gradients)


@ONE_REGION
@pytest.mark.parametrize('kwargs', [{'compute_gradients': False},
                                    {'compute_gradients': True, 'normalization': True}],
                         ids=['without gradients', 'normalized'])
def test_other_energies_are_standard(geometry, no_qref, kwargs):
    model, standard, params = geometry
    qref_manager = manager.digester(model, standard, params)
    sites_cart = model.get_sites_cart()
    assert_same(qref_manager.energies_sites(sites_cart=sites_cart, **kwargs),
                standard.energies_sites(sites_cart=sites_cart, **kwargs))


@ONE_REGION
def test_weight_calculation_is_standard(geometry, no_qref):
    model, standard, params = geometry
    qref_manager = manager.digester(model, standard, params)
    sites_cart = model.get_sites_cart()
    namespace = {}
    exec(compile(WEIGHT_CALCULATION, os.path.join('mmtbx', 'refinement', 'weights.py'),
                 'exec'), namespace)
    assert_same(namespace['energies'](qref_manager, sites_cart),
                standard.energies_sites(sites_cart=sites_cart, compute_gradients=True))


@ONE_REGION
def test_selection_is_standard(geometry, no_qref):
    model, standard, params = geometry
    qref_manager = manager.digester(model, standard, params)
    selected = qref_manager.select(selection=flex.bool(model.get_number_of_atoms(), True))
    assert type(selected) is standard_manager
    sites_cart = model.get_sites_cart()
    assert_same(selected.energies_sites(sites_cart=sites_cart, compute_gradients=True),
                standard.energies_sites(sites_cart=sites_cart, compute_gradients=True))


@ONE_REGION
def test_macro_cycle_of_the_caller_is_passed_on(geometry, monkeypatch):
    model, standard, params = geometry
    qref_manager = manager.digester(model, standard, params)
    passed = []
    def run(**kwargs):
        passed.append(kwargs['macro_cycle'])
        return kwargs['mm_residual_sum'], kwargs['mm_gradients']
    monkeypatch.setattr(qref, 'run', run)
    sites_cart = model.get_sites_cart()
    qref_manager.energies_sites(sites_cart=sites_cart, compute_gradients=True)
    MacroCycle(3).energies(qref_manager, sites_cart)
    assert passed == [None, 3]
