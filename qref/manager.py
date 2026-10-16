"""QRef as part of the geometry restraints of the whole model."""
from __future__ import absolute_import, division

import os
import sys

from cctbx.geometry_restraints.manager import manager as standard_manager

from qref import qref
from qref.common import Settings


def in_weight_calculation():
    """Whether the call comes from the calculation of the X-ray/geometry weight
    for coordinates in mmtbx/refinement/weights.py."""
    frame = sys._getframe(1)
    while frame is not None:
        if frame.f_code.co_filename.replace(os.sep, '/').endswith('mmtbx/refinement/weights.py'):
            return True
        frame = frame.f_back
    return False


def macro_cycle():
    """The integer macro_cycle of the first caller's self that has one, or None."""
    frame = sys._getframe(1)
    while frame is not None:
        number = getattr(frame.f_locals.get('self'), 'macro_cycle', None)
        if isinstance(number, int) and not isinstance(number, bool):
            return number
        frame = frame.f_back
    return None


class QRefManager(standard_manager):
    """The geometry restraints of the whole model, with QRef's terms added to
    every energy computed with gradients outside the weight calculation."""

    def energies_sites(self, sites_cart, *args, **kwargs):
        result = standard_manager.energies_sites(self, sites_cart, *args, **kwargs)
        # without normalization, target is residual_sum
        if (result.gradients is not None and not result.normalization
                and os.path.exists('qref.dat') and not in_weight_calculation()):
            result.residual_sum, result.gradients = qref.run(sites_cart=sites_cart,
                mm_residual_sum=result.residual_sum, mm_gradients=result.gradients,
                pdb_interpretation_params=self.pdb_interpretation_params,
                macro_cycle=macro_cycle())
            result.target = result.residual_sum
        return result


def digester(model, geometry, pdb_interpretation_params):
    """The geometry restraints manager for model: a QRefManager for the whole
    model when qref.dat is present, otherwise geometry itself."""
    if not os.path.exists('qref.dat') or model.get_number_of_atoms() != Settings.read().n_atoms:
        return geometry
    qref_manager = QRefManager.__new__(QRefManager)
    qref_manager.__dict__.update(vars(geometry))
    qref_manager.pdb_interpretation_params = pdb_interpretation_params
    return qref_manager
