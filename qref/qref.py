from __future__ import absolute_import, division
import os
import subprocess
import pickle

import numpy as np

from iotbx.data_manager import DataManager

from qref.common import Region
from qref.common import parse_serials
from qref.common import read_settings

KCAL_PER_HARTREE = 627.509474063112
KJ_PER_HARTREE = 2625.499639479950
ANGSTROM_PER_BOHR = 0.529177249


class OrcaFailed(Exception):
    """Orca failed with no earlier gradient to carry on with."""


def read_engrad(infile):
    gradients = list()
    with open(infile, 'r') as handle:
        line = handle.readline()
        while line:
            if line == '# Number of atoms\n':
                handle.readline()
                n_atoms = int(handle.readline().strip())
            if line == '# The current total energy in Eh\n':
                handle.readline()
                energy = float(handle.readline().strip())
                break
            line = handle.readline()
        for i in range(3): handle.readline()
        for i in range(n_atoms):
            gradient = (float(handle.readline()), float(handle.readline()), float(handle.readline()))
            gradients.append(gradient)
    return energy, np.array(gradients)


def rescale_qm_gradients(qm_gradients, w):
    return w*qm_gradients


def write_coordinates(infile, sites_cart):
    with open(infile, 'r') as handle:
        model = handle.readlines()
    records = {'ATOM', 'HETATM'}
    width = 8
    with open(infile, 'w') as handle:
        for line in model:
            if line[0:6].strip() in records:
                serial = int(line[6:11].strip())
                coords = ''
                for i in range(3): coords += '{:.3f}'.format(round(sites_cart[serial-1][i], 3)).rjust(width)
                line = line[0:30] + coords + line[54:]
            handle.write(line)


class DistanceRestraint(object):
    """Harmonic restraint force_constant*(r - distance)**2 on the distance r between
    atom1 and atom2, in Angstrom."""

    def __init__(self, atom1, atom2, distance, force_constant):
        self.atom1 = atom1
        self.atom2 = atom2
        self.distance = distance
        self.force_constant = force_constant

    def add_to(self, sites_cart, gradients):
        """Add the gradient of the restraint to gradients, and return its energy."""
        r_ij = np.array(sites_cart[self.atom1-1]) - np.array(sites_cart[self.atom2-1])
        r = np.linalg.norm(r_ij)
        delta_r = r - self.distance
        d_U_ij_d_r = 2*self.force_constant*delta_r
        d_r_d_r_ij = r_ij/r
        gradients[self.atom1-1] += d_U_ij_d_r*d_r_d_r_ij
        gradients[self.atom2-1] += -d_U_ij_d_r*d_r_d_r_ij
        return self.force_constant*delta_r**2


class AngleRestraint(object):
    """Harmonic restraint force_constant*(alpha - angle)**2 on the angle alpha
    atom1-atom2-atom3, in degrees, with atom2 at the tip."""

    def __init__(self, atom1, atom2, atom3, angle, force_constant):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.angle = angle
        self.force_constant = force_constant

    def add_to(self, sites_cart, gradients):
        """Add the gradient of the restraint to gradients, and return its energy.
        i, j and k are atom1, atom2 and atom3."""
        site2 = np.array(sites_cart[self.atom2-1])
        r_ij = np.array(sites_cart[self.atom1-1]) - site2
        r_kj = np.array(sites_cart[self.atom3-1]) - site2
        norm_r_ij = np.linalg.norm(r_ij)
        norm_r_kj = np.linalg.norm(r_kj)
        u_ij = r_ij/norm_r_ij
        u_kj = r_kj/norm_r_kj
        cos_alpha = np.dot(r_ij, r_kj)/(norm_r_ij * norm_r_kj)
        alpha = np.degrees(np.arccos(cos_alpha))
        delta_alpha = alpha - self.angle
        d_U_d_alpha = 2*self.force_constant*delta_alpha
        inv_sin_alpha = 1.0/np.sqrt(1 - cos_alpha**2)
        d_alpha_d_r_i = np.degrees(inv_sin_alpha * (1.0/norm_r_ij) * (cos_alpha * u_ij - u_kj))
        d_alpha_d_r_k = np.degrees(inv_sin_alpha * (1.0/norm_r_kj) * (cos_alpha * u_kj - u_ij))
        d_alpha_d_r_j = - d_alpha_d_r_i - d_alpha_d_r_k
        gradients[self.atom1-1] += d_U_d_alpha*d_alpha_d_r_i
        gradients[self.atom2-1] += d_U_d_alpha*d_alpha_d_r_j
        gradients[self.atom3-1] += d_U_d_alpha*d_alpha_d_r_k
        return self.force_constant*delta_alpha**2


class QMRegion(Region):
    """A Region during refinement: its restraints, the MM1 and QM energies and
    gradients that replace each other, and its log."""

    def __init__(self, number, syst1, settings):
        super(QMRegion, self).__init__(number, syst1, settings)
        self.restraints = ([DistanceRestraint(*restraint) for restraint in self.restraints_distance]
                           + [AngleRestraint(*restraint) for restraint in self.restraints_angle])
        self.cif = settings['cif']
        self.orca_binary = settings['orca_binary']
        self.w_qm = settings['w_qm']

    def _add_qm_mm(self, gradient, qm_gradients, mm1_gradients):
        """Fold this region's QM and MM1 gradients into the whole model's."""
        for atom in self.serials:
            gradient[atom-1] -= np.array(mm1_gradients[self.position(atom)])
            qm = qm_gradients[self.position(atom)]
            if atom not in self.link_partners:
                gradient[atom-1] += qm
            else:
                gradient[self.link_partners[atom]-1] += (1-self.g_bond[atom])*qm
                gradient[atom-1] += self.g_bond[atom]*qm
        return gradient

    def _rotate_back(self, gradients):
        """Undo the transforms, so the gradients match the untransformed model."""
        for transform in self.transforms:
            R_inv = np.linalg.inv(np.array(transform['R']))
            for atom in parse_serials(transform['atoms']):
                gradients[self.position(atom)] = np.matmul(R_inv, gradients[self.position(atom)])
        return gradients

    def _log(self, qm_energy, mm_energy, mm1_energy):
        """Append a row to this region's log; qm_energy is None for a failed run,
        whose energies belong to the previous iteration and are left out."""
        logfile = self.log_file
        # macro_cycle_width = 12
        iter_width = 5
        width = 25
        if not os.path.exists(logfile):
            header = ''
            # header += 'macro cycle'.rjust(macro_cycle_width)
            header += 'iter'.rjust(iter_width)
            # header += 'QM energy (Ha)'.rjust(width)
            # header += 'QM energy (kcal/mol)'.rjust(width)
            header += 'QM energy (kJ/mol)'.rjust(width)
            header += 'MM energy ("kJ/mol")'.rjust(width)
            header += 'MM1 energy ("kJ/mol")'.rjust(width)
            header += 'QM/MM energy ("kJ/mol")'.rjust(width)
            header += '\n'
            with open(logfile, 'w') as handle:
                handle.write(header)
        # The iteration counter continues from the last row written.
        with open(logfile, 'r') as handle:
            log = handle.readlines()
        line = ''
        # line += str(macro_cycle).rjust(macro_cycle_width)
        # if len(log) == 1 or int(log[-1].split()[0].strip()) != macro_cycle:
        #     iter = '1'.rjust(iter_width)
        # else:
        if len(log) == 1:
            iteration = '1'
        else:
            iteration = str(int(log[-1].split()[0]) + 1)
        line += iteration.rjust(iter_width)
        if qm_energy is None:
            line += 'Failed'.rjust(width)
        else:
            scale = self.w_qm*KCAL_PER_HARTREE
            # line += '{:.12f}'.format(round(qm_energy, 12)).rjust(width)
            # line += '{:.12f}'.format(round(qm_energy*harkcal, 12)).rjust(width)
            line += '{:.10f}'.format(round(qm_energy*KJ_PER_HARTREE, 10)).rjust(width)
            line += '{:.10f}'.format(round((mm_energy/scale)*KJ_PER_HARTREE, 10)).rjust(width)
            line += '{:.10f}'.format(round((mm1_energy/scale)*KJ_PER_HARTREE, 10)).rjust(width)
            line += '{:.10f}'.format(round((qm_energy + mm_energy/scale - mm1_energy/scale)*KJ_PER_HARTREE, 10)).rjust(width)
        line += '\n'
        with open(logfile, 'a') as handle:
            handle.write(line)

    def mm1_energies(self, sites_cart):
        """The model for this region, and its MM energy and gradients unscaled."""
        # at this point we need a model object for syst1 :( but we have sites_cart for model_real
        write_coordinates(infile=self.mm1_file, sites_cart=sites_cart)
        dm = DataManager()
        if self.cif is not None:
            for cif in self.cif:
                dm.process_restraint_file(str(cif))
        dm.process_model_file(self.mm1_file)
        model_mm1 = dm.get_model(filename=self.mm1_file)

        with open('settings.pickle', 'rb') as handle:
            params = pickle.load(handle)
        model_mm1.process(pdb_interpretation_params=params, make_restraints=True)
        residuals = model_mm1.restraints_manager_energies_sites(compute_gradients=True)
        residuals.target = residuals.target*(1.0/residuals.normalization_factor)
        residuals.gradients = residuals.gradients*(1.0/residuals.normalization_factor)

        # this needs to come after model_mm1.process()
        self.restore_serials(model_mm1)
        return model_mm1, residuals

    def _orca_succeeded(self, status):
        """Whether Orca exited cleanly and ran to completion."""
        if status != 0 or not os.path.exists(self.engrad_file):
            return False
        with open(self.output_file, 'r') as handle:
            return any('ORCA TERMINATED NORMALLY' in line for line in handle)

    def qm_energy_and_gradient(self, model_mm1):
        """Run Orca on this region and return its energy, rescaled gradients, and
        whether Orca finished.  A failed run writes no gradient, so the one from
        the previous iteration is used instead."""
        self.apply_transforms(model_mm1)
        self.write_qm_model(model_mm1)
        with open(self.output_file, 'w') as out:
            status = subprocess.call([self.orca_binary, self.input_file],
                stdout=out, stderr=subprocess.STDOUT)
        succeeded = self._orca_succeeded(status)
        if not succeeded and not os.path.exists(self.engrad_file):
            raise OrcaFailed('%s did not finish, see %s'
                             % (self.orca_binary, self.output_file))
        qm_energy, qm_gradients = read_engrad(self.engrad_file)
        qm_gradients = rescale_qm_gradients(qm_gradients, self.w_qm*KCAL_PER_HARTREE/ANGSTROM_PER_BOHR)
        qm_gradients = self._rotate_back(qm_gradients)
        return qm_energy, qm_gradients, succeeded

    def contribute(self, sites_cart, gradient, target, mm_residual_sum):
        """Add this region to the gradient and target, and log what it added."""
        model_mm1, mm1 = self.mm1_energies(sites_cart)
        qm_energy, qm_gradients, succeeded = self.qm_energy_and_gradient(model_mm1)

        target = target - mm1.target + self.w_qm*KCAL_PER_HARTREE*qm_energy
        gradient = self._add_qm_mm(gradient, qm_gradients, mm1.gradients)
        for restraint in self.restraints:
            target += restraint.add_to(sites_cart, gradient)

        self._log(qm_energy if succeeded else None, mm_residual_sum, mm1.target)
        return gradient, target


def run(sites_cart, mm_gradients, mm_residual_sum):
    settings = read_settings('qref.dat')

    if len(sites_cart) != settings['n_atoms']:
        return mm_gradients, mm_residual_sum

    with open('qm.lock', 'w'):
        pass

    try:
        # if specified, update coordinates of restart file
        if settings['restart'] is not None:
            write_coordinates(infile=settings['restart'], sites_cart=sites_cart)

        target = mm_residual_sum
        total_gradient = mm_gradients

        # loop over all the definitions of syst1 and process (order matters)
        for number, syst1 in enumerate(settings['syst1_files'], 1):
            region = QMRegion(number, syst1, settings)
            total_gradient, target = region.contribute(sites_cart, total_gradient,
                target, mm_residual_sum)
    finally:
        os.remove('qm.lock')

    return total_gradient, target
