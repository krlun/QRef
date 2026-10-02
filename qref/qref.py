from __future__ import absolute_import, division
import os
import sys
import subprocess
import pickle
import re
import json

import numpy as np

from iotbx.data_manager import DataManager

from qref.common import apply_transforms
from qref.common import convert_serial_to_index
from qref.common import parse_atoms_line
from qref.common import read_dat
from qref.common import read_syst1
from qref.common import restore_serial_in_model
from qref.common import write_pdb_h

harkcal = 627.509474063112
harkJ = 2625.499639479950
bohrang = 0.529177249


def read_energy_and_gradient_from_orca(infile):
    gradients = list()
    with open(infile, 'r') as file:
        line = file.readline()
        while line:
            if line == '# Number of atoms\n':
                file.readline()
                n_atoms = int(file.readline().strip())
            if line == '# The current total energy in Eh\n':
                file.readline()
                energy = float(file.readline().strip())
                break
            line = file.readline()
        for i in range(3): file.readline()
        for i in range(n_atoms):
            g = (float(file.readline()), float(file.readline()), float(file.readline()))
            gradients.append(g)
    return energy, gradients


def rescale_qm_gradients(qm_gradients, w):
    return [tuple([w*component for component in gradient]) for gradient in qm_gradients]


def update_file_coordinates(infile, sites_cart):
    with open(infile, 'r') as file:
        model = file.readlines()
    records = {'ATOM', 'HETATM'}
    width = 8
    with open(infile, 'w') as file:
        for line in model:
            if line[0:6].strip() in records:
                serial = int(line[6:11].strip())
                coords = ''
                for i in range(3): coords += '{:.3f}'.format(round(sites_cart[serial-1][i], 3)).rjust(width)
                line = line[0:30] + coords + line[54:]
            file.write(line)


class Region(object):
    """One syst1 definition: its atoms, its junctions and the files named after it."""

    def __init__(self, index, syst1, dat):
        self.index = index
        self.qm_atoms = read_syst1(syst1)[0]
        self.serial_to_index = convert_serial_to_index(self.qm_atoms)
        self.link_pairs = dat[syst1]['link_pairs']
        self.g = dat[syst1]['g']
        self.transforms = dat[syst1]['transforms']
        self.restraints_distance = dat[syst1]['restraints_distance']
        self.restraints_angle = dat[syst1]['restraints_angle']
        self.cif = dat['cif']
        self.orca_binary = dat['orca_binary']
        self.w_qm = dat['w_qm']

    @property
    def mm1_file(self):
        return 'mm_%d_c.pdb' % self.index

    @property
    def qm_file(self):
        return 'qm_%d_h.pdb' % self.index

    @property
    def input_file(self):
        return 'qm_%d.inp' % self.index

    @property
    def output_file(self):
        return 'qm_%d.out' % self.index

    @property
    def engrad_file(self):
        return 'qm_%d.engrad' % self.index

    @property
    def log_file(self):
        return 'qref_%d.log' % self.index

    def _add_qm_mm(self, gradient, qm_gradients, mm1_gradients):
        """Fold this region's QM and MM1 gradients into the whole model's."""
        for atom in self.qm_atoms:
            gradient[atom-1] -= np.array(mm1_gradients[self.serial_to_index[atom]])
            if atom not in self.link_pairs.keys():
                gradient[atom-1] += np.array(qm_gradients[self.serial_to_index[atom]])
            else:
                gradient[self.link_pairs[atom]-1] += (1-self.g[atom])*np.array(qm_gradients[self.serial_to_index[atom]])
                gradient[atom-1] += self.g[atom]*np.array(qm_gradients[self.serial_to_index[atom]])
        return gradient

    def _rotate_back(self, gradients):
        """Undo the transforms, so the gradients match the untransformed model."""
        for transform in self.transforms:
            R_inv = np.linalg.inv(np.array(transform['R']))
            for atom in parse_atoms_line(transform['atoms']):
                gradients[self.serial_to_index[atom]] = np.matmul(R_inv, gradients[self.serial_to_index[atom]])
        return gradients

    def _apply_restraints_distance(self, sites_cart, gradients, target):
        # restraint[0] = atom1_serial, restraint[1] = atom2_serial, restraint[2] = desired distance in Angstrom, restraint[3] = force constant
        for restraint in self.restraints_distance:
            r_ij = np.array(sites_cart[restraint[0]-1]) - np.array(sites_cart[restraint[1]-1])
            r = np.sqrt(np.sum(r_ij**2))
            delta_r = r - restraint[2]
            d_U_ij_d_r = 2*restraint[3]*delta_r
            d_r_d_r_ij = r_ij/r
            gradients[restraint[0]-1] += d_U_ij_d_r*d_r_d_r_ij
            gradients[restraint[1]-1] += -d_U_ij_d_r*d_r_d_r_ij
            target += restraint[3]*delta_r**2
        return gradients, target

    def _apply_restraints_angle(self, sites_cart, gradients, target):
        # restraint[0] = atom1_serial (i), restraint[1] = atom2_serial ("middle" atom) (j), restraint[2] = atom3_serial (k), restraint[3] = desired angle in degrees, restraint[4] = force constant
        for restraint in self.restraints_angle:
            r_ij = np.array(sites_cart[restraint[0]-1]) - np.array(sites_cart[restraint[1]-1])
            r_kj = np.array(sites_cart[restraint[2]-1]) - np.array(sites_cart[restraint[1]-1])
            norm_r_ij = np.linalg.norm(r_ij)
            norm_r_kj = np.linalg.norm(r_kj)
            cos_alpha = np.dot(r_ij, r_kj)/(norm_r_ij * norm_r_kj)
            alpha = np.degrees(np.arccos(cos_alpha))
            delta_alpha = alpha - restraint[3]
            d_U_d_alpha = 2*restraint[4]*delta_alpha
            d_alpha_d_r_i = (1.0/np.sqrt(1 - cos_alpha**2)) * (1.0/norm_r_ij) * (cos_alpha * (r_ij/norm_r_ij) - (r_kj/norm_r_kj))
            d_alpha_d_r_k = (1.0/np.sqrt(1 - cos_alpha**2)) * (1.0/norm_r_kj) * (cos_alpha * (r_kj/norm_r_kj) - (r_ij/norm_r_ij))
            d_alpha_d_r_j = - d_alpha_d_r_i - d_alpha_d_r_k
            gradients[restraint[0]-1] += d_U_d_alpha*d_alpha_d_r_i
            gradients[restraint[1]-1] += d_U_d_alpha*d_alpha_d_r_j
            gradients[restraint[2]-1] += d_U_d_alpha*d_alpha_d_r_k
            target += restraint[4]*delta_alpha**2
        return gradients, target

    def _log(self, qm_energy, mm_energy, mm1_energy):
        """Append a row of energies to this region's log."""
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
            with open(logfile, 'w') as file:
                file.write(header)
        # The iteration counter continues from the last row written.
        with open(logfile, 'r') as file:
            log = file.readlines()
        line = ''
        # line += str(macro_cycle).rjust(macro_cycle_width)
        # if len(log) == 1 or int(log[-1].split()[0].strip()) != macro_cycle:
        #     iter = '1'.rjust(iter_width)
        # else:
        if len(log) == 1:
            iter = '1'
        else:
            iter = str(int(log[-1].split()[0]) + 1)
        line += iter.rjust(iter_width)
        with open(self.output_file, 'r') as file:
            terminated = any('ORCA TERMINATED NORMALLY' in out_line for out_line in file)
        if not terminated:
            line += 'Failed'.rjust(width)
        else:
            scale = self.w_qm*harkcal
            # line += '{:.12f}'.format(round(qm_energy, 12)).rjust(width)
            # line += '{:.12f}'.format(round(qm_energy*harkcal, 12)).rjust(width)
            line += '{:.10f}'.format(round(qm_energy*harkJ, 10)).rjust(width)
            line += '{:.10f}'.format(round((mm_energy/scale)*harkJ, 10)).rjust(width)
            line += '{:.10f}'.format(round((mm1_energy/scale)*harkJ, 10)).rjust(width)
            line += '{:.10f}'.format(round((qm_energy + mm_energy/scale - mm1_energy/scale)*harkJ, 10)).rjust(width)
        line += '\n'
        with open(logfile, 'a') as file:
            file.write(line)

    def mm1_energies(self, sites_cart):
        """The model for this region, and its MM energy and gradients unscaled."""
        # at this point we need a model object for syst1 :( but we have sites_cart for model_real
        update_file_coordinates(infile=self.mm1_file, sites_cart=sites_cart)
        dm = DataManager()
        if self.cif is not None:
            for cif in self.cif:
                dm.process_restraint_file(str(cif))
        dm.process_model_file(self.mm1_file)
        model_mm1 = dm.get_model(filename=self.mm1_file)

        with open('settings.pickle', 'rb') as file:
            params = pickle.load(file)
        model_mm1.process(pdb_interpretation_params=params, make_restraints=True)
        residuals = model_mm1.restraints_manager_energies_sites(compute_gradients=True)
        residuals.target = residuals.target*(1.0/residuals.normalization_factor)
        residuals.gradients = residuals.gradients*(1.0/residuals.normalization_factor)

        # this needs to come after model_mm1.process()
        restore_serial_in_model(model_mm1, self.serial_to_index)
        return model_mm1, residuals

    def qm_energy_and_gradient(self, model_mm1):
        """Run Orca on this region and return its energy and rescaled gradients."""
        apply_transforms(model_mm1, self.transforms, self.serial_to_index)
        write_pdb_h(self.qm_file, model_mm1, link_pairs=self.link_pairs, g=self.g,
            serial_to_index=self.serial_to_index)
        subprocess.check_call([self.orca_binary, self.input_file],
            stdout=open(self.output_file, 'w'), stderr=subprocess.STDOUT)
        qm_energy, qm_gradients = read_energy_and_gradient_from_orca(self.engrad_file)
        qm_gradients = rescale_qm_gradients(qm_gradients, self.w_qm*harkcal/bohrang)
        qm_gradients = self._rotate_back(qm_gradients)
        return qm_energy, qm_gradients

    def contribute(self, sites_cart, gradient, target, mm_residual_sum):
        """Add this region to the gradient and target, and log what it added."""
        model_mm1, mm1 = self.mm1_energies(sites_cart)
        qm_energy, qm_gradients = self.qm_energy_and_gradient(model_mm1)

        target = target - mm1.target + self.w_qm*harkcal*qm_energy
        gradient = self._add_qm_mm(gradient, qm_gradients, mm1.gradients)
        gradient, target = self._apply_restraints_distance(sites_cart, gradient, target)
        gradient, target = self._apply_restraints_angle(sites_cart, gradient, target)

        self._log(qm_energy, mm_residual_sum, mm1.target)
        return gradient, target


def run(sites_cart, mm_gradients, mm_residual_sum):
    dat = read_dat('qref.dat')

    if not len(sites_cart) == dat['n_atoms']:
        return mm_gradients, mm_residual_sum

    # establish lock
    with open('qm.lock', 'w'):
        pass

    # if specified, update coordinates of restart file
    if dat['restart'] is not None:
        update_file_coordinates(infile=dat['restart'], sites_cart=sites_cart)

    target = mm_residual_sum
    total_gradient = mm_gradients
    
    # loop over all the definitions of syst1 and process (order matters)
    for index, syst1 in enumerate(dat['syst1_files'], 1):
        region = Region(index, syst1, dat)
        total_gradient, target = region.contribute(sites_cart, total_gradient,
            target, mm_residual_sum)

    # unlock
    os.remove('qm.lock')

    return total_gradient, target
