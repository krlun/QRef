from __future__ import division
import os
import sys
import subprocess
import pickle
import re
import json

import numpy as np

from iotbx.data_manager import DataManager

harkcal = 627.509474063112
harkJ = 2625.499639479950
bohrang = 0.529177249


def parse_atoms_line(line):
    comments = '[#!]'
    delimiters = '[^,\\s]+'
    atoms = set()
    line = re.findall(delimiters, re.split(comments, line)[0])
    for interval in line:
        interval = [int(x) for x in interval.split('-')]
        for i in range(min(interval), max(interval) + 1):
            atoms.add(i)
    return atoms


def read_syst1(infile):
    qm_atoms = set()
    link_atoms = set()
    with open(infile, 'r') as file:
        line = file.readline()
        while line:
            atoms = parse_atoms_line(line)
            for atom in atoms:
                link_atoms.add(atom) if atom in qm_atoms else qm_atoms.add(atom)
            line = file.readline()
    return qm_atoms, link_atoms


def convert_serial_to_index(qm):
    qm_sorted = sorted(qm)
    indices = dict()
    for i in range(len(qm_sorted)):
        indices[qm_sorted[i]] = i
    return indices


def write_pdb_h(outfile, model, link_pairs, g, serial_to_index):
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    for atom in atoms:
        atom_serial = int(atom.serial.strip())
        if atom.element_is_hydrogen() or atom_serial in link_pairs.keys():
            atom.element = ' H'
        if atom_serial in link_pairs.keys():
            c_qm = atoms[serial_to_index[link_pairs[atom_serial]]]
            atom.xyz = (c_qm.xyz[0] + g[atom_serial]*(atom.xyz[0] - c_qm.xyz[0]), 
                c_qm.xyz[1] + g[atom_serial]*(atom.xyz[1] - c_qm.xyz[1]),
                c_qm.xyz[2] + g[atom_serial]*(atom.xyz[2] - c_qm.xyz[2]))
    hierarchy.write_pdb_file(file_name=outfile)


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


def calculate_total_gradient(qm_gradients, mm1_gradients, mm_gradients, qm_atoms, g, link_pairs, serial_to_index):
    for atom in qm_atoms:
        mm_gradients[atom-1] -= np.array(mm1_gradients[serial_to_index[atom]])
        if atom not in link_pairs.keys():
            mm_gradients[atom-1] += np.array(qm_gradients[serial_to_index[atom]])
        else:
            mm_gradients[link_pairs[atom]-1] += (1-g[atom])*np.array(qm_gradients[serial_to_index[atom]])
            mm_gradients[atom-1] += g[atom]*np.array(qm_gradients[serial_to_index[atom]])
    return mm_gradients


def rescale_qm_gradients(qm_gradients, w):
    return [tuple([w*component for component in gradient]) for gradient in qm_gradients]


def calculate_target(qm_energy, mm1_target, mm_target, w_qm):
    return mm_target - mm1_target + w_qm*qm_energy


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


def restore_serial_in_model(model, serial_to_index):
    index_to_serial = {value: key for key, value in serial_to_index.items()}
    atoms = model.get_hierarchy().atoms()
    width = 5
    for atom in atoms: atom.serial = str(index_to_serial[int(atom.serial) - 1]).rjust(width)


def logging(index, w_qm, qm_energy, mm_energy, mm1_energy):
    logfile = 'qref_' + str(index) + '.log'
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
    with open(logfile, 'r') as file:
        log = file.readlines()
    with open(logfile, 'w') as file:
        file.writelines(log)
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
        if subprocess.call(['grep', '-q', 'ORCA TERMINATED NORMALLY', 'qm_' + str(index) + '.out'], shell=False) != 0:
            line += 'Failed'.rjust(width)
        else:
            # line += '{:.12f}'.format(round(qm_energy, 12)).rjust(width)
            # line += '{:.12f}'.format(round(qm_energy*harkcal, 12)).rjust(width)
            line += '{:.10f}'.format(round(qm_energy*harkJ, 10)).rjust(width)
            line += '{:.10f}'.format(round((mm_energy/(w_qm*harkcal))*harkJ, 10)).rjust(width)
            line += '{:.10f}'.format(round((mm1_energy/(w_qm*harkcal))*harkJ, 10)).rjust(width)
            line += '{:.10f}'.format(round((qm_energy + mm_energy/(w_qm*harkcal) - mm1_energy/(w_qm*harkcal))*harkJ, 10)).rjust(width)
        line += '\n'
        file.write(line)
    return


def read_dat():
    with open('qref.dat', 'r') as file:
        dat = json.load(file, object_hook=lambda d: {int(key) if key.isdigit() else key: value for key, value in d.items()})
    return dat


def apply_restraints_distance(sites_cart, gradients, target, restraints):
    # restraint[0] = atom1_serial, restraint[1] = atom2_serial, restraint[2] = desired distance in Angstrom, restraint[3] = force constant
    for restraint in restraints:
        r_ij = np.array(sites_cart[restraint[0]-1]) - np.array(sites_cart[restraint[1]-1])
        r = np.sqrt(np.sum(r_ij**2))
        delta_r = r - restraint[2]
        d_U_ij_d_r = 2*restraint[3]*delta_r
        d_r_d_r_ij = r_ij/r
        gradients[restraint[0]-1] += d_U_ij_d_r*d_r_d_r_ij
        gradients[restraint[1]-1] += -d_U_ij_d_r*d_r_d_r_ij
        target += restraint[3]*delta_r**2
    return gradients, target


def apply_restraints_angle(sites_cart, gradients, target, restraints):
    # restraint[0] = atom1_serial (i), restraint[1] = atom2_serial ("middle" atom) (j), restraint[2] = atom3_serial (k), restraint[3] = desired angle in degrees, restraint[4] = force constant
    for restraint in restraints:
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


def apply_transforms(model, transforms, serial_to_index):
    for transform in transforms:
        R = np.array(transform['R'])
        t = np.array(transform['t'])
        atoms_model = model.get_hierarchy().atoms()
        for atom in parse_atoms_line(transform['atoms']):
            atoms_model[serial_to_index[atom]].xyz = np.matmul(R, atoms_model[serial_to_index[atom]].xyz) + t


def rotate_gradients(gradients, transforms, serial_to_index):
    for transform in transforms:
        R_inv = np.linalg.inv(np.array(transform['R']))
        for atom in parse_atoms_line(transform['atoms']):
            gradients[serial_to_index[atom]] = np.matmul(R_inv, gradients[serial_to_index[atom]])
    return gradients


def run(sites_cart, mm_gradients, mm_residual_sum):
    dat = read_dat()

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
        # read syst1 file, which containts QM system + link atoms
        qm_atoms, link_atoms = read_syst1(syst1)

        # construct a dict {serial:index} for the indices of the link atoms in the qm system
        serial_to_index = convert_serial_to_index(qm_atoms)

        # at this point we need a model object for syst1 :( but we have sites_cart for model_real
        mm1_file = 'mm_' + str(index) + '_c.pdb'
        qm_file = 'qm_' + str(index) + '_h.pdb'
        update_file_coordinates(infile=mm1_file, sites_cart=sites_cart)
        dm = DataManager()
        if dat['cif'] is not None:
            for cif in dat['cif']:
                dm.process_restraint_file(str(cif))
        dm.process_model_file(mm1_file)
        model_mm1 = dm.get_model(filename=mm1_file)

        # calculate mm gradients and target for model system
        with open('settings.pickle', 'rb') as file:
            params = pickle.load(file)
        model_mm1.process(pdb_interpretation_params=params, make_restraints=True)
        model_mm1_residuals = model_mm1.restraints_manager_energies_sites(compute_gradients=True)

        # unscale target and gradients in model_mm1
        model_mm1_residuals.target = model_mm1_residuals.target*(1.0/model_mm1_residuals.normalization_factor)
        model_mm1_residuals.gradients = model_mm1_residuals.gradients*(1.0/model_mm1_residuals.normalization_factor)

        # restore original serial to model object; this needs to come after model_mm1.process()
        restore_serial_in_model(model_mm1, serial_to_index)

        # read link_pairs and g
        link_pairs = dat[syst1]['link_pairs']
        g = dat[syst1]['g']

        # prepare for orca, transform model_mm1 if necessary
        apply_transforms(model_mm1, dat[syst1]['transforms'], serial_to_index)
        write_pdb_h(qm_file, model_mm1, link_pairs=link_pairs, g=g, serial_to_index=serial_to_index)

        # run orca
        orca_string = dat['orca_binary'] + ' qm_' + str(index) + '.inp > qm_' + str(index) + '.out'
        os.system(orca_string)

        # read results from orca
        qmengrad = 'qm_' + str(index) + '.engrad'
        qm_energy, qm_gradients = read_energy_and_gradient_from_orca(qmengrad)

        # rescale qm gradients
        qm_gradients = rescale_qm_gradients(qm_gradients, dat['w_qm']*harkcal/bohrang)

        # rotate gradients back if needed
        qm_gradients = rotate_gradients(qm_gradients, dat[syst1]['transforms'], serial_to_index)

        # update target
        target = target - model_mm1_residuals.target + dat['w_qm']*harkcal*qm_energy

        # update gradient (QM/MM)
        total_gradient = calculate_total_gradient(qm_gradients, model_mm1_residuals.gradients, total_gradient, qm_atoms, g, link_pairs, serial_to_index)

        # update gradient (restraints)
        total_gradient, target = apply_restraints_distance(sites_cart, total_gradient, target, dat[syst1]['restraints_distance'])
        total_gradient, target = apply_restraints_angle(sites_cart, total_gradient, target, dat[syst1]['restraints_angle'])

        # perform logging
        logging(index=index, w_qm=dat['w_qm'], qm_energy=qm_energy, mm_energy=mm_residual_sum, mm1_energy=model_mm1_residuals.target)

    # unlock
    os.remove('qm.lock')

    return total_gradient, target
