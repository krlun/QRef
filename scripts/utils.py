import re
import json

import numpy as np

from cctbx.array_family import flex


def apply_transforms(model, transforms, serial_to_index):
    for transform in transforms:
        R = np.array(transform['R'])
        t = np.array(transform['t'])
        atoms_model = model.get_hierarchy().atoms()
        for atom in parse_atoms_line(transform['atoms']):
            atoms_model[serial_to_index[atom]].xyz = np.matmul(R, atoms_model[serial_to_index[atom]].xyz) + t


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


def read_junc_factors(junc_factor_file):
    delimiters = '[#!]'
    junc_factors = dict()
    with open(junc_factor_file, 'r') as file:
        line = file.readline()
        while line:
            line = re.split(delimiters, line)[0].strip()
            if len(line) > 0:
                line = line.split()
                res = line[0]
                if len(res) == 3:
                    if res not in junc_factors.keys():
                        junc_factors[res] = dict()
                    # bond = '-'.join(sorted(line[1:3]))
                    bond = '-'.join(line[1:3])
                    if bond not in junc_factors[res].keys():
                        junc_factors[res][bond] = dict()
                    line = file.readline()
                    while line != '\n' and line:
                        line = re.split(delimiters, line)[0].strip().split()
                        ltype = int(line[0])
                        distance = float(line[1].split('d')[0])
                        junc_factors[res][bond][ltype] = distance
                        line = file.readline()
            line = file.readline()
    return junc_factors


def read_dat(infile):
    with open(infile, 'r') as file:
        dat = json.load(file, object_hook=lambda d: {int(key) if key.isdigit() else key: value for key, value in d.items()})
    return dat


def select_qm_model(model, qm):
    hierarchy = model.get_hierarchy()
    sel = flex.bool(hierarchy.atoms_size())
    for atom in qm:
        sel[atom-1] = True
    return model.select(sel)


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
    hierarchy.write_pdb_file(file_name=outfile, crystal_symmetry=model.crystal_symmetry())