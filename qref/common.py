"""Functions shared by the QRef interface and the supporting scripts."""
from __future__ import absolute_import, division

import json
import re

import numpy as np


def parse_serials(line):
    comments = '[#!]'
    delimiters = '[^,\\s]+'
    atoms = set()
    line = re.findall(delimiters, re.split(comments, line)[0])
    for interval in line:
        interval = [int(x) for x in interval.split('-')]
        for i in range(min(interval), max(interval) + 1):
            atoms.add(i)
    return atoms


def read_region(infile):
    serials = set()
    link_atoms = set()
    with open(infile, 'r') as handle:
        line = handle.readline()
        while line:
            atoms = parse_serials(line)
            for atom in atoms:
                link_atoms.add(atom) if atom in serials else serials.add(atom)
            line = handle.readline()
    return serials, link_atoms


def convert_serial_to_index(qm):
    qm_sorted = sorted(qm)
    indices = dict()
    for i in range(len(qm_sorted)):
        indices[qm_sorted[i]] = i
    return indices


def write_qm_model(outfile, model, link_partners, g_bond, serial_to_index):
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    for atom in atoms:
        atom_serial = int(atom.serial.strip())
        if atom.element_is_hydrogen() or atom_serial in link_partners:
            atom.element = ' H'
        if atom_serial in link_partners:
            c_qm = atoms[serial_to_index[link_partners[atom_serial]]]
            atom.xyz = (c_qm.xyz[0] + g_bond[atom_serial]*(atom.xyz[0] - c_qm.xyz[0]), 
                c_qm.xyz[1] + g_bond[atom_serial]*(atom.xyz[1] - c_qm.xyz[1]),
                c_qm.xyz[2] + g_bond[atom_serial]*(atom.xyz[2] - c_qm.xyz[2]))
    hierarchy.write_pdb_file(file_name=outfile, crystal_symmetry=model.crystal_symmetry(), anisou=False)


def restore_serials(model, serial_to_index):
    index_to_serial = {value: key for key, value in serial_to_index.items()}
    atoms = model.get_hierarchy().atoms()
    width = 5
    for atom in atoms: atom.serial = str(index_to_serial[int(atom.serial) - 1]).rjust(width)


def read_settings(infile):
    with open(infile, 'r') as handle:
        settings = json.load(handle, object_hook=lambda d: {int(key) if key.isdigit() else key: value for key, value in d.items()})
    return settings


def apply_transforms(model, transforms, serial_to_index):
    for transform in transforms:
        R = np.array(transform['R'])
        t = np.array(transform['t'])
        atoms_model = model.get_hierarchy().atoms()
        for atom in parse_serials(transform['atoms']):
            atoms_model[serial_to_index[atom]].xyz = np.matmul(R, atoms_model[serial_to_index[atom]].xyz) + t
