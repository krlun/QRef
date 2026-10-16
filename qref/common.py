"""Functions and classes shared by the QRef interface and the supporting scripts."""
from __future__ import absolute_import, division

import json
import re

import numpy as np

from cctbx.array_family import flex


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


class Region(object):
    """One syst1 file: its atoms, how its link atoms are capped, its restraints,
    and the files named after it.  Without settings, the rest is left for
    qref_prep.py to fill in."""

    def __init__(self, number, syst1, settings=None):
        self.number = number
        self.syst1 = syst1
        self.serials, self.link_atoms = read_region(syst1)
        self._positions = dict((serial, i) for i, serial in enumerate(sorted(self.serials)))
        self.link_partners = dict()
        self.g_bond = dict()
        self.transforms = list()
        self.restraints_distance = list()
        self.restraints_angle = list()
        if settings is not None:
            entry = settings.entries[syst1]
            self.link_partners = entry['link_pairs']
            self.g_bond = entry['g']
            # missing from older qref.dat files
            self.transforms = entry.get('transforms', [])
            self.restraints_distance = entry['restraints_distance']
            self.restraints_angle = entry['restraints_angle']

    def position(self, serial):
        """Where the atom is in this region's files: mm_i_c.pdb, qm_i_h.pdb and qm_i.engrad."""
        return self._positions[serial]

    @property
    def mm1_file(self):
        return 'mm_%d_c.pdb' % self.number

    @property
    def qm_file(self):
        return 'qm_%d_h.pdb' % self.number

    @property
    def input_file(self):
        return 'qm_%d.inp' % self.number

    @property
    def output_file(self):
        return 'qm_%d.out' % self.number

    @property
    def engrad_file(self):
        return 'qm_%d.engrad' % self.number

    @property
    def gbw_file(self):
        return 'qm_%d.gbw' % self.number

    @property
    def log_file(self):
        return 'qref_%d.log' % self.number

    def select(self, model):
        """The atoms of the whole model that are in this region."""
        hierarchy = model.get_hierarchy()
        sel = flex.bool(hierarchy.atoms_size())
        for atom in self.serials:
            sel[atom-1] = True
        return model.select(sel)

    def restore_serials(self, model):
        """Give the atoms of a model of this region their serials in the whole model."""
        serials = sorted(self.serials)
        atoms = model.get_hierarchy().atoms()
        width = 5
        for atom in atoms: atom.serial = str(serials[int(atom.serial) - 1]).rjust(width)

    def apply_transforms(self, model):
        for transform in self.transforms:
            R = np.array(transform['R'])
            t = np.array(transform['t'])
            atoms_model = model.get_hierarchy().atoms()
            for atom in parse_serials(transform['atoms']):
                atoms_model[self.position(atom)].xyz = np.matmul(R, atoms_model[self.position(atom)].xyz) + t

    def write_qm_model(self, model, path=None):
        """Write the model with each link atom replaced by a hydrogen, to qm_file
        unless another path is given."""
        hierarchy = model.get_hierarchy()
        atoms = hierarchy.atoms()
        for atom in atoms:
            atom_serial = int(atom.serial.strip())
            if atom.element_is_hydrogen() or atom_serial in self.link_partners:
                atom.element = ' H'
            if atom_serial in self.link_partners:
                c_qm = atoms[self.position(self.link_partners[atom_serial])]
                atom.xyz = (c_qm.xyz[0] + self.g_bond[atom_serial]*(atom.xyz[0] - c_qm.xyz[0]),
                    c_qm.xyz[1] + self.g_bond[atom_serial]*(atom.xyz[1] - c_qm.xyz[1]),
                    c_qm.xyz[2] + self.g_bond[atom_serial]*(atom.xyz[2] - c_qm.xyz[2]))
        hierarchy.write_pdb_file(file_name=path or self.qm_file, crystal_symmetry=model.crystal_symmetry(), anisou=False)

    def settings_entry(self):
        """This region's entry in qref.dat."""
        return {'link_pairs': self.link_partners,
                'g': self.g_bond,
                'transforms': self.transforms,
                'restraints_distance': self.restraints_distance,
                'restraints_angle': self.restraints_angle}


class Settings(object):
    """qref.dat: what qref_prep.py writes, and QRef and the scripts read back.
    entries holds each region's entry, keyed by its syst1 file."""

    def __init__(self, syst1_files, n_atoms=None, restart=None, w_qm=None, cif=None,
                 orca_binary=None, ltype=None):
        self.syst1_files = syst1_files
        self.n_atoms = n_atoms
        self.restart = restart
        self.w_qm = w_qm
        self.cif = cif
        self.orca_binary = orca_binary
        self.ltype = ltype
        self.entries = dict()

    @classmethod
    def read(cls, path='qref.dat'):
        # the serial keys of link_pairs and g back to int
        with open(path, 'r') as handle:
            data = json.load(handle, object_hook=lambda d: {int(key) if key.isdigit() else key: value for key, value in d.items()})
        settings = cls(data['syst1_files'], n_atoms=data['n_atoms'], restart=data['restart'],
                       w_qm=data['w_qm'], cif=data['cif'], orca_binary=data['orca_binary'],
                       ltype=data['ltype'])
        for syst1 in settings.syst1_files:
            settings.entries[syst1] = data[syst1]
        return settings

    def write(self, path='qref.dat'):
        data = dict(self.entries)
        data.update({'syst1_files': self.syst1_files, 'n_atoms': self.n_atoms,
                     'restart': self.restart, 'w_qm': self.w_qm, 'cif': self.cif,
                     'orca_binary': self.orca_binary, 'ltype': self.ltype})
        with open(path, 'w') as handle:
            json.dump(data, handle, indent=4, sort_keys=True)

    def add(self, region):
        self.entries[region.syst1] = region.settings_entry()

    def regions(self, region_class=Region):
        """The regions, in the order of syst1_files."""
        return [region_class(number, syst1, self) for number, syst1 in enumerate(self.syst1_files, 1)]
