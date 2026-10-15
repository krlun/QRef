#!/usr/bin/env cctbx.python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
import json
import subprocess

from iotbx.data_manager import DataManager

from utils import apply_transforms
from utils import convert_serial_to_index
from utils import read_syst1
from utils import read_junc_factors
from utils import restore_serial_in_model
from utils import select_qm_model
from utils import write_pdb_h


def identify_link_pairs(model, link_atoms, serial_to_index):
    link_pairs = dict()
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    bond_params_table = model.restraints_manager.geometry.bond_params_table
    for link in link_atoms:
        min_distance = float('inf')
        closest = None
        link_atom = atoms[serial_to_index[link]]
        for atom in atoms:
            serial = int(atom.serial.strip())
            if serial != link and bond_params_table.lookup(serial_to_index[link], serial_to_index[serial]) is not None:
                distance = atom.distance(link_atom)
                if distance < min_distance:
                    closest = serial
                    min_distance = distance
        if closest is None:
            raise SystemExit('No atom bonded to link atom ' + str(link) + ' in the QM system. Exiting...')
        link_pairs[link] = closest
    return link_pairs


def check_overlap(syst1_files):
    # E_MM counts the restraints among atoms that two QM regions share once, but
    # each region subtracts them again in E_MM1. Link atoms may be shared: there
    # are no restraints among them, and Orca only sees the hydrogens replacing them.
    regions = dict()
    qm_only = set()
    for syst1 in syst1_files:
        qm_atoms, link_atoms = read_syst1(syst1)
        shared = sorted(atom for atom in qm_atoms if atom in regions and (atom in qm_only or atom not in link_atoms))
        if shared:
            others = sorted(set(other for atom in shared for other in regions[atom]))
            raise SystemExit(syst1 + ' shares atoms ' + ', '.join(str(atom) for atom in shared) + ' with ' + ', '.join(others) + '. QM regions must not share atoms other than link atoms; duplicate the QM region into separate altlocs instead. Exiting...')
        for atom in qm_atoms:
            regions.setdefault(atom, []).append(syst1)
        qm_only |= qm_atoms - link_atoms


def calculate_g_factor(model, link_pairs, junc_factors, ltype, serial_to_index):
    g = dict()
    atoms = model.get_hierarchy().atoms()
    for key, value in link_pairs.items():
        c_qm = atoms[serial_to_index[value]]
        link_atom = atoms[serial_to_index[key]]
        resname = link_atom.parent().resname.strip()
        bond = '-'.join([c_qm.name.strip(), link_atom.name.strip()])
        if resname not in junc_factors:
            raise SystemExit(resname + ' missing in junction factor file. Exiting...')
        if bond not in junc_factors[resname]:
            raise SystemExit(bond + ' (between atoms ' + str(value) + ' and ' + str(key) + ') missing for ' + resname + ' in junction factor file. Exiting...')
        if ltype not in junc_factors[resname][bond]:
            raise SystemExit('Missing ltype ' + str(ltype) + ' for bond ' + bond + ' in residue ' + resname + ' in junction factor file. Exiting...')
        c_h_ideal = junc_factors[resname][bond][ltype] # eg. print(junc_factors['PHE']['CB-CG'][12])
        g[key] = c_h_ideal/model.restraints_manager.geometry.bond_params_table.lookup(serial_to_index[key], serial_to_index[value]).distance_ideal
    return g


def suggest_selection_string(model):
    hierarchy = model.get_hierarchy()
    selection_string = '('
    for chain in hierarchy.only_model().chains():
        selection_string += '('
        for residue in chain.residue_groups():
            selection_string += 'resseq ' + residue.resseq.strip() + ' or '
        selection_string = selection_string[:-4] + ') and chain ' + chain.id + ') or ('
    selection_string = selection_string[:-5]
    print()
    print('Suggested selection string (reciprocal space):')
    print()
    print('\"' + selection_string + '\"')
    print()
    print('Suggested selection string (real space):')
    print()
    print('\"all and not (' + selection_string + ')\"')
    print()


def check_altlocs(model, link_atoms):
    # link atoms can be shared between the regions of different altlocs
    hierarchy = model.get_hierarchy()
    altlocs = set()
    for atom in hierarchy.atoms():
        if int(atom.serial) not in link_atoms: altlocs.add(atom.pdb_label_columns()[4])
    altlocs = sorted(list(altlocs))
    if len(altlocs) != 1:
        if altlocs[0] == ' ': altlocs[0] = '<blank>'
        print('Warning, multiple altlocs for subsystem found:  ' + ', '.join(altlocs))
    elif altlocs[0] != ' ':
        print('Altloc:  \"and altloc ' + altlocs[0] + '\"')


def write_dat(dat):
    with open('qref.dat', 'w') as file:
        json.dump(dat, file, indent=4, sort_keys=True)


def prepare_restart(infile, outfile):
    records = {'ATOM', 'HETATM', 'ANISOU', 'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3', 'TER', 'END'}
    with open(infile, 'r') as file, open(outfile, 'w') as out:
        line = file.readline()
        while line:
            if line[0:6].strip() in records:
                out.write(line)
            line = file.readline()


def locate_binary(binary):
    proc = subprocess.Popen(['which', binary], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    (out, err) = proc.communicate()
    if not isinstance(out, str):    # bytes under python 3, already str under 2
        out = out.decode()
    out = out.rstrip()
    return out if out else binary


def parse_args():
    parser = argparse.ArgumentParser(
        description='Prepares for quantum refinement: writes the model of each QM region (mm_i_c.pdb and qm_i_h.pdb) and the settings for QRef (qref.dat).'
    )
    parser.add_argument('pdb', type=str, help='PDB file describing the entire model')
    parser.add_argument('-c', '--cif', nargs='+', type=str, metavar='CIF', help='CIF files with restraints for residues and ligands')
    parser.add_argument('-s', '--syst1', nargs='+', default=['syst1'], type=str, metavar='SYST1', help='files defining the QM regions, one per region (default: syst1)')
    parser.add_argument('-u', '--skip_h', action='store_true', help='only write mm_i_c.pdb')
    parser.add_argument('-j', '--junctfactor', default='junctfactor', type=str, metavar='FILE', help='file with the ideal C-H distances for the link atoms (default: junctfactor)')
    parser.add_argument('-l', '--ltype', default=12, type=int, metavar='LTYPE', help='level of theory for the C-H distances in the junctfactor file (default: 12, TPSS(RI)/def2-SV(P))')
    parser.add_argument('-w', '--w_qm', default=7.5, type=float, metavar='W_QM', help='weight of the QM energy and gradients, in mol/kcal (default: 7.5)')
    parser.add_argument('-r', '--restart', nargs='?', const='restart.pdb', type=str, metavar='FILE', help='write a restart file (default name: restart.pdb)')
    parser.add_argument('-rd', '--restraint_distance', nargs=5, action='append', type=str, metavar=('I', 'ATOM1', 'ATOM2', 'DISTANCE', 'k'), help='harmonic distance restraint in QM region I between atoms ATOM1 and ATOM2 (serial numbers), with the ideal DISTANCE in Å and force constant k; can be repeated')
    parser.add_argument('-ra', '--restraint_angle', nargs=6, action='append', type=str, metavar=('I', 'ATOM1', 'ATOM2', 'ATOM3', 'ANGLE', 'k'), help='harmonic angle restraint in QM region I on the angle ATOM1-ATOM2-ATOM3 (serial numbers), with the ideal ANGLE in degrees and force constant k; can be repeated')
    parser.add_argument('-t', '--transform', nargs=14, action='append', type=str, metavar=('I', 'ATOMS', 'R11', 'R12', 'R13', 'R21', 'R22', 'R23', 'R31', 'R32', 'R33', 'T1', 'T2', 'T3'), help='symmetry transform of ATOMS in QM region I, given in quotes in the syntax of syst1 files, by the rotation matrix R, row-wise, and the translation T, both Cartesian; can be repeated')
    parser.add_argument('-v', '--version', action="version", version="%(prog)s 0.3.1")
    return parser.parse_args()


def main():
    args = parse_args()
    print('*** This is qref_prep.py ***')
    print()
    model_file = args.pdb
    syst1_files = args.syst1
    junc_factor_file = args.junctfactor
    ltype = args.ltype
    if not args.skip_h:
        junc_factors = read_junc_factors(junc_factor_file=junc_factor_file)
        check_overlap(syst1_files)
    dm = DataManager()
    if args.cif is not None:
        for cif in args.cif:
            dm.process_restraint_file(cif)
    dm.process_model_file(model_file)
    model_mm = dm.get_model(filename=model_file)
    model_mm.add_crystal_symmetry_if_necessary()

    dat = dict()

    for index, syst1 in enumerate(syst1_files, 1):
        print('----------')
        print(syst1.center(10))
        print('----------')
        print()
        qm_atoms, link_atoms = read_syst1(syst1)
        serial_to_index = convert_serial_to_index(qm_atoms)
        model_mm1 = select_qm_model(model=model_mm, qm=qm_atoms)
        model_mm1.process(pdb_interpretation_params=model_mm1.get_current_pdb_interpretation_params(), make_restraints=True)
        restore_serial_in_model(model_mm1, serial_to_index=serial_to_index)
        mm1_file = 'mm_' + str(index) + '_c.pdb'
        print('Writing file:  ' + mm1_file)
        model_mm1.get_hierarchy().write_pdb_file(file_name=mm1_file, crystal_symmetry=model_mm.crystal_symmetry(), anisou=False)

        if not args.skip_h:
            dat[syst1] = dict()
            link_pairs = identify_link_pairs(model_mm1, link_atoms, serial_to_index)
            g = calculate_g_factor(model_mm1, link_pairs=link_pairs, junc_factors=junc_factors, ltype=ltype, serial_to_index=serial_to_index)
            dat[syst1]['link_pairs'] = link_pairs
            dat[syst1]['g'] = g
            dat[syst1]['restraints_distance'] = list()
            dat[syst1]['restraints_angle'] = list()
            dat[syst1]['transforms'] = list()
            if args.restraint_distance is not None:
                for restraint in args.restraint_distance:
                    if int(restraint[0]) == index:
                        dat[syst1]['restraints_distance'].append([int(restraint[1]), int(restraint[2]), float(restraint[3]), float(restraint[4])])
            if args.restraint_angle is not None:
                for restraint in args.restraint_angle:
                    if int(restraint[0]) == index:
                        dat[syst1]['restraints_angle'].append([int(restraint[1]), int(restraint[2]), int(restraint[3]), float(restraint[4]), float(restraint[5])])
            if args.transform is not None:
                for transform in args.transform:
                    if int(transform[0]) == index:
                        dat[syst1]['transforms'].append(dict())
                        dat[syst1]['transforms'][-1]['atoms'] = transform[1]
                        dat[syst1]['transforms'][-1]['R'] = [[float(x) for x in transform[3*i+2:3*i+5]] for i in range(3)]
                        dat[syst1]['transforms'][-1]['t'] = [float(x) for x in transform[11:14]]
                apply_transforms(model_mm1, dat[syst1]['transforms'], serial_to_index)
            qm_file = 'qm_' + str(index) + '_h.pdb'
            print('Writing file:  ' + qm_file)
            write_pdb_h(qm_file, model_mm1, link_pairs, g, serial_to_index)

        suggest_selection_string(model_mm1)
        check_altlocs(model_mm1, link_atoms)
        print('----------')

    if not args.skip_h:
        dat['n_atoms'] = model_mm.get_number_of_atoms()
        dat['restart'] = args.restart
        dat['w_qm'] = args.w_qm
        dat['cif'] = args.cif
        dat['orca_binary'] = locate_binary('orca')
        dat['ltype'] = ltype
        dat['syst1_files'] = args.syst1
        print('----------')
        print()
        print('Writing file:  qref.dat')
        print()
        print('----------')
        write_dat(dat)

    if args.restart is not None:
        print()
        print('Writing file:  ' + args.restart)
        print()
        print('----------')
        prepare_restart(model_file, args.restart)

    # cleanup
    temp_files = ['settings.pickle']
    for file in temp_files:
        if os.path.exists(file):
            os.remove(file)


if (__name__ == "__main__"):
    main()
