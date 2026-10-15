#!/usr/bin/env cctbx.python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import os
import re
import sys
import argparse
import subprocess

from iotbx.data_manager import DataManager

from utils import Region
from utils import Settings
from utils import read_region


def read_junctfactor(junc_factor_file):
    # a line 'residue atom1 atom2 type' starts a bond, followed by lines
    # 'ltype distance' up to a blank line; distances are written as Fortran
    # doubles, e.g. 1.1070d0
    junctfactor = dict()
    distances = None
    with open(junc_factor_file, 'r') as handle:
        for line in handle:
            fields = re.split('[#!]', line)[0].split()
            if not line.strip():
                distances = None
            elif len(fields) == 4:
                resname, bond = fields[0], '-'.join(fields[1:3])
                distances = junctfactor.setdefault(resname, dict()).setdefault(bond, dict())
            elif distances is not None and fields:
                distances[int(fields[0])] = float(fields[1].split('d')[0])
    return junctfactor


def find_link_partners(model, region):
    link_partners = dict()
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    bond_params_table = model.restraints_manager.geometry.bond_params_table
    for link in region.link_atoms:
        min_distance = float('inf')
        closest = None
        link_atom = atoms[region.position(link)]
        for atom in atoms:
            serial = int(atom.serial.strip())
            if serial != link and bond_params_table.lookup(region.position(link), region.position(serial)) is not None:
                distance = atom.distance(link_atom)
                if distance < min_distance:
                    closest = serial
                    min_distance = distance
        if closest is None:
            raise SystemExit('No atom bonded to link atom ' + str(link) + ' in the QM system. Exiting...')
        link_partners[link] = closest
    return link_partners


def check_overlap(syst1_files):
    # E_MM counts the restraints among atoms that two QM regions share once, but
    # each region subtracts them again in E_MM1. Link atoms may be shared: there
    # are no restraints among them, and Orca only sees the hydrogens replacing them.
    regions = dict()
    qm_atoms = set()
    for syst1 in syst1_files:
        serials, link_atoms = read_region(syst1)
        shared = sorted(atom for atom in serials if atom in regions and (atom in qm_atoms or atom not in link_atoms))
        if shared:
            others = sorted(set(other for atom in shared for other in regions[atom]))
            raise SystemExit(syst1 + ' shares atoms ' + ', '.join(str(atom) for atom in shared) + ' with ' + ', '.join(others) + '. QM regions must not share atoms other than link atoms; duplicate the QM region into separate altlocs instead. Exiting...')
        for atom in serials:
            regions.setdefault(atom, []).append(syst1)
        qm_atoms |= serials - link_atoms


def junction_factors(model, region, junctfactor, ltype):
    g_bond = dict()
    atoms = model.get_hierarchy().atoms()
    for key, value in region.link_partners.items():
        c_qm = atoms[region.position(value)]
        link_atom = atoms[region.position(key)]
        resname = link_atom.parent().resname.strip()
        bond = '-'.join([c_qm.name.strip(), link_atom.name.strip()])
        if resname not in junctfactor:
            raise SystemExit(resname + ' missing in junction factor file. Exiting...')
        if bond not in junctfactor[resname]:
            raise SystemExit(bond + ' (between atoms ' + str(value) + ' and ' + str(key) + ') missing for ' + resname + ' in junction factor file. Exiting...')
        if ltype not in junctfactor[resname][bond]:
            raise SystemExit('Missing ltype ' + str(ltype) + ' for bond ' + bond + ' in residue ' + resname + ' in junction factor file. Exiting...')
        c_h_ideal = junctfactor[resname][bond][ltype] # eg. print(junctfactor['PHE']['CB-CG'][12])
        g_bond[key] = c_h_ideal/model.restraints_manager.geometry.bond_params_table.lookup(region.position(key), region.position(value)).distance_ideal
    return g_bond


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


def prepare_restart(infile, outfile):
    records = {'ATOM', 'HETATM', 'ANISOU', 'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3', 'TER', 'END'}
    with open(infile, 'r') as handle, open(outfile, 'w') as out:
        line = handle.readline()
        while line:
            if line[0:6].strip() in records:
                out.write(line)
            line = handle.readline()


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


def add_options(region, args):
    # the restraints and transforms given for this region on the command line
    for restraint in args.restraint_distance or []:
        if int(restraint[0]) == region.number:
            region.restraints_distance.append([int(restraint[1]), int(restraint[2]), float(restraint[3]), float(restraint[4])])
    for restraint in args.restraint_angle or []:
        if int(restraint[0]) == region.number:
            region.restraints_angle.append([int(restraint[1]), int(restraint[2]), int(restraint[3]), float(restraint[4]), float(restraint[5])])
    for transform in args.transform or []:
        if int(transform[0]) == region.number:
            region.transforms.append({'atoms': transform[1],
                                      'R': [[float(x) for x in transform[3*i+2:3*i+5]] for i in range(3)],
                                      't': [float(x) for x in transform[11:14]]})


def prepare_region(number, syst1, model_mm, args, junctfactor):
    """Write mm_i_c.pdb for one region and, unless -u, qm_i_h.pdb, and return the region."""
    print('----------')
    print(syst1.center(10))
    print('----------')
    print()
    region = Region(number, syst1)
    model_mm1 = region.select(model_mm)
    model_mm1.process(pdb_interpretation_params=model_mm1.get_current_pdb_interpretation_params(), make_restraints=True)
    region.restore_serials(model_mm1)
    print('Writing file:  ' + region.mm1_file)
    model_mm1.get_hierarchy().write_pdb_file(file_name=region.mm1_file, crystal_symmetry=model_mm.crystal_symmetry(), anisou=False)

    if not args.skip_h:
        region.link_partners = find_link_partners(model_mm1, region)
        region.g_bond = junction_factors(model_mm1, region, junctfactor, args.ltype)
        add_options(region, args)
        region.apply_transforms(model_mm1)
        print('Writing file:  ' + region.qm_file)
        region.write_qm_model(model_mm1)

    suggest_selection_string(model_mm1)
    check_altlocs(model_mm1, region.link_atoms)
    print('----------')
    return region


def main():
    args = parse_args()
    print('*** This is qref_prep.py ***')
    print()
    junctfactor = None
    if not args.skip_h:
        junctfactor = read_junctfactor(args.junctfactor)
        check_overlap(args.syst1)
    dm = DataManager()
    if args.cif is not None:
        for cif in args.cif:
            dm.process_restraint_file(cif)
    dm.process_model_file(args.pdb)
    model_mm = dm.get_model(filename=args.pdb)
    model_mm.add_crystal_symmetry_if_necessary()

    regions = [prepare_region(number, syst1, model_mm, args, junctfactor)
               for number, syst1 in enumerate(args.syst1, 1)]

    if not args.skip_h:
        settings = Settings(args.syst1, n_atoms=model_mm.get_number_of_atoms(), restart=args.restart,
                            w_qm=args.w_qm, cif=args.cif, orca_binary=locate_binary('orca'),
                            ltype=args.ltype)
        for region in regions:
            settings.add(region)
        print('----------')
        print()
        print('Writing file:  qref.dat')
        print()
        print('----------')
        settings.write()

    if args.restart is not None:
        print()
        print('Writing file:  ' + args.restart)
        print()
        print('----------')
        prepare_restart(args.pdb, args.restart)

    # cleanup
    temp_files = ['settings.pickle']
    for path in temp_files:
        if os.path.exists(path):
            os.remove(path)


if (__name__ == "__main__"):
    main()
