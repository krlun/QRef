#!/Applications/phenix-1.20.1-4487/build/bin/cctbx.python

import os
import sys
import argparse
import json
import subprocess

from iotbx.data_manager import DataManager
from cctbx.array_family import flex

from utils import read_qm_and_link_atoms
from utils import read_junc_factors


def identify_link_pairs(model, link_atoms, serial_to_index):
    link_pairs = dict()
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    for link in link_atoms:
        min_distance = float('inf')
        link_atom = atoms[serial_to_index[link]]
        for atom in atoms:
            if atom != link_atom and atom.element.strip() == 'C':
                distance = atom.distance(link_atom)
                if distance < min_distance:
                    closest = int(atom.serial.strip())
                    min_distance = distance
        link_pairs[link] = closest
    return link_pairs


def calculate_g_factor(model, link_pairs, junc_factors, ltype, serial_to_index):
    g = dict()
    atoms = model.get_hierarchy().atoms()
    for key, value in link_pairs.items():
        c_qm = atoms[serial_to_index[value]]
        link_atom = atoms[serial_to_index[key]]
        resname = link_atom.parent().resname
        bond = '-'.join([c_qm.name.strip(), link_atom.name.strip()])
        if resname not in junc_factors.keys():
            raise SystemExit(resname + ' missing in junction factor file. Exiting...')
        if bond not in junc_factors[resname].keys():
            raise SystemExit(bond + ' (between atoms ' + str(value) + ' and ' + str(key) + ') missing for ' + resname + ' in junction factor file. Exiting...')
        if ltype not in junc_factors[resname][bond].keys():
            raise SystemExit('Missing ltype ' +  ltype + ' for bond ' + bond + ' in residue ' + resname + ' in junction factor file. Exiting...')
        c_h_ideal = junc_factors[resname][bond][ltype] # eg. print(junc_factors['PHE']['CB-CG'][12])
        g[key] = c_h_ideal/model.restraints_manager.geometry.bond_params_table.lookup(serial_to_index[key], serial_to_index[value]).distance_ideal
    return g


def convert_serial_to_index(qm):
    qm_sorted = sorted(qm)
    indices = dict()
    for i in range(len(qm_sorted)):
        indices[qm_sorted[i]] = i
    return indices


def restore_serial_in_model(model, serial_to_index):
    index_to_serial = {value: key for key, value in serial_to_index.items()}
    atoms = model.get_hierarchy().atoms()
    width = 5
    for atom in atoms: atom.serial = str(index_to_serial[int(atom.serial) - 1]).rjust(width)


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


def suggest_selection_string(model):
    hierarchy = model.get_hierarchy()
    selection_string = '('
    for chain in hierarchy.only_model().chains():
        selection_string += '('
        for residue in chain.residue_groups():
            selection_string += 'resseq ' + residue.resseq.strip() + ' or '
        selection_string = selection_string[:-4] + ') and chain ' + chain.id + ') or ('
    selection_string = selection_string[:-5]
    print
    print('Suggested selection string (reciprocal space):')
    print
    print('\"' + selection_string + '\"')
    print
    print('Suggested selection string (real space):')
    print
    print('\"all and not (' + selection_string + ')\"')
    print


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
    return binary if len(out) < len(binary) else out.rstrip()


def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Prepares for quantum refinement by extracting the QM system(s) as well as saving a settings file (qref.dat).'
    )
    parser.add_argument('pdb', type=str, help='PDB file describing the entire model.')
    parser.add_argument('-c', '--cif', nargs='+', type=str, help='CIF files describing the restraints for novel ligands')
    parser.add_argument('-s', '--syst1', nargs='+', default=['syst1'], type=str, help='name of the file(s) describing what should constitute the QM system(s) (default: \'syst1\')')
    parser.add_argument('-u', '--skip_h', action='store_true')
    parser.add_argument('-j', '--junctfactor', nargs=1, default='junctfactor', type=str, help='name of file containing the junction factors (default: \'junctfactor\')')
    parser.add_argument('-l', '--ltype', nargs=1, default=12, type=int, help='link type (default: 12)')
    parser.add_argument('-w', '--w_qm', nargs=1, default=7.5, type=float, help='scaling factor for QM energy and gradients (default: 7.5)')
    parser.add_argument('-r', '--restart', nargs='?', const='restart.pdb', type=str, help='if specified creates a restart file (optional: name)')
    parser.add_argument('-rd', '--restraint_distance', nargs=5, action='append', type=str, help='if specified applies a harmonic (bond) distance restraint according to \'i atom1_serial atom2_serial desired_distance force_constant\'')
    parser.add_argument('-ra', '--restraint_angle', nargs=6, action='append', type=str, help='if specified applies a harmonic (bond) angle restraint according to \'i atom1_serial atom2_serial atom3_serial desired_angle force_constant\'')
    #parser.add_argument('-t', '--transform', nargs=4, action='append', type=str, help='if specified applies the specified transformations according to \'i atoms R t\'')
    parser.add_argument('-v', '--version', action="version", version="%(prog)s 0.2.0")
    return parser.parse_args(args)


def main(args):
    print('*** This is qref_prep.py ***')
    print
    args = parse_args(args)

    model_real_file = args.pdb

    syst1_files = args.syst1
    junc_factor_file = args.junctfactor
    ltype = args.ltype
    if args.skip_h is not True:
        junc_factors = read_junc_factors(junc_factor_file=junc_factor_file)
    dm = DataManager()
    if args.cif is not None:
        for cif in args.cif:
            dm.process_restraint_file(cif)
    dm.process_model_file(model_real_file)
    model_real = dm.get_model(filename=model_real_file)
    model_real.add_crystal_symmetry_if_necessary()

    dat = dict()

    for index, syst1 in enumerate(syst1_files, 1):
        print('----------')
        print(syst1.center(10))
        print('----------')
        print
        qm_atoms, link_atoms = read_qm_and_link_atoms(syst1)
        serial_to_index = convert_serial_to_index(qm_atoms)
        model_model = select_qm_model(model=model_real, qm=qm_atoms)
        model_model.process(pdb_interpretation_params=model_model.get_current_pdb_interpretation_params(), make_restraints=True)
        restore_serial_in_model(model_model, serial_to_index=serial_to_index)
        name_c = 'qm_' + str(index) + '_c.pdb'
        print('Writing file:  ' + name_c)
        model_model.get_hierarchy().write_pdb_file(file_name=name_c, crystal_symmetry=model_real.crystal_symmetry())

        if args.skip_h is not True:
            dat[syst1] = dict()
            link_pairs = identify_link_pairs(model_model, link_atoms, serial_to_index)
            dat[syst1]['link_pairs'] = link_pairs
            g = calculate_g_factor(model_model, link_pairs=link_pairs, junc_factors=junc_factors, ltype=ltype, serial_to_index=serial_to_index)
            dat[syst1]['g'] = g
            dat[syst1]['restraint_distance'] = list()
            dat[syst1]['restraint_angle'] = list()
            #dat[syst1]['transform'] = list() # to be implemented
            if args.restraint_distance is not None:
                for restraint in args.restraint_distance:
                    if int(restraint[0]) == index:
                        dat[syst1]['restraint_distance'].append([int(restraint[1]), int(restraint[2]), float(restraint[3]), float(restraint[4])])
            if args.restraint_angle is not None:
                for restraint in args.restraint_angle:
                    if int(restraint[0]) == index:
                        dat[syst1]['restraint_angle'].append([int(restraint[1]), int(restraint[2]), int(restraint[3]), float(restraint[4]), float(restraint[5])])
            # to be implemented
            # if args.transform is not None:
            #     for transform in args.transform:
            #         if int(transform[0]) == index:
            #             dat[syst1]['transform'].append(dict())
            #             dat[syst1]['transform'][-1]['atoms'] = transform[1]
            #             dat[syst1]['transform'][-1]['R'] = transform[2]
            #             dat[syst1]['transform'][-1]['t'] = transform[3]
            name_h = 'qm_' + str(index) + '_h.pdb'
            print('Writing file:  ' + name_h)
            write_pdb_h(name_h, model_model, link_pairs, g, serial_to_index)


        suggest_selection_string(model_model)

    print('----------')
    print

    if args.skip_h is not True:
        dat['n_atoms'] = model_real.get_number_of_atoms()
        dat['restart'] = args.restart
        dat['w_qm'] = args.w_qm
        dat['cif'] = args.cif
        dat['orca_binary'] = locate_binary('orca')
        dat['ltype'] = ltype
        dat['syst1_files'] = args.syst1
        print('Writing file:  qref.dat')
        write_dat(dat)

    if args.restart is not None:
        print('Writing file:  ' + args.restart)
        prepare_restart(model_real_file, args.restart)

    print
    print('----------')

    # cleanup
    temp_files = ['settings.pickle']
    for file in temp_files:
        if os.path.exists(file):
            os.remove(file)


if (__name__ == "__main__"):
    main(sys.argv[1:])
