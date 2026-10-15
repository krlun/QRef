#!/usr/bin/env cctbx.python

from __future__ import absolute_import, division, print_function

import argparse
import os
import shutil

from iotbx.data_manager import DataManager

from utils import apply_transforms
from utils import read_settings
from utils import read_region
from utils import select_region
from utils import write_qm_model


def modify_qm_input(path, link_atoms, serial_to_index):
    # To be implemented
    # theory = {6: '! B3LYP 6-31G(d) D4 UKS Opt',
    #           7: '! BP ',
    #           12: '! TPSS DEF2-SV(P) D4 UKS Opt'
    # }

    # 5: experimental (Amber)
    # 6: B3LYP/6-31G*
    # 7: BP(RI)/6-31G*
    # 8: BP(RI)/SVP
    # 9: BP(RI)/def2-SV(P)
    #10: PBE(RI)/def2-SVP 
    #11: B3LYP(RI)/def2-SV(P) 
    #12: TPSS(RI)/def2-SV(P) 
    #13: B97-D(RI)/def2-SV(P) 

    with open(path, 'r') as f:
        old = f.readlines()
    new = list()
    for line in old:
        if 'engrad' not in line.lower():
            new.append(line)
    new.append('! MORead\n%moinp "qm.gbw"\n')
    new.append('! Opt\n')
    if link_atoms:
        new.append('%geom\n  Constraints\n')
        for link_atom in link_atoms:
            new.append('    {C ' + str(serial_to_index[link_atom]) + ' C}\n')
        new.append('  end\nend\n')
    with open(path, 'w') as f:
        f.writelines(new)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Prepares for geometry optimization after a QRef run, with fixed junction atoms.'
    )
    parser.add_argument('pdb', type=str, help='PDB file describing the entire refined model.')
    return parser.parse_args()


def main():
    if not os.path.isfile('qref.dat'):
        raise SystemExit('\'qref.dat\' missing, exiting...')
    args = parse_args()
    model_file = args.pdb
    settings = read_settings('qref.dat')
    dm = DataManager()
    dm.process_model_file(model_file)
    model_mm = dm.get_model(filename=model_file)
    model_mm.add_crystal_symmetry_if_necessary()
    print('----------')
    print()
    for index, syst1 in enumerate(settings['syst1_files'], 1):
        serials, link_atoms = read_region(syst1)
        serial_to_index = {value: key for key, value in dict(enumerate(sorted(serials))).items()}
        model_mm1 = select_region(model=model_mm, serials=serials)
        geo_opt_dir = 'geo_opt_qm_constrained_' + str(index)
        print('Generating: ' + geo_opt_dir)
        qm_pdb_file = 'qm_' + str(index) + '_h.pdb'
        qm_inp_file = 'qm_' + str(index) + '.inp'
        if not os.path.isdir(geo_opt_dir):
            os.mkdir(geo_opt_dir)
        if 'transforms' in settings[syst1]: # backwards compatibility
            apply_transforms(model_mm1, settings[syst1]['transforms'], serial_to_index)
        write_qm_model(geo_opt_dir + '/' + qm_pdb_file, model_mm1, settings[syst1]['link_pairs'], settings[syst1]['g'], serial_to_index)
        shutil.copy(qm_inp_file, geo_opt_dir)
        shutil.copy('qm_' + str(index) + '.gbw', geo_opt_dir + '/qm.gbw')
        modify_qm_input(geo_opt_dir + '/' + qm_inp_file, sorted(link_atoms), serial_to_index)
    print()
    print('----------')


if (__name__ == "__main__"):
    main()
