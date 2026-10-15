#!/usr/bin/env cctbx.python

from __future__ import absolute_import, division, print_function

import argparse
import os
import shutil

from iotbx.data_manager import DataManager

from utils import Settings


def modify_qm_input(path, region):
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
    if region.link_atoms:
        new.append('%geom\n  Constraints\n')
        for link_atom in sorted(region.link_atoms):
            new.append('    {C ' + str(region.position(link_atom)) + ' C}\n')
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
    settings = Settings.read()
    dm = DataManager()
    dm.process_model_file(model_file)
    model_mm = dm.get_model(filename=model_file)
    model_mm.add_crystal_symmetry_if_necessary()
    print('----------')
    print()
    for region in settings.regions():
        model_mm1 = region.select(model_mm)
        geo_opt_dir = 'geo_opt_qm_constrained_' + str(region.number)
        print('Generating: ' + geo_opt_dir)
        if not os.path.isdir(geo_opt_dir):
            os.mkdir(geo_opt_dir)
        region.apply_transforms(model_mm1)
        region.write_qm_model(model_mm1, geo_opt_dir + '/' + region.qm_file)
        shutil.copy(region.input_file, geo_opt_dir)
        shutil.copy(region.gbw_file, geo_opt_dir + '/qm.gbw')
        modify_qm_input(geo_opt_dir + '/' + region.input_file, region)
    print()
    print('----------')


if (__name__ == "__main__"):
    main()
