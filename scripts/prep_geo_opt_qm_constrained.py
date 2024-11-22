#!/Applications/phenix-1.21.2-5419/build/bin/cctbx.python

import argparse
import os
import shutil

from iotbx.data_manager import DataManager

from utils import apply_transforms
from utils import read_dat
from utils import read_syst1
from utils import select_qm_model
from utils import write_pdb_h


def modify_qm_input(file, link_atoms, serial_to_index):
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

    with open(file, 'r') as f:
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
    with open(file, 'w') as f:
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
    dat = read_dat('qref.dat')
    dm = DataManager()
    dm.process_model_file(model_file)
    model_mm = dm.get_model(filename=model_file)
    print('----------')
    print
    for index, syst1 in enumerate(dat['syst1_files'], 1):
        qm_atoms, link_atoms = read_syst1(syst1)
        serial_to_index = {value: key for key, value in dict(enumerate(sorted(qm_atoms))).items()}
        model_mm1 = select_qm_model(model=model_mm, qm=qm_atoms)
        geo_opt_dir = 'geo_opt_qm_constrained_' + str(index)
        print('Generating: ' + geo_opt_dir)
        qm_pdb_file = 'qm_' + str(index) + '_h.pdb'
        qm_inp_file = 'qm_' + str(index) + '.inp'
        if not os.path.isdir(geo_opt_dir):
            os.mkdir(geo_opt_dir)
        if 'transforms' in dat[syst1].keys(): # backwards compatability
            apply_transforms(model_mm1, dat[syst1]['transforms'], serial_to_index)
        write_pdb_h(geo_opt_dir + '/' + qm_pdb_file, model_mm1, dat[syst1]['link_pairs'], dat[syst1]['g'], serial_to_index)
        shutil.copy(qm_inp_file, geo_opt_dir)
        shutil.copy('qm_' + str(index) + '.gbw', geo_opt_dir + '/qm.gbw')
        modify_qm_input(geo_opt_dir + '/' + qm_inp_file, sorted(link_atoms), serial_to_index)
    print
    print('----------')


if (__name__ == "__main__"):
    main()
