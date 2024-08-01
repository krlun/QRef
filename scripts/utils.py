import re
import json


def read_qm_and_link_atoms(infile):
    # reads QM and link atoms from syst1 file, assuming that the numbering in the syst1 file corresponds to serial in pdb for the whole protein
    # lines can be commented out with # and !, multiple atoms and intervals can be specified per line (separated by either , or blank)
    # first occurance of atom will be treated as part of qm system, second occurance as link atom
    qm_atoms = set()
    link_atoms = set()
    comments = '[#!]'
    delimiters = '[^,\s]+'
    with open(infile, 'r') as file:
        line = file.readline()
        while line:
            line = re.findall(delimiters, re.split(comments, line)[0])
            for interval in line:
                interval = interval.split('-')
                for i in range(min([int(x) for x in interval]), max([int(x) for x in interval]) + 1):
                    link_atoms.add(i) if i in qm_atoms else qm_atoms.add(i)
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