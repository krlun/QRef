import re
import json


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


def read_qm_and_link_atoms(infile):
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