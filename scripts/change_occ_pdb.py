#!/opt/homebrew/Caskroom/miniconda/base/bin/python

import sys
import argparse
from utils import read_syst1


def parse_and_update(occupancy, infile, outfile, atoms=None, residues=None):
    records = {'ATOM', 'HETATM'}
    with open(infile, 'r') as file:
        with open(outfile, 'w') as out:
            line = file.readline()
            while line:
                if line[0:6].strip() in records:
                    if (atoms is not None and int(line[6:11].strip()) in atoms) or (residues is not None and line[16:26] in residues):
                        line = line[0:54] + occupancy.rjust(6) + line[60:]
                out.write(line)
                line = file.readline()

def read_residues(infile):
    residues = set()
    with open(infile, 'r') as file:
        line = file.readline().rstrip()
        while line:
            residues.add(line)
            line = file.readline().rstrip()
    return residues


def parse_args():
    parser = argparse.ArgumentParser(
        description='Rudimentary script that changes occupancies in a PDB file based on either a syst1 file or a list of residues (columns 17-26 in the PDB file are used as identifiers).'
    )
    parser.add_argument('pdb', type=str, help='PDB file')
    parser.add_argument('occupancy', type=float, help='desired occupancy')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--syst1', type=str, help='name of the syst1 file')
    group.add_argument('-r', '--residues', type=str, help='name of the residue file')
    parser.add_argument('-v', '--version', action="version", version="%(prog)s 0.0.1")
    return parser.parse_args()

def main():
    args = parse_args(args)
    occupancy = '{:.2f}'.format(round(float(args.occupancy), 2))
    atoms = None
    residues = None
    if args.syst1 is not None:
        atoms, _ = read_syst1(args.syst1)
        outfile = args.pdb[:-4] + '_' + args.syst1 + '_' + occupancy + args.pdb[-4:]
    elif args.residues is not None:
        residues = read_residues(args.residues)
        outfile = args.pdb[:-4] + '_' + args.residues + '_' + occupancy + args.pdb[-4:]
    parse_and_update(occupancy, args.pdb, outfile, atoms=atoms, residues=residues)

if __name__ == "__main__":
    main()