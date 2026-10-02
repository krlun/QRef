#!/usr/bin/env python3
"""Stands in for Orca, writing a .engrad for an analytic energy surface."""
import math
import os
import sys


def pdb_from_input(path):
    with open(path) as handle:
        for line in handle:
            if line.strip().startswith('*pdbfile'):
                return line.split()[-1]
    raise SystemExit(f'stub_orca: no *pdbfile line in {path}')


def coordinates(path):
    sites = []
    with open(path) as handle:
        for line in handle:
            if line[0:6].strip() in ('ATOM', 'HETATM'):
                sites.append((float(line[30:38]), float(line[38:46]),
                              float(line[46:54])))
    return sites


BOHR = 0.529177249              # Angstrom per bohr, as QRef uses


# gradient() must stay the derivative of energy(), in Eh/bohr, or
# test_gradient_matches_finite_difference measures the stub's error, not QRef's
def energy(sites):
    return -sum(math.sin(x) + math.cos(y) + math.sin(z) for x, y, z in sites)


def gradient(sites):
    return [(-math.cos(x) * BOHR, math.sin(y) * BOHR, -math.cos(z) * BOHR)
            for x, y, z in sites]


def write_engrad(path, sites):
    # read_energy_and_gradient_from_orca() matches these header lines whole, so
    # they have to read exactly as Orca writes them
    with open(path, 'w') as out:
        out.write(f'#\n# Number of atoms\n#\n {len(sites)}\n')
        out.write('#\n# The current total energy in Eh\n#\n')
        out.write(f'{energy(sites):22.12f}\n')
        out.write('#\n# The current gradient in Eh/bohr\n#\n')
        for components in gradient(sites):
            for value in components:
                out.write(f'{value:22.12f}\n')


def main(argv):
    if len(argv) != 1:
        raise SystemExit('usage: stub_orca.py <input>')
    stem = os.path.splitext(argv[0])[0]
    sites = coordinates(pdb_from_input(argv[0]))
    write_engrad(stem + '.engrad', sites)
    print(f'stub_orca: {len(sites)} atoms')
    # qref's logging() reads this out of qm_i.out to tell success from failure
    print('ORCA TERMINATED NORMALLY')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
