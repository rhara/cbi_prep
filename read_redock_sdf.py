from rdkit import Chem
from collections import namedtuple
import sys, os, argparse
import numpy as np

ELEMENTS = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53}

def get_lines_from_sdf(fname, no=1):
    Line = namedtuple('Line', ['no', 'section', 'linecount', 'cont'])
    Atom = namedtuple('Atom', ['x', 'y', 'z', 'el', 'hyb', 'gr', 'gr_name', 'charge'])
    Bond = namedtuple('Bond', ['idx', 'atom1', 'atom2', 'degree'])
    f = open(fname)
    molno = 1
    count = 0
    section = 'HEADER'
    lines = []
    atomcount = 0
    bondcount = 0
    for line in f:
        line = line.strip()
        if line == '$$$$':
            molno += 1
            count = 0
            section = 'HEADER'
            continue
        count += 1
        if count < 4:
            section = 'HEADER'
            lines.append(Line(molno, section, count, line))
            continue
        if count == 4:
            section = 'COUNT'
            it = line.split()
            natoms, nbonds = int(it[0]), int(it[1])
            lines.append(Line(molno, section, count, line))
            section = 'ATOM'
            atomcount = 0
            continue
        if section == 'ATOM':
            atomcount += 1
            lines.append(Line(molno, section, count, line))
            if atomcount == natoms:
                section = 'BOND'
                bondcount = 0
            continue
        if section == 'BOND':
            bondcount += 1
            lines.append(Line(molno, section, count, line))
            if bondcount == nbonds:
                section = 'END'
                bondcount = 0
            continue
        if section == 'END':
            lines.append(Line(molno, section, count, line))
            section = 'TAG'
            continue
        if section == 'TAG':
            lines.append(Line(molno, section, count, line))
            continue
        section = 'OTHER'
        lines.append(Line(molno, section, count, line))

    for line in lines:
        if line.no != no:
            continue
        if line.section == 'ATOM':
            print(line)
    

def main(args):
    iname = args.iname
    get_lines_from_sdf(iname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('iname', type=str)
    args = parser.parse_args()
    main(args)
