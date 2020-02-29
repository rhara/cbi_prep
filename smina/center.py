#!/usr/bin/env python

import argparse
from rdkit import Chem


def main(args):
    iname = args.iname
    coords = []
    mol = Chem.MolFromMolFile(iname, sanitize=False)
    conf = mol.GetConformer(0)
    coords = conf.GetPositions()
    center = [round(x, 3) for x in coords.mean(axis=0)]
    x, y, z = center
    print(f'--center_x {x} --center_y {y} --center_z {z} --size_x 25 --size_y 25 --size_z 25')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('iname', type=str)
    args = parser.parse_args()
    main(args)
