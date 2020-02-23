#!/usr/bin/env python

import argparse
import gzip
import numpy as np
import math

"""
0         1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890123456789
ATOM    133  CA  LYS A  43      15.845  21.095  -4.234  1.00 49.99           C
"""
def main(args):
    pdb_name = args.pdb_name
    ligand_name = args.ligand_name
    oname = args.oname

    openf = gzip.open if pdb_name.endswith('.gz') else open

    coords = []
    for line in open(ligand_name, 'rt'):
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append([x, y, z])
    ligand_coords = np.array(coords)

    chains = {}
    for line in openf(pdb_name, 'rt'):
        line = line.rstrip()
        try:
            pdb_tag = line[:6]
            pdb_chainId = line[21]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except:
            continue
        if pdb_tag == 'ATOM  ':
            if pdb_chainId not in chains:
                chains[pdb_chainId] = []
            chains[pdb_chainId].append([x, y, z])
    dist = {}
    for k in chains:
        chains[k] = np.array(chains[k])
        ds = []
        for i in range(len(chains[k])):
            for j in range(len(ligand_coords)):
                v = chains[k][i] - ligand_coords[j]
                ds.append(sum(v*v))
        dist[k] = math.sqrt(min(ds))
    for k in list(dist.keys()):
        if 6 < dist[k]:
            del dist[k]

    chains = list(dist)

    openf2 = gzip.open if oname.endswith('.gz') else open

    with openf2(oname, 'wt') as out:
        for line in openf(pdb_name, 'rt'):
            line = line.rstrip()
            try:
                pdb_tag = line[:6]
                pdb_chainId = line[21]
            except:
                continue
            if pdb_tag == 'ATOM  ' and pdb_chainId in chains:
                out.write(line + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_name', type=str)
    parser.add_argument('ligand_name', type=str)
    parser.add_argument('oname', type=str)
    args = parser.parse_args()
    main(args)

