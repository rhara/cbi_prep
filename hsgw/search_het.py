#!/usr/bin/env python

import argparse, gzip, os

def main(args):
    pdb_name = args.pdb_name
    openf = gzip.open if pdb_name.endswith('.gz') else open
    print(os.path.basename(pdb_name)[:4])

    hetatms = {}
    for line in openf(pdb_name, 'rt'):
        line = line.rstrip()
        try:
            pdb_tag = line[:6]
            pdb_ligname = line[17:20]
            pdb_chainId = line[21]
        except:
            continue
        if pdb_tag == 'HETATM':
            if (pdb_chainId, pdb_ligname) not in hetatms:
                hetatms[pdb_chainId, pdb_ligname] = 0
            hetatms[pdb_chainId, pdb_ligname] += 1
    info = []
    for (chainId, ligname), count in hetatms.items():
        if ligname == 'HOH':
            continue
        if ligname.startswith(' '):
            continue
        info.append((ligname, chainId, count))
    info.sort(key=lambda x: (-x[2], x[0], x[1]))
    for ligname, chainId, count in info:
        if 10 < count:
            print(ligname, chainId, count, '*')
        else:
            print(ligname, chainId, count)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_name', type=str)
    args = parser.parse_args()
    main(args)
