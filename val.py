#!/usr/bin/env python

import sys, os, re, argparse
sys.path.insert(0, os.getcwd())
from cbilib.catalog import get_catalog

def main(args):
    pat_Kd = re.compile('([-0-9.]+)(mM|uM|nM|pM)')
    factor = {'mM': 1000000, 'uM': 1000, 'nM': 1, 'pM': 0.001}

    csv_name = args.csv
    odir = args.out

    df = get_catalog(csv_name)
    D = {}
    for i in df.index:
        pdbid, Kd = df.loc[i, ['PDB_code', 'Kd']]
        m = pat_Kd.match(Kd)
        D[pdbid] = float(m.group(1))*factor[m.group(2)]
    for pdbid in D:
        if os.path.exists(f'{odir}/{pdbid}'):
            v = D[pdbid]
            print(pdbid, v)
            print(v, file=open(f'{odir}/{pdbid}/value', 'wt'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', '-c', type=str, required=True)
    parser.add_argument('--out', '-o', type=str, required=True)
    args = parser.parse_args()
    main(args)
