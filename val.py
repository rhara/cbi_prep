#!/usr/bin/env python

import sys, os, re, argparse
sys.path.insert(0, os.getcwd())
from cbilib.catalog import get_catalog
from math import log10

def main(args):
    pat_Kd = re.compile('([-0-9.]+)(mM|uM|nM|pM|fM)')
    factor = {'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}

    csv_name = args.csv
    odir = args.out

    df = get_catalog(csv_name)
    D = {}
    for i in df.index:
        pdbid, Kd = df.loc[i, ['PDB_code', 'Kd']]
        m = pat_Kd.match(Kd)
        D[pdbid] = -log10(float(m.group(1))*factor[m.group(2)])
    for pdbid in D:
        if os.path.exists(f'{odir}/{pdbid}'):
            v = round(D[pdbid], 3)
            print(pdbid, v)
            print(v, file=open(f'{odir}/{pdbid}/value', 'wt'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', '-c', type=str, required=True)
    parser.add_argument('--out', '-o', type=str, required=True)
    args = parser.parse_args()
    main(args)
