#!/usr/bin/env python

import os, sys, argparse
sys.path.insert(0, os.getcwd())
import prody
from cbilib.catalog import get_catalog
import multiprocessing as mp

def worker(args):
    count, pdbid = args
    prody.fetchPDB(pdbid)
    return count, pdbid

def main(args):
    csv_name = args.csv
    out_name = args.out
    df = get_catalog(csv_name)
    print(df)

    os.makedirs(out_name, exist_ok=True)
    os.chdir(out_name)

    def gen():
        count = 0
        for i in df.index:
            count += 1
            pdbid = df.loc[i, 'PDB_code']
            yield count, pdbid

    pool = mp.Pool(mp.cpu_count())
    for count, pdbid in pool.imap_unordered(worker, gen()):
        print(count, pdbid)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', '-i', type=str, required=True)
    parser.add_argument('--out', '-o', type=str, required=True)
    args = parser.parse_args()
    main(args)
