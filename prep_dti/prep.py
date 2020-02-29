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
    csv_name = args.csv_name
    odir = args.odir
    df = get_catalog(csv_name)
    print(df)

    os.makedirs(odir, exist_ok=True)
    os.chdir(odir)

    def gen():
        count = 0
        for i in df.index:
            count += 1
            pdbid = df.loc[i, 'pdbid']
            yield count, pdbid

    pool = mp.Pool(mp.cpu_count())
    for count, pdbid in pool.imap_unordered(worker, gen()):
        print(f'{count:5} {pdbid}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_name', type=str)
    parser.add_argument('odir', type=str)
    args = parser.parse_args()
    main(args)
