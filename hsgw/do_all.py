#!/usr/bin/env python

import argparse
from collections import namedtuple
import subprocess as sp
import pandas as pd

Record = namedtuple('Record', ['pdbid', 'ligand_name'])
AMINOACIDS = 'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()
EXCLUDES = 'NH4'.split()

def get_catalog(csv_name):
    df = pd.read_csv(csv_name)
    records = []
    for i in df.index:
        r = df.iloc[i]
        pdbid = r.loc['PDB_code']
        ligand_name = r.loc['ligand_name'][1:-1].strip()
        if ligand_name in AMINOACIDS:
            continue
        if  ligand_name in EXCLUDES:
            continue
        if len(ligand_name) != 3:
            continue
        if ligand_name.startswith('_'):
            ligand_name = ligand_name[1:]
        records.append(Record(pdbid, ligand_name))
    return {r.pdbid: r.ligand_name for r in records}

def main(args):
    csv_name = args.csv_name
    dir = args.dir

    cat = get_catalog(csv_name)
    command = 'python smina_run.py'
    if dir:
        command += f' --dir {dir}'
    for pdbid, ligand_name in cat.items():
            sp.call(f'{command} {pdbid} {ligand_name}', shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_name', type=str)
    parser.add_argument('--dir', '-d', type=str, required=False)
    args = parser.parse_args()
    main(args)

