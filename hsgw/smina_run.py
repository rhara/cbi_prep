#!/usr/bin/env python

import argparse, sys
import subprocess as sp
from os import cpu_count

class FatalError(Exception): pass

def run(command, capture=None, **kwargs):
    print(command)
    if capture:
        try:
            output = sp.check_output(command, shell=True, stderr=sp.STDOUT, **kwargs).decode()
        except sp.CalledProcessError:
            print(f'Error: {command}', file=sys.stderr)
            raise FatalError(command)
        if isinstance(capture, str):
            open(capture, 'wt').write(output)
        return output
    else:
        ret = sp.call(command, shell=True, **kwargs)
        if ret == 0:
            return ret
        else:
            print(f'Error: {command}', file=sys.stderr)
            raise FatalError(command)

def prep_from_pdb(pdbid, ligand_name, dir=None):
    if dir is not None:
        pdb_fname = dir + '/' + pdbid + '.pdb.gz'
    else:
        pdb_fname = pdbid + '.pdb.gz'
    run(f'python extract_lig.py {pdb_fname} {ligand_name}')
    run(f'obabel {pdbid}_{ligand_name}.pdb -O {pdbid}_{ligand_name}.mol2', capture=True)
    run(f'python make_apo.py {pdb_fname} {pdbid}_{ligand_name}.pdb {pdbid}_apo.pdb')
    
def prep_protein(pdbid):
    print('prpe_protein')
    output = run(f'python write_tleaprc.py {pdbid} {pdbid}_apo.pdb {pdbid}_apo_H.pdb {pdbid}_apo_H_ref.mol2', capture=f'{pdbid}_tleaprc')
    run(f'tleap -s -f {pdbid}_tleaprc', capture=f'{pdbid}_tleaplog')
    run(f'obabel {pdbid}_apo_H.pdb -O {pdbid}_apo_H_nocharge.mol2', capture=f'{pdbid}_obabellog')
    run(f'python make_charged_protein.py {pdbid}_apo_H_nocharge.mol2 {pdbid}_apo_H_ref.mol2 {pdbid}_apo_H_charged.mol2')

def run_smina(pdbid, ligand_name, ncpu=8, num_modes=4, seed=0):
    boxpars = run(f'python center.py {pdbid}_{ligand_name}.pdb', capture=True).strip()
    receptor = f'{pdbid}_apo_H_charged.mol2'
    ligand = f'{pdbid}_{ligand_name}.mol2'
    docked = f'{pdbid}_{ligand_name}_redock.sdf'
    log = f'{pdbid}_smina.log'
    run(f'smina -r {receptor} -l {ligand} {boxpars} --cpu {ncpu} --num_modes {num_modes} --seed 0 -o {docked} --log {log}')

def get_rmsd(pdbid, ligand_name):
    ref_mol = f'{pdbid}_{ligand_name}.pdb'
    fit_mol = f'{pdbid}_{ligand_name}_redock.sdf'
    out_fname = f'{pdbid}_{ligand_name}_rmsd'
    output = run(f'python rmsd.py {ref_mol} {fit_mol}', capture=out_fname)
    print()
    print('RMSD vs Xray')
    print(output)

def main(args):
    try:
        ncpu = cpu_count()
        ncpu = max(ncpu-2, 1)
        pdbid = args.pdbid
        ligand_name = args.ligand_name
        dir = args.dir
        print('='*10, 'Smina self-docking driver', pdbid, ligand_name, '='*10)
        prep_from_pdb(pdbid, ligand_name, dir)
        prep_protein(pdbid)
        run_smina(pdbid, ligand_name, ncpu=ncpu)
        get_rmsd(pdbid, ligand_name)
    except:
        return -1
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdbid', type=str)
    parser.add_argument('ligand_name', type=str)
    parser.add_argument('--dir', '-d', type=str, required=False)
    args = parser.parse_args()
    sys.exit(main(args))
