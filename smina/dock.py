#!/usr/bin/env python

import os, gzip, uuid, argparse, shutil
import openbabel as ob
import subprocess as sp

def log(cont):
    textfs.write(cont)

def prep_protein(uid, protein_iname):
    tmpdir = f'tmp/{uid}.tmpdir'
    if protein_iname.endswith('.gz'):
        openf = gzip.open
    else:
        openf = open
    protein_cont = openf(protein_iname, 'rt').read()
    protein_oname = f'{tmpdir}/protein.pdb'
    open(protein_oname, 'wt').write(protein_cont)

def prep_ligand(uid, ligand_iname):
    tmpdir = f'tmp/{uid}.tmpdir'
    ligand_cont = open(ligand_iname, 'rt').read()
    ligand_oname = f'{tmpdir}/ligand.sdf'
    open(ligand_oname, 'wt').write(ligand_cont)

def tleap(uid):
    tmpdir = f'tmp/{uid}.tmpdir'
    cont = open(f'{THISDIR}/tleaprc.template', 'rt').read()
    cont = cont.replace('{{tmpdir}}', tmpdir)
    open(f'{tmpdir}/tleaprc', 'wt').write(cont)
    sp.call(f'tleap -s -f {tmpdir}/tleaprc', shell=True)
    sp.call(f'obabel {tmpdir}/protein.mol2 -O {tmpdir}/protein_charged.mol2 > /dev/null 2>&1', shell=True)

def smina(uid, ncpu=None, num_modes=4, seed=0):
    tmpdir = f'tmp/{uid}.tmpdir'
    if ncpu is None:
        ncpu = os.cpu_count()
    receptor = f'{tmpdir}/protein_charged.mol2'
    ligand = f'{tmpdir}/ligand.sdf'
    oname = f'{tmpdir}/docked.sdf'
    logname = f'{tmpdir}/smina.log'
    boxpars = sp.check_output(f'python {THISDIR}/center.py {ligand}', shell=True).decode().strip()
    sp.call(f'smina -r {receptor} -l {ligand} {boxpars} --cpu {ncpu} --num_modes {num_modes} --seed {seed} -o {oname} --log {logname}', shell=True)

def rmsd(uid):
    tmpdir = f'tmp/{uid}.tmpdir'
    ref = f'{tmpdir}/ligand.sdf'
    fit = f'{tmpdir}/docked.sdf'
    oname = f'{tmpdir}/rmsd'
    cont = sp.check_output(f'python {THISDIR}/rmsd.py {ref} {fit}', shell=True).decode().strip()
    print(cont)
    open(oname, 'wt').write(cont + '\n')

def retrieve(uid):
    tmpdir = f'tmp/{uid}.tmpdir'
    sp.call(f'cp -av {tmpdir}/smina.log ./', shell=True)
    sp.call(f'cp -av {tmpdir}/docked.sdf ./', shell=True)
    sp.call(f'cp -av {tmpdir}/tleap.log ./', shell=True)
    sp.call(f'cp -av {tmpdir}/rmsd ./', shell=True)


THISDIR = None
textfs = None

def main(args):
    global THISDIR, textfs

    protein_iname = args.protein_iname
    ligand_iname = args.ligand_iname

    assert protein_iname.endswith('.pdb') or protein_iname.endswith('.pdb.gz')
    assert ligand_iname.endswith('.sdf')

    uid = uuid.uuid4().hex
    print(uid)

    THISDIR = os.path.abspath(os.path.dirname(__file__))

    tmpdir = f'tmp/{uid}.tmpdir'
    os.makedirs(tmpdir, mode=0o775, exist_ok=True)
    textfs = open(f'{tmpdir}/info.txt', 'wt')

    log(f'uid={uid}\n')
    log(f'protein_iname={protein_iname}\n')
    log(f'ligand_iname={ligand_iname}\n')

    prep_protein(uid, protein_iname)
    prep_ligand(uid, ligand_iname)
    tleap(uid)
    smina(uid)
    rmsd(uid)
    retrieve(uid)
    # shutil.rmtree(uid)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protein_iname', type=str)
    parser.add_argument('ligand_iname', type=str)
    args = parser.parse_args()
    main(args)
