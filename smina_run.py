import argparse, os, sys, uuid
import multiprocessing as mp
from cbilib.io import read_molecule, write_molecule
import multiprocessing as mp

leap_in_template = open('leap_in.template', 'rt').read()

def worker(args):
    global leap_in_template
    print(args)
    idir, odir, pdbid, pdb_name, sdf_name = args
    uid = uuid.uuid4().hex
    tempdir = f'tmp/{uid}'
    os.makedirs(tempdir, exist_ok=True)
    tt_iname = f'{tempdir}/{pdbid}_apo.pdb'
    os.system(f'gunzip -c {pdb_name} > {tt_iname}')
    tt_mol2 = tempdir + '/' + pdbid + '_apo_ref.mol2'
    tt_pdb = tempdir + '/' + pdbid + '_apo_H.pdb'
    leap_in = leap_in_template.replace('{iname}', tt_iname).\
                               replace('{oname_mol2}', tt_mol2).\
                               replace('{oname_pdb}', tt_pdb).\
                               replace('{pdbid}', pdbid)
    open(tempdir + '/leap_in', 'wt').write(leap_in)
    os.system(f'tleap -s -f {tempdir}/leap_in')
    os.system(f'python make_charged_protein.py {tempdir}/{pdbid}_apo_H.pdb {tempdir}/{pdbid}_apo_ref.mol2')
    os.system(f'obabel {tempdir}/{pdbid}_apo_H_charged.mol2 -O {tempdir}/{pdbid}_apo_H_charged.pdbqt')
    os.system(f'obabel {sdf_name} -O {tempdir}/{pdbid}_ligand.pdbqt')
    os.system(f'python get_center.py {sdf_name} > {tempdir}/center')
    x, y, z = eval(open(f'{tempdir}/center', 'rt').read())
    x = round(x, 3)
    y = round(y, 3)
    z = round(z, 3)
    ncpu = mp.cpu_count()
    ncpu -= 2
    os.system(f'smina -r {tempdir}/{pdbid}_apo_H_charged.pdbqt -l {tempdir}/{pdbid}_ligand.pdbqt --center_x {x} --center_y {y} --center_z {z} --size_x 25 --size_y 25 --size_z 25 --cpu {ncpu} --num_modes 4 -o {tempdir}/{pdbid}_redock.sdf')
    os.system(f'python rmsd_norm.py {sdf_name} {tempdir}/{pdbid}_redock.sdf > {tempdir}/rmsd')

    return idir, odir, pdbid, pdb_name, sdf_name

def main(args):
    odir = args.odir
    idir = args.idir

    def gen():
        count = 0
        for root, dirs, files in os.walk(idir):
            pdbid = root[-4:]
            sdf_name = None
            pdb_name = None
            for f in files:
                if f.endswith('.sdf'):
                    sdf_name = f
                if f.endswith('_apo.pdb.gz'):
                    pdb_name = f
            if sdf_name and pdb_name:
                count += 1
                pdb_name = root + '/' + pdb_name
                sdf_name = root + '/' + sdf_name
                yield idir, odir, pdbid, pdb_name, sdf_name
                # if count == 5:
                #     break

    pool = mp.Pool(mp.cpu_count())
    for idir, odir, pdbid, pdb_name, sdf_name in pool.imap_unordered(worker, gen()):
        print('Done:', idir, odir, pdbid, pdb_name, sdf_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('idir', type=str)
    parser.add_argument('odir', type=str)
    args = parser.parse_args()
    main(args)
