import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import math

"""
0         1         2         3         4         5         6         7
01234567890123456789012345678901234567890123456789012345678901234567890123456789
HETATM10147  C7  PNM A1007      25.104  27.308  86.289  1.00 42.92           C
"""

ATOMS = dict(H=1, C=6, N=7, O=8, F=9, P=15, S=16, CL=17, BR=35, I=53)

def mk_mols_atom_only(atoms):
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for i, (el, x, y, z) in enumerate(atoms):
        atom = Chem.Atom(el)
        mol.AddAtom(atom)
        conf.SetAtomPosition(i, [x, y, z])
    mol.AddConformer(conf)
    return Chem.Mol(mol)

def read_pdb(fname):
    atoms = []
    for line in open(fname, 'rt'):
        line = line.strip()
        if not line.startswith('HETATM'):
            continue
        el = line[76:78].strip().upper()
        if el not in ATOMS:
            raise Exception(f'{el} not in ATOMS')
        if el == 'H':
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atoms.append((ATOMS[el], x, y, z))
    return mk_mols_atom_only(atoms)

def read_sdf_first_mol(fname):
    f = open(fname, 'rt')
    next(f)
    next(f)
    next(f)
    countline = next(f)
    it = countline.split()
    natoms = int(it[0])
    atoms = []
    for i in range(natoms):
        it = next(f).strip().split()
        x, y, z, el = float(it[0]), float(it[1]), float(it[2]), it[3].upper()
        if el not in ATOMS:
            raise Exception(f'{el} not in ATOMS')
        if el == 'H':
            continue
        atoms.append((ATOMS[el], x, y, z))
    return mk_mols_atom_only(atoms)

def print_mol(mol):
    for atom in mol.GetAtoms():
        print(atom.GetSymbol(), end=' ')
    print()

def get_rmsd(mol1, mol2):
    order1 = [x[1] for x in sorted(list(zip(Chem.CanonicalRankAtoms(mol1), range(mol1.GetNumAtoms()))))]
    order2 = [x[1] for x in sorted(list(zip(Chem.CanonicalRankAtoms(mol2), range(mol2.GetNumAtoms()))))]
    conf1 = mol1.GetConformer(0)
    conf2 = mol2.GetConformer(0)
    s = 0
    for i, j in zip(order1, order2):
        c1 = np.array(conf1.GetAtomPosition(i))
        c2 = np.array(conf2.GetAtomPosition(j))
        d2 = np.sum((c1 - c2)*(c1 - c2))
        s += d2
    rmsd = math.sqrt(s/mol1.GetNumAtoms())
    return rmsd

def connect_nearest(mol):
    rwmol = Chem.RWMol(mol)
    conf = rwmol.GetConformer(0)
    coords = np.array(([list(conf.GetAtomPosition(i)) for i in range(rwmol.GetNumAtoms())]))
    thres = 1.7
    thres2 = thres*thres
    connected = set()
    for i in range(coords.shape[0]):
        r = coords[i]
        for j in range(coords.shape[0]):
            if j == i:
                continue
            c = coords[j]
            v = np.sum((c - r)*(c - r))
            if v < thres2:
                if (i, j) in connected:
                    continue
                rwmol.AddBond(i, j, Chem.BondType.SINGLE)
                connected.add((i, j))
                connected.add((j, i))
    mol = Chem.Mol(rwmol)
    Chem.SanitizeMol(mol)
    return mol

def main(args):
    ref_name = args.ref_name
    fit_name = args.fit_name
    mol1 = read_pdb(ref_name)
    mol1 = connect_nearest(mol1)
    mol2 = read_sdf_first_mol(fit_name)
    mol2 = connect_nearest(mol2)
    rmsd = get_rmsd(mol1, mol2)
    print(rmsd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_name', type=str)
    parser.add_argument('fit_name', type=str)
    args = parser.parse_args()
    main(args)

