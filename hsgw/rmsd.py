import sys, argparse, math
from rdkit import Chem
import numpy as np

def read_mol2(fname):
    mol = Chem.MolFromMol2File(fname, sanitize=False)
    Chem.SanitizeMol(mol, Chem.SANITIZE_CLEANUP | Chem.SANITIZE_CLEANUPCHIRALITY)
    mol = Chem.RemoveHs(mol)
    return mol

def read_pdb(fname):
    mol = Chem.MolFromPDBFile(fname, sanitize=False)
    Chem.SanitizeMol(mol, Chem.SANITIZE_CLEANUP | Chem.SANITIZE_CLEANUPCHIRALITY)
    mol = Chem.RemoveHs(mol)
    return mol

def read_sdf(fname):
    for mol in Chem.SDMolSupplier(fname, sanitize=False):
        Chem.SanitizeMol(mol, Chem.SANITIZE_CLEANUP | Chem.SANITIZE_CLEANUPCHIRALITY)
        mol = Chem.RemoveHs(mol)
        yield mol

def print_mol(mol):
    for atom in mol.GetAtoms():
        print(atom.GetSymbol(), end=' ')
    print()

def get_rmsd(mol1, mol2):
    smi1 = Chem.MolToSmiles(mol1)
    smi2 = Chem.MolToSmiles(mol2)
    if smi1 != smi2:
        print(smi1, 'ref', file=sys.stderr)
        print(smi2, 'docked did not match', file=sys.stderr)
        return -1.0
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

def main(args):
    ref_name = args.ref_name
    fit_name = args.fit_name
    #for mol1 in read_sdf(ref_name):
    #    break
    mol1 = read_pdb(ref_name)
    count = 0
    for mol2 in read_sdf(fit_name):
        count += 1
        rmsd = get_rmsd(mol1, mol2)
        rmsd = round(rmsd, 3)
        print(count, rmsd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_name', type=str)
    parser.add_argument('fit_name', type=str)
    args = parser.parse_args()
    main(args)

