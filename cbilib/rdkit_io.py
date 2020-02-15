from rdkit import Chem

def read_sdf(fname, sanitize=False, removeHs=True):
    mol = Chem.MolFromMolFile(fname, sanitize=sanitize, removeHs=removeHs)
    return mol

def read_pdb(fname, sanitize=False, removeHs=True):
    mol = Chem.MolFromPDBFile(fname, sanitize=sanitize, removeHs=removeHs)
    return mol
