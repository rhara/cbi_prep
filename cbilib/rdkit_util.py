import pickle
from rdkit import Chem

def make_pair(lig_iname, pocket_iname, oname):
    ligand = Chem.MolFromMolFile(lig_iname, sanitize=False)
    pocket = Chem.MolFromPDBFile(pocket_iname, sanitize=False)
    pickle.dump((ligand, pocket), open(oname, 'wb'), protocol=4)
