from openeye.oechem import *
import sys
import numpy as np

ifs = oemolistream(sys.argv[1])
mol = OEGraphMol()
OEReadMolecule(ifs, mol)
ifs.close()

coords = np.array(list(mol.GetCoords().values()))
print(list(np.mean(coords, axis=0)))
