from openeye.oechem import *

def read_pdb(fname):
    ifs = oemolistream(fname)
    pdbmol = OEGraphMol()
    OEReadMolecule(ifs, pdbmol)
    return pdbmol

def write_molecule(mol, oname):
    ofs = oemolostream(oname)
    OEWriteMolecule(ofs, mol)
    ofs.close()
