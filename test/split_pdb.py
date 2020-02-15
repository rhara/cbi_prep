from openeye.oechem import *
import os, argparse, pickle
import numpy as np

def read_pdb(fname):
    ifs = oemolistream(fname)
    pdbmol = OEGraphMol()
    OEReadMolecule(ifs, pdbmol)
    return pdbmol

def write_molecule(mol, oname):
    ofs = oemolostream(oname)
    OEWriteMolecule(ofs, mol)
    ofs.close()

def split_ligands(ligand):
    n, parts = OEDetermineComponents(ligand)
    pred = OEPartPredAtom(parts)
    submols = []
    for i in range(1, max(parts)+1):
        pred.SelectPart(i)
        submol = OEGraphMol()
        OESubsetMol(submol, ligand, pred, True)
        submols.append(submol)
    return submols

def compact_protein(protein, ligand):
    chains = set()
    hv = OEHierView(protein)
    for chain in hv.GetChains():
        for frag in chain.GetFragments():
            for res in frag.GetResidues():
                for atom in res.GetAtoms():
                    for lig_atom in ligand.GetAtoms():
                        if OEGetDistance(protein, atom, ligand, lig_atom) < 6.0:
                            chains.add(chain.GetChainID())
    compact_protein = OEGraphMol()
    for chain in hv.GetChains():
        if chain.GetChainID() not in chains:
            continue
        for frag in chain.GetFragments():
            for res in frag.GetResidues():
                for atom in res.GetAtoms():
                    compact_protein.NewAtom(atom)
    return compact_protein

def get_ligand_and_protein(pdbmol, ligandname):
    AMINO_ACIDS = 'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(' ')
    ligand = OEGraphMol()
    hv = OEHierView(pdbmol)
    found = False
    for chain in hv.GetChains():
        if found:
            break
        for frag in chain.GetFragments():
            for res in frag.GetResidues():
                if res.GetResidueName() == ligandname:
                    found = True
                    for atom in res.GetAtoms():
                        ligand.NewAtom(atom)
    if not found:
        raise Exception('Ligand {ligandname} was not found')
    ligand.SetDimension(3)
    OEDetermineConnectivity(ligand)
    OEFindRingAtomsAndBonds(ligand)
    OEPerceiveBondOrders(ligand)
    OEAssignImplicitHydrogens(ligand)
    OEAssignFormalCharges(ligand)
    protein = OEGraphMol()
    for chain in hv.GetChains():
        for frag in chain.GetFragments():
            for res in frag.GetResidues():
                if res.GetResidueName() in AMINO_ACIDS:
                    for atom in res.GetAtoms():
                        protein.NewAtom(atom)
    smi = OEMolToSmiles(ligand)
    if '.' in smi:
        ligands = split_ligands(ligand)
        smis = set()
        for i in range(len(ligands)):
            submol = ligands[i]
            smis.add(OECreateSmiString(submol))
        if 1 < len(smis):
            raise Exception('Ligand was not resolved')
        ligand = ligands[0]

    protein = compact_protein(protein, ligand)
    return found, ligand, protein

def get_center(mol):
    xyz = np.array(list(mol.GetCoords().values()))
    center = [round(x, 3) for x in list(np.mean(xyz, axis=0))]
    return center

def main(args):
    pdb_name = args.pdb
    pdbid = os.path.basename(pdb_name)[:4]
    ligandname = args.ligand
    pdbmol = read_pdb(pdb_name)
    found, ligand, protein = get_ligand_and_protein(pdbmol, ligandname)
    center = get_center(ligand)
    print(OEMolToSmiles(ligand), center)
    D = {}
    D['center'] = center
    write_molecule(ligand, f'{pdbid}_ligand.mol2')
    write_molecule(protein, f'{pdbid}_apo.pdb')
    pickle.dump(D, open(f'{pdbid}_info.pkl', 'wb'), protocol=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str)
    parser.add_argument('ligand', type=str)
    args = parser.parse_args()
    main(args)
