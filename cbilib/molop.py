from openeye.oechem import *
import numpy as np

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
        return False, None, None
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
    return found, ligand, protein

def make_pocket(ligand, apo):
    ligand_coords = np.array([v[1] for v in sorted(ligand.GetCoords().items())])
    protein_coords = np.array([v[1] for v in sorted(apo.GetCoords().items())])
    #print(ligand_coords.shape)
    #print(protein_coords.shape)
    D = np.zeros((len(protein_coords), len(ligand_coords)))
    for i in range(len(protein_coords)):
        for j in range(len(ligand_coords)):
            #print(protein_coords[i])
            #print(ligand_coords[j])
            v = protein_coords[i] - ligand_coords[j]
            D[i, j] = np.sum(v*v)
    D = D < 100
    v = []
    for i in range(len(protein_coords)):
        if np.any(D[i]):
            v.append(i)
    pocket_info = set()
    for i in v:
        res = OEAtomGetResidue(apo.GetAtom(OEHasAtomIdx(i)))
        chain_id = res.GetChainID()
        res_num = res.GetResidueNumber()
        pocket_info.add((chain_id, res_num))
    hv = OEHierView(apo)
    pocket = OEGraphMol()
    for chain in hv.GetChains():
        chain_id = chain.GetChainID()
        for frag in chain.GetFragments():
            for res in frag.GetResidues():
                resnum = res.GetResidueNumber()
                if (chain_id, resnum) in pocket_info:
                    for atom in res.GetAtoms():
                        pocket.NewAtom(atom)
    return pocket                

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
