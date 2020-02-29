import prody
import os
import numpy as np

def fetch_pdb(pdbid, dir):
    cwd = os.getcwd()
    os.chdir(dir)
    prody.fetchPDB(pdbid)
    os.chdir(cwd)

def get_distance_matrix(ag1, ag2):
    a = ag1.getCoords()
    b = ag2.getCoords()
    x = np.tile(a, (1, b.shape[0])).reshape(a.shape[0], b.shape[0], 3)
    y = np.tile(b, (a.shape[0], 1)).reshape(a.shape[0], b.shape[0], 3)
    d = x - y
    D = np.sqrt(np.sum(d*d, axis=2))
    return D

def get_chains(protein, ligand, thres=5.0):
    dmat = get_distance_matrix(protein, ligand)
    cmat = dmat < thres
    atomidxs = []
    for i in range(cmat.shape[0]):
        if np.any(cmat[i]):
            atomidxs.append(i)
    chains = set()
    for i in atomidxs:
        atom = protein[i]
        chainid = atom.getChid()
        chains.add(chainid)
    return sorted(chains)

def ligand_pick_one(ligand, thres=2.2):
    dmat = get_distance_matrix(ligand, ligand)
    conn = []
    for i, j in zip(*np.where(dmat < thres)):
        if i != j:
            conn.append((i, j))
    conn.sort()
    grps = []
    while conn:
        m, n = conn.pop(0)
        if len(grps) == 0:
            grps.append(set())
            g = grps[-1]
            g.add(m)
            g.add(n)
        else:
            found = False
            for g in grps:
                if m in g or n in g:
                    g.add(m)
                    g.add(n)
                    found = True
            if not found:
                grps.append(set())
                g = grps[-1]
                g.add(m)
                g.add(n)
    grps.sort(key=lambda x: (-len(x), x))
    for i in range(len(grps)):
        grps[i] = sorted(grps[i])
    if len(grps) == 0:
        return None
    head = grps[0]
    if len(head) == 0:
        return None
    ligand = ligand[head]
    return ligand

def apo_near_ligand(apo, ligand, thres=5.0):
    chains = get_chains(apo, ligand, thres=thres)
    atomidxs = []
    for i in range(apo.numAtoms()):
        atom = apo[i]
        chid = atom.getChid()
        if chid in chains:
            atomidxs.append(i)
    if len(atomidxs) == 0:
        return None
    apo = apo[atomidxs]
    return apo

def make_pocket(apo, ligand, thres=5.0):
    dmat = get_distance_matrix(apo, ligand)
    cmat = dmat < thres
    
    atomidxs = []
    for i in range(cmat.shape[0]):
        if np.any(cmat[i]):
            atomidxs.append(i)
    ress = set()
    for i in atomidxs:
        atom = apo[i]
        resnum = atom.getResnum()
        resname = atom.getResname()
        chainid = atom.getChid()
        ress.add((resnum, resname, chainid))
    ress = sorted(ress)

    atomidxs = []
    for i in range(apo.numAtoms()):
        atom = apo[i]
        resnum = atom.getResnum()
        resname = atom.getResname()
        chainid = atom.getChid()
        if (resnum, resname, chainid) in ress:
            atomidxs.append(i)
    if len(atomidxs) == 0:
        return None
    pocket = apo[atomidxs]
    return pocket
