import pandas as pd

def get_catalog(csv_name):
    df = pd.read_csv(csv_name)
    return df

def get_ligandname_from_catalog(df, pdbid):
    ref_ligandname = df[df['PDB_code'] == pdbid]['ligand_name'].item()[1:-1]
    return ref_ligandname
