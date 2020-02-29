import pandas as pd

class Catalog:
    def __init__(self, csv_name):
        self.csv_name = csv_name
        self.df = pd.read_csv(csv_name)

    def get_ligname(self, pdbid):
        ligname = self.df[self.df['pdbid'] == pdbid]['ligname'].item()
        return ligname
