import pandas as pd
from rdkit.Chem import DataStructs

from moonshot_project import Molecule

class DistanceCalculator:
    """Calculates Tanimoto distance matrix given a dict of SMILES strings
    keyed by compound_id
    """
    def __init__(self, smiles_dict: dict[str, str]):
        """

        :param smiles_dict: SMILES strings keyed by compound_id
        """
        self.smiles_dict = smiles_dict
        self.molecule_dict = self._create_molecule_dict()
        self.tanimoto_matrix = None

    def _create_molecule_dict(self):
        return {_id: Molecule(self.smiles_dict[_id])
                for _id in self.smiles_dict}

    @staticmethod
    def tanimoto_distance(mol1: Molecule, mol2: Molecule):
        return 1 - DataStructs.TanimotoSimilarity(mol1.morgan_fp,
                                              mol2.morgan_fp)

    def get_tanimoto_matrix(self):
        df = pd.DataFrame(columns=list(self.smiles_dict.keys()),
                          index=list(self.smiles_dict.keys()))
        for id_1 in self.molecule_dict.keys():
            for id_2 in self.molecule_dict.keys():
                mol_1 = self.molecule_dict[id_1]
                mol_2 = self.molecule_dict[id_2]
                df.at[id_1, id_2] = self.tanimoto_distance(mol_1, mol_2)
        self.tanimoto_matrix = df
        return df

    def write_to_csv(self):
        matrix = self.get_tanimoto_matrix()
        matrix.to_csv("tanimoto.csv")
