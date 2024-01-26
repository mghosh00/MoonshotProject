import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


class Molecule:
    """A class to represent Molecules and their Morgan fingerprints
    """

    def __init__(self, smiles: str):
        self.smiles = smiles
        self.molecule = self._create_molecule()
        self.bit_info = {}
        self.morgan_fp = None
        self.morgan = self.create_morgan()

    def _create_molecule(self) -> rdkit.Chem.Mol:
        """Produces a molecule from a SMILES string and returns an error if
        the string is invalid

        :return: The molecule
        """
        molecule = Chem.MolFromSmiles(self.smiles, sanitize=True)
        if molecule:
            return molecule
        else:
            raise ValueError("SMILES string is invalid")

    def create_morgan(self, radius: int = 2) -> np.array:
        """Creates the morgan fingerprint of the molecule and returns it
        in np array form

        :param radius: The radius for the fingerprint (optional)
        :return: The molecule as an np.array
        """
        morgan_fp = AllChem.GetMorganFingerprintAsBitVect(self.molecule,
                                                          useChirality=True,
                                                          radius=radius,
                                                          nBits = 2048,
                                                          bitInfo=self.bit_info)
        morgan_np = np.array(morgan_fp)
        self.morgan = morgan_np
        self.morgan_fp = morgan_fp
        return morgan_np

    def get_locs_of_fingerprints(self):
        """Gets the locations of the Morgan fingerprints
        """
        if self.morgan is not None:
            return np.nonzero(self.morgan)


molecule = Molecule('C(C[C@@H](C(=O)O)N)CNC(=N)N')
molecule.create_morgan()
