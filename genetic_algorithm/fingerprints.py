import numpy as np
from rdkit.Chem import DataStructs
from rdkit.Chem.AtomPairs import Pairs


class get_ap_fingerprintfn:
    def __call__(self, mols):
        return [Pairs.GetAtomPairFingerprint(x) for x in mols]

    def __repr__(self) -> str:
        return f"AtomPairFingerprint: get_ap_fingerprintfn()"


class get_ap_distancefn:
    def __call__(self, a, b):
        return 1 - DataStructs.DiceSimilarity(a,b)

    def __repr__(self) -> str:
        return f"AtomPairFingerprint Dice Distance: get_ap_distancefn()"
