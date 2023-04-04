from typing import Union
from typing import List

import datamol as dm
import numpy as np
import padelpy

from rdkit.Chem import rdchem
from molfeat.utils.commons import requires_standardization
from molfeat.utils.datatype import to_numpy
from molfeat.calc.base import SerializableCalculator


class PadelDescriptors(SerializableCalculator):
    """
    Compute Padel descriptors for a molecule using padelpy (see https://github.com/ecrl/padelpy)
    For the original padel descriptor software, see: http://www.yapcwsoft.com/dd/padeldescriptor/
    Padel descriptors are under the GNU AGPL license v3.0.
    The descriptor calculator does not mask errors in featurization and will propagate them.

    !!! note:
        Padel Descriptors can be very slow to compute, as an exception, this calculator supports batch computation.
        It's recommended to process smiles in batch for large datasets and save them to disk.

    """

    def __init__(
        self,
        descriptors: bool = True,
        fingerprints: bool = True,
        timeout: int = 30,
        replace_nan: bool = False,
        do_not_standardize: bool = False,
    ):
        """
        Padel descriptor computation
        Args:
            descriptors: Whether to compute descriptors
            fingerprints: Whether to compute fingerprints
            timeout: Timeout for the computation
            do_not_standardize: Whether to force standardize molecules or keep it the same
        """
        self.descriptors = descriptors
        self.fingerprints = fingerprints
        self.timeout = timeout
        self.replace_nan = replace_nan
        self.do_not_standardize = do_not_standardize
        self._length = None
        self._columns = None
        self.__init_meta_data()

    def __getstate__(self):
        """Serialize the class for pickling."""
        state = {}
        state["replace_nan"] = self.replace_nan
        state["fingerprints"] = self.fingerprints
        state["timeout"] = self.timeout
        state["descriptors"] = self.descriptors
        state["do_not_standardize"] = getattr(self, "do_not_standardize", False)
        return state

    def __setstate__(self, state: dict):
        """Reload the class from pickling."""
        self.__dict__.update(state)
        self.__init_meta_data()

    def __init_meta_data(self):
        """
        Initialize the padel descriptors meta data for the current object
        """
        toy_smiles = "CCCC"
        out = padelpy.from_smiles(
            toy_smiles,
            fingerprints=self.fingerprints,
            descriptors=self.descriptors,
            timeout=self.timeout,
        )
        self._length = len(out)
        self._columns = [f"PaDEL_{name}" for name in out.keys()]

    @property
    def columns(self):
        """
        Get the name of all the descriptors of this calculator
        """
        return self._columns

    def __len__(self):
        """Return the length of the calculator"""
        return self._length

    @requires_standardization(disconnect_metals=True, remove_salt=True)
    def __call__(self, mol: Union[rdchem.Mol, str, List[rdchem.Mol], List[str]]):
        r"""
        Get rdkit basic descriptors for a molecule

        Args:
            mol: the molecule(s) of interest

        Returns:
            props (np.ndarray): list of computed mordred molecular descriptors
        """
        if isinstance(mol, (str, rdchem.Mol)):
            mol = [mol]
        smiles = [dm.to_smiles(dm.to_mol(m)) for m in mol]
        vals = padelpy.from_smiles(
            smiles,
            fingerprints=self.fingerprints,
            descriptors=self.descriptors,
            timeout=self.timeout,
        )

        def _as_float(x):
            try:
                x = float(x)
            except ValueError:
                x = float("nan")
            return x

        vals = to_numpy([[_as_float(x) for x in v.values()] for v in vals])
        if self.replace_nan:
            vals = np.nan_to_num(vals)
        if len(smiles) == 1:
            return vals[0]
        return vals
