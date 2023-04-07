from typing import Optional
from typing import Union
from typing import Dict
from typing import Any
from typing import Callable
from typing import List

import numpy as np
import datamol as dm
from collections.abc import Iterable
from molfeat.trans.base import MoleculeTransformer
from molfeat_padel.calc.padel import PadelDescriptors


class PadelTransformer(MoleculeTransformer):
    """
    A faster implementation of PadelTransformer that does not relies on
    individual calls to padelpy for each molecule
    """

    def __init__(
        self,
        n_jobs: int = 1,
        verbose: bool = False,
        dtype: Optional[Union[str, Callable]] = None,
        parallel_kwargs: Optional[Dict[str, Any]] = None,
        **params,
    ):
        """
        Batch padel descriptor computation

        Args:
            n_jobs: Number of jobs to use for the computation. -1 means all cores.
            dtype: The dtype of the output
            parallel_kwargs: Extra keyword arguments to pass to the parallel backend of datamol
            **params: Extra keyword arguments to pass to the padel descriptor calculator.
        """

        featurizer = PadelDescriptors(**params)
        super().__init__(
            featurizer, n_jobs=n_jobs, verbose=verbose, dtype=dtype, parallel_kwargs=parallel_kwargs
        )

    def transform(
        self,
        mols: List[Union[dm.Mol, str]],
        ignore_errors: bool = False,
        **kwargs,
    ):
        r"""
        Compute the features for a set of molecules.

        Args:
            mols: a list containing smiles or mol objects
            ignore_errors (bool, optional): Whether to silently ignore errors

        Returns:
            features: a list of features for each molecule in the input set
        """
        # Convert single mol to iterable format
        if isinstance(mols, (str, dm.Mol)) or not isinstance(mols, Iterable):
            mols = [mols]

        parallel_kwargs = getattr(self, "parallel_kwargs", {})
        features = self.featurizer.batch_compute(
            mols, n_jobs=self.n_jobs, parallel_kwargs=parallel_kwargs
        )
        if not ignore_errors:
            for ind, feat in enumerate(features):
                if feat is None or np.any(np.isnan(feat)):
                    raise ValueError(f"Cannot transform molecule at index {ind}")
        return features
