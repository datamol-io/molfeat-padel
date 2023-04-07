import unittest as ut
import datamol as dm
import numpy as np

from molfeat.plugins import load_registered_plugins

load_registered_plugins()
from molfeat.trans import PadelTransformer

class TestPadel(ut.TestCase):
    r"""Test cases for descriptors and pharmacophore generation"""

    def test_padel_fp(self):
        smiles = dm.freesolv()["smiles"].values[:50]
        mol_transf = PadelTransformer(dtype="float")
        fps = mol_transf(smiles)
        self.assertEqual(fps.shape, (50, 2756))

        from molfeat_padel.calc.padel import PadelDescriptors as PD
        calc = PD()
        fps2 = calc.batch_compute(smiles)
        np.testing.assert_allclose(fps.values, fps2)


if __name__ == "__main__":
    ut.main()
