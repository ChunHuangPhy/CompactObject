import unittest

import numpy as np

import EOSgenerators.RMF_DDH as DDH
import EOSgenerators.RMF_EOS as RMF
from TOVsolver.maxium_central_density import maxium_central_density
from TOVsolver.maximum_central_density import maximum_central_density


class RMFValidationTests(unittest.TestCase):
    def test_ddh_theta_is_rejected_by_rmf_solver(self):
        gsf, gwf, grf, dgsf, dgwf, dgrf = DDH.Function(type="Typel99")
        theta_ddh = np.array(
            [2.77, 3.97, 3.87, gsf, gwf, grf, dgsf, dgwf, dgrf, 0.153],
            dtype=object,
        )

        with self.assertRaisesRegex(TypeError, "length-10 numeric RMF theta"):
            RMF.compute_EOS(np.array([1.0]), np.array([1.0]), theta_ddh)

    def test_misspelled_central_density_import_still_aliases_canonical_function(self):
        self.assertIs(maxium_central_density, maximum_central_density)


if __name__ == "__main__":
    unittest.main()
