import unittest

import numpy as np

import TOVsolver.main as main
from TOVsolver.unit import Msun, dyn_cm_2, g_cm_3, km


class OutputMRFallbackTests(unittest.TestCase):
    def test_outputmr_uses_tidal_fallback_for_nonphysical_primary_results(self):
        density = np.linspace(1.0e14, 2.0e15, 30) * g_cm_3
        pressure = np.linspace(1.0e32, 1.0e35, 30) * dyn_cm_2

        original_solve_tov = main.TOV_solver.solveTOV
        original_solve_tov_tidal = main.TOV_solver.solveTOV_tidal

        def nonphysical_solve_tov(*args, **kwargs):
            return 1.0e-15 * Msun, 1.0e-4 * km

        def physical_solve_tov_tidal(center_density, energy_density, pressure):
            mass = (1.0 + 1.0e-16 * float(center_density / g_cm_3)) * Msun
            radius = 12.0 * km
            return mass, radius, 0.0

        try:
            main.TOV_solver.solveTOV = nonphysical_solve_tov
            main.TOV_solver.solveTOV_tidal = physical_solve_tov_tidal

            mr = main.OutputMR("", density, pressure)
        finally:
            main.TOV_solver.solveTOV = original_solve_tov
            main.TOV_solver.solveTOV_tidal = original_solve_tov_tidal

        self.assertEqual(mr.shape, (50, 2))
        self.assertTrue(np.all((mr[:, 0] / Msun) > 1.0))
        self.assertTrue(np.allclose(mr[:, 1] / km, 12.0))


if __name__ == "__main__":
    unittest.main()
