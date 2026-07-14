import unittest

import numpy as np

import EOSgenerators.RMF_DDH as DDH
from TOVsolver.unit import MeV, fm


class DDHChar23Tests(unittest.TestCase):
    def test_builtin_char23_matches_user_defined_form(self):
        oneoverfm_mev = 197.33
        m_sig = 550.0 / oneoverfm_mev
        m_w = 782.6 / oneoverfm_mev
        m_rho = 769.0 / oneoverfm_mev
        rho0_solver = 0.153

        (
            char_gs,
            char_gw,
            char_gr,
            char_dgs,
            char_dgw,
            char_dgr,
        ) = DDH.Function(type="Char23")

        (
            as_user,
            av_user,
            ar_user,
            bs_user,
            bv_user,
            br_user,
            cs_user,
            cv_user,
            cr_user,
            ds_user,
            dv_user,
            dr_user,
            rho0_user,
        ) = [
            8.225494,
            10.426752,
            0.64584657,
            2.7079569,
            1.6468675,
            5.2033131,
            2.4776689,
            6.8349408,
            0.4262597,
            3.8630221,
            1.4458185,
            -0.1824181,
            0.16194209,
        ]
        n0 = 0.16
        coupling_functions = [
            f"{as_user} + ({bs_user} + {ds_user}*(x/{n0})**(3))*exp(-{cs_user}*x/{n0})",
            f"{av_user} + ({bv_user} + {dv_user}*(x/{n0})**(3))*exp(-{cv_user}*x/{n0})",
            f"({ar_user} + ({br_user} + {dr_user}*(x/{n0})**(3))*exp(-{cr_user}*x/{n0}))*2",
        ]
        (
            user_gs,
            user_gw,
            user_gr,
            user_dgs,
            user_dgw,
            user_dgr,
        ) = DDH.Function(type="UserDefined", couplings=coupling_functions)

        rho_grid = np.array([0.08, 0.16, 0.32, 0.64])
        np.testing.assert_allclose(char_dgs(rho_grid), user_dgs(rho_grid))
        np.testing.assert_allclose(char_dgw(rho_grid), user_dgw(rho_grid))
        np.testing.assert_allclose(char_dgr(rho_grid), user_dgr(rho_grid))

        theta_char23 = np.array(
            [m_sig, m_w, m_rho, char_gs, char_gw, char_gr, char_dgs, char_dgw, char_dgr, rho0_solver],
            dtype=object,
        )
        theta_user = np.array(
            [m_sig, m_w, m_rho, user_gs, user_gw, user_gr, user_dgs, user_dgw, user_dgr, rho0_user],
            dtype=object,
        )

        eos_char23 = DDH.compute_eos([0], [0], theta_char23)
        eos_user = DDH.compute_eos([0], [0], theta_user)

        np.testing.assert_allclose(
            eos_char23[1] / (MeV * fm ** (-3)),
            eos_user[1] / (MeV * fm ** (-3)),
            rtol=1.0e-11,
            atol=1.0e-11,
        )
        np.testing.assert_allclose(
            eos_char23[2] / (MeV * fm ** (-3)),
            eos_user[2] / (MeV * fm ** (-3)),
            rtol=1.0e-11,
            atol=1.0e-11,
        )


if __name__ == "__main__":
    unittest.main()
