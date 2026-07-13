import numpy as np
from pathlib import Path
import EOSgenerators.RMF_DDH as DDH
import TOVsolver.main as tov
import EOSgenerators.crust_EOS as crust
from TOVsolver.unit import g_cm_3, dyn_cm_2, km, Msun, MeV, fm
from TOVsolver.constant import oneoverfm_MeV,G,c
from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline, interp1d
import pandas as pd


def _find_tolos_crust_file():
    """Return the bundled or local Tolos crust table path."""
    filename = "Tolos_crust_out.txt"
    module_dir = Path(__file__).resolve().parent
    repo_root = module_dir.parent
    candidates = [
        Path.cwd() / filename,
        Path.cwd() / "Test_Case" / filename,
        repo_root / "Test_Case" / filename,
        module_dir / "data" / filename,
    ]

    checked = []
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate in checked:
            continue
        checked.append(candidate)
        if candidate.exists():
            return candidate

    paths = ", ".join(str(path) for path in checked)
    raise FileNotFoundError(f"Could not find {filename}. Checked: {paths}")


class Likelihood():
   
    """
    A class designed for the DDH (Density-Dependent Coupling for Hadronic EOS) model, used to calculate 
    log-likelihoods for nuclear saturation properties and the maximum mass of neutron stars. 
    
    Nuclear saturation properties are calculated for quantities such as:
    - Nuclear binding energy (e0)
    - Nuclear matter incompressibility (k0)
    - Nuclear symmetry energy (jsym0) at saturation density
    
    Default values for these properties are set in `self.NMP_Base`, defined as:
        {"e0": [-16, 0.2, 'g'], "k0": [230, 40, 'g'], "jsym0": [32.5, 1.8, 'g'], 
         "de0": [0, 0.3, 'g'], "M": [2.0, 0.001, 's']}
    
    - The third entry in each property ("g" or "sg") defines the likelihood function type:
      - "g" for Gaussian likelihood, with the first value as mean (μ) and the second as standard deviation (σ).
      - "sg" for Super Gaussian, allowing definition of lower and upper limits.
    - The derivative of e0 (denoted by "de0") should be zero at saturation, as symmetric matter reaches a minimum here; 
      the default value is set to 0 with a tolerance of 0.3.

    Maximum neutron star mass (denoted by "M") is constrained with a default of 2.0 solar masses. For any equation of 
    state (EOS) yielding a maximum mass above 2.0, the log-likelihood cost is set to zero.

    ### Key Methods:
    - `update_base`: Updates the base values for nuclear saturation properties.
    - `update_quantities`: Toggles calculation of specific quantities, such as "lsym0", "ksym0", etc., which are 
      off by default but can be activated as needed.

    Parameters:
    - theta (array): Array of coupling parameters and other physical properties.
    - crust (array): Array defining the crust EOS, required for the full EOS calculation.

    Attributes:
    - `NMP_Base` (dict): Default values and likelihood types for nuclear saturation properties.
    - `NMP_Quantities_To_Compute` (dict): Dict to toggle on/off specific quantities for computation.
    - `Crust` (array): Crust EOS data.
    - `EoS` (array): Combined EOS after crust interpolation.
    - `MR` (array): Mass-radius relation output, used for neutron star mass calculations.

    ### Example Usage:
    ```
    ddh_likelihood = DDH_Likelihood(theta, crust)
    log_likelihood_value = ddh_likelihood.compute_NMP()
    max_mass_likelihood = ddh_likelihood.compute_MR()
    ```
    """
    
    # Define class attributes to store crust data
    eps_com = None
    pres_com = None

    def __init__(self, theta):
        self.theta    = theta
        self.NMP_Base = {"e0":[-16, 0.2, 'g'], "k0":[230, 40, 'g'], 
                         "jsym0":[32.5, 1.8, 'g'], "de0":[0, 0.3, 'g'], "M":[2.0, 0.001, "s"]}
        self.NMP_Quantities_To_Compute = {"de0" :True, "e0"   :True , "k0"   :True,
                                          "jsym0":True, "lsym0":False, "ksym0":False, "qsym0":False, "zsym0":False }

        self.EoS   = None
        self.MR    = None
        
        if Likelihood.eps_com is None or Likelihood.pres_com is None:
            Likelihood.initialize_crust_data()

    @classmethod
    def initialize_crust_data(cls):
        """Load and interpolate crust data only once."""
        Tolos_crust_out = np.loadtxt(_find_tolos_crust_file())
        eps_crust_T_out = Tolos_crust_out[:, 3] * g_cm_3
        pres_crust_T_out = Tolos_crust_out[:, 4] * dyn_cm_2
        cls.eps_com, cls.pres_com = crust.PolyInterpolate(eps_crust_T_out, pres_crust_T_out)


    def compute_eos(self):
        self.EoS = DDH.compute_eos(self.eps_com, self.pres_com, self.theta)
        self.EoS = [np.array([*self.eps_com, *self.EoS[1]]), np.array([*self.pres_com, *self.EoS[2]])]
      
    def update_base(self, quantity, new_base):
        self.NMP_Base[quantity] = new_base
        return
    
    def update_quantities(self, quantity, boolean):
        self.NMP_Quantities_To_Compute[quantity] = boolean
    
    def gaussian(self, quantity, x):
        mu, sigma, key = self.NMP_Base[quantity]
        # Calculate the normalization term: -0.5 * log(2 * pi * sigma^2)
        log_of_norm = -0.5 * np.log(2 * np.pi * sigma**2)
        
        # Calculate the exponent term: -0.5 * ((x - mu) / sigma)^2
        log_of_exp = -0.5 * ((x - mu) / sigma) ** 2
        
        # Return the sum of both terms, representing the log of the Gaussian PDF
        return log_of_norm + log_of_exp
    
    def super_gaussian(self, quantity, x):
        mu, sigma, key = self.NMP_Base[quantity]
        
        # Compute the first denominator term, incorporating the enlarged standard deviation
        denom_1 = 2 * sigma ** 2
        
        # Compute the second denominator term, which dampens large deviations from the mean
        denom_2 = 1 + np.exp((np.abs(x - mu) - sigma) / 0.015)
        
        # Calculate the Super-Gaussian PDF
        fx = 1 / (denom_1 * denom_2)
    
        log_fx = np.log(fx)
    
        return np.clip(log_fx, -1e21, np.inf)

    def smooth_step(self, quantity, x):
        mu, sigma, key = self.NMP_Base[quantity]

        operator = (mu - x)/sigma
        value    = (1 + np.exp(operator) )**(-1)
        
        return np.log(value)

    
    def compute_NMP(self):
        
        NMP = {}

        rho0=self.theta[-1]

        theta   = np.concatenate((self.theta, np.array([0.5])))

        EOS_SNM = DDH.compute_eos_alpha(theta)

        EA = (EOS_SNM[1]/EOS_SNM[0])*oneoverfm_MeV - 939.0

        eos1 = pd.DataFrame({
        'rho': EOS_SNM[0],
        'e'  : EOS_SNM[1]*oneoverfm_MeV,
        'p'  : EOS_SNM[2]*oneoverfm_MeV,
        'EA' : EA
        })

        fun1 = InterpolatedUnivariateSpline(eos1.rho, eos1.EA, k=3,ext=0)

        h_vec = 1e-4
        fhx   = fun1(eos1.rho)
        fhp   = fun1(eos1.rho+h_vec)
        fhm   = fun1(eos1.rho-h_vec)

        dedr         = (fhp-fhm)/(2.0*h_vec)
        eos1["dedr"] = dedr
        d2ed2r       = (fhp-2.0*fhx+fhm)/h_vec**2
        eos1["K"]    = 9*d2ed2r

        df_slice=eos1.query("0.05<rho<0.4")

        fun2        = InterpolatedUnivariateSpline(df_slice.rho, df_slice.dedr,k=3,ext=1)
        de0         = float(fun2(rho0))
        NMP["de0"]  = de0

        fun_e0    = InterpolatedUnivariateSpline(df_slice.rho,df_slice.EA, k=3,ext=0)
        e0        = float(fun_e0(rho0))
        NMP["e0"] = e0

        fun_k0    = InterpolatedUnivariateSpline(df_slice.rho,df_slice.K, k=3,ext=0)
        k0        = fun_k0(rho0)*rho0**2
        NMP["k0"] = k0



        ## Calculate PNM prop ........
        theta[-1] = 0.
        EOS_PNM   = DDH.compute_eos_alpha(theta)
        EA_pnm    = (EOS_PNM[1]/EOS_PNM[0])*oneoverfm_MeV - 939.0

        eos2 = pd.DataFrame({
        'rho': EOS_PNM[0],
        'e'  : EOS_PNM[1]*oneoverfm_MeV,
        'p'  : EOS_PNM[2]*oneoverfm_MeV,
        'EA' : EA_pnm
        })

        eos2["S"] = eos2.EA-eos1.EA
        fun1      = InterpolatedUnivariateSpline(eos2.rho, eos2.S, k=3,ext=0)

        h_vec = 1e-4
        fhx   = fun1(eos2.rho)
        fhp   = fun1(eos2.rho+h_vec)
        fhm   = fun1(eos2.rho-h_vec)

        dedr   = (fhp-fhm)/(2.0*h_vec)
        d2ed2r = (fhp-2.0*fhx+fhm)/h_vec**2

        eos2["L"]    = 3*rho0*dedr
        eos2["Ksym"] = 9*rho0**2*d2ed2r

        jsym0        = float(fun1(rho0))
        NMP["jsym0"] = jsym0


        fun_lsym0    = InterpolatedUnivariateSpline(eos2.rho,eos2.L, k=3,ext=0)
        lsym0        = float(fun_lsym0(rho0))
        NMP["lsym0"] = lsym0


        fun_ksym0    = InterpolatedUnivariateSpline(eos2.rho,eos2.Ksym, k=3,ext=0)
        ksym0        = float(fun_ksym0(rho0))
        NMP["ksym0"] = ksym0

        log_likelihood_value = 0
        for quantity in self.NMP_Quantities_To_Compute.keys():
            if self.NMP_Quantities_To_Compute[quantity]:
                if   self.NMP_Base[quantity][2] == 'g':
                    log_likelihood_value += self.gaussian(quantity , NMP[quantity])
                elif self.NMP_Base[quantity][2] == 'sg':
                    log_likelihood_value += self.super_gaussian(quantity , NMP[quantity])

        return log_likelihood_value 

    def compute_MR(self):
        if type(self.MR) == type(None):
            self.compute_eos()
            self.MR    = tov.OutputMR("", self.EoS[0], self.EoS[1]).T
            #self.M_max = np.max(self.MR[0])/Msun
            #self.R_max = self.MR[1][np.argwhere(self.MR[0]/Msun == self.M_max)][0][0] / km
            if self.MR[0].size > 0:
                self.M_max = np.max(self.MR[0]) / Msun
                self.R_max = self.MR[1][np.argwhere(self.MR[0] / Msun == self.M_max)][0][0] / km
            else:
                self.M_max = 0  
                self.R_max = 0     
        log_likelihood_value = self.smooth_step("M", self.M_max)

   
        return log_likelihood_value
