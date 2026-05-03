from NumbaMinpack import lmdif, hybrd, minpack_sig
from TOVsolver.unit import g_cm_3, dyn_cm_2
from numba import njit, jit, cfunc
from scipy import optimize
import numpy as np
import math

c = 2.99792458e10
G = 6.67428e-8
Msun = 1.989e33

dyncm2_to_MeVfm3 = 1.0 / (1.6022e33)
gcm3_to_MeVfm3 = 1.0 / (1.7827e12)
oneoverfm_MeV = 197.327053

m_e = 2.5896041 * 10**-3
m_mu = 0.5354479981
m_n = 4.7583690772
m_p = 4.7583690772

J_B = 1 / 2.0
b_B = 1

m_l = [m_e, m_mu]
m_b = [m_p, m_n]

Matrix_b = np.array(
    [[1.0, 1.0, 1 / 2.0, 1.0, 1.0, 1.0], [1.0, 0.0, -1 / 2.0, 1.0, 1.0, 1.0]]
)

Matrix_l = np.array([[0.0, -1.0, 1 / 2.0], [0.0, -1.0, 1 / 2.0]])


@jit
def initial_values(rho, theta):
    """Outputs the the sigma, omega, rho term and chemical potential of electron and neutron at
    given initial density.
    
    (Compiled Version)
    
    Args:
        rho (float): given nuclear density
        theta (array): parameters of determine a RMF model in Lagrangian; here, we have 10 parameters.

    Returns:
        sigma (float): sigma term in Lagrangian
        omega (float): omega term in Lagrangian
        rho_03 (float): rho term in Lagrangian
        mu_n (float): chemical potential of neutron matter
        mu_e (float): chemical potential of electron portion

    """
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w = theta

    m_e = 2.5896041 * 10**-3
    m_mu = 0.5354479981
    m_n = 4.7583690772
    m_p = 4.7583690772

    rho_0 = 0.1505

    sigma = g_sigma * rho / (m_sig**2)
    rho_03 = -g_rho * rho / (2.0 * (m_rho**2))
    omega = rho * (
        (((m_w**2) / g_omega) + (2.0 * Lambda_w * ((g_rho * rho_03) ** 2) * g_omega))
        ** (-1.0)
    )
    m_eff = m_n - (g_sigma * sigma)
    mu_n = m_eff + g_omega * omega + g_rho * rho_03 * Matrix_b[1, 2]
    mu_e = 0.12 * m_e * (rho / rho_0) ** (2 / 3.0)

    return sigma, omega, rho_03, mu_n, mu_e


@cfunc(minpack_sig)
def myfunc(x, fvec, args):
    """Outputs the the sigma, omega, rho term and chemical potential of electron and neutron at
    a given initial density,
    
    (Faster Version Using Numba)
    
    Args:
        rho (float): given nuclear density.
        theta (array): parameters of specific a RMF model in Lagrangian; here, we have 10 parameters.

    Returns:
        sigma (float): sigma term in Lagrangian.
        omega (float): omega term in Lagrangian.
        rho_03 (float): rho term in Lagrangian.
        mu_n (float): chemical potential of neutron matter.
        mu_e (float): chemical potential of electron portion.

    """
    m_sig = args[0]
    m_w = args[1]
    m_rho = args[2]
    g_sigma = args[3]
    g_omega = args[4]
    g_rho = args[5]
    kappa = args[6]
    lambda_0 = args[7]
    zeta = args[8]
    Lambda_w = args[9]
    rho = args[10]

    m_e = 2.5896041 * 10**-3
    m_mu = 0.5354479981
    m_n = 4.7583690772
    m_p = 4.7583690772

    J_B = 1 / 2.0
    b_B = 1

    m_l = np.array([m_e, m_mu])
    m_b = np.array([m_p, m_n])

    Matrix_b = np.array(
        [[1.0, 1.0, 1 / 2.0, 1.0, 1.0, 1.0], [1.0, 0.0, -1 / 2.0, 1.0, 1.0, 1.0]]
    )

    Matrix_l = np.array([[0.0, -1.0, 1 / 2.0], [0.0, -1.0, 1 / 2.0]])

    sigma = x[0]
    omega = x[1]
    rho_03 = x[2]
    mu_n = x[3]
    mu_e = x[4]

    rho_B_list = []
    rho_SB_list = []
    q_list = []

    m_eff = m_n - (g_sigma * sigma)
    # i = 0
    # j = 0
    for i in range(len(Matrix_b)):
        # while i < len(Matrix_b):

        mu_b = Matrix_b[i, 0] * mu_n - Matrix_b[i, 1] * mu_e

        E_fb = mu_b - g_omega * omega - g_rho * rho_03 * Matrix_b[i, 2]

        k_fb_sq = E_fb**2 - m_eff**2
        if k_fb_sq < 0:
            k_fb_sq = 0.0
            E_fb = m_eff

        k_fb = math.sqrt(k_fb_sq)

        rho_B = ((2 * J_B) + 1) * b_B * k_fb**3 / (6.0 * math.pi**2)
        rho_SB = (m_eff / (2.0 * math.pi**2)) * (
            E_fb * k_fb - (m_eff ** (2)) * np.log((E_fb + k_fb) / m_eff)
        )

        rho_B_list.append(rho_B)
        rho_SB_list.append(rho_SB)

        Q_B = ((2.0 * J_B) + 1.0) * Matrix_b[i, 1] * k_fb**3 / (6.0 * math.pi**2)
        q_list.append(Q_B)
        # i += 1

    for j in range(len(Matrix_l)):
        # while j < len(Matrix_l):

        mu_l = Matrix_l[j, 0] * mu_n - Matrix_l[i, 1] * mu_e
        E_fl = mu_l

        k_fl_sq = E_fl**2 - m_l[j]**2
        
        if k_fl_sq < 0.0:
            k_fl_sq = 0.0
    
        k_fl = math.sqrt(k_fl_sq)

        Q_L = ((2.0 * J_B) + 1.0) * Matrix_l[i, 1] * (k_fl**3) / (6.0 * (math.pi**2))
        q_list.append(Q_L)

        # rho_l = k_fl**3 / (3.*(math.pi**2))
        # rho_list.append(rho_l)
        
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    sum4 = 0.0
    sum5 = 0.0
    for i in range(len(rho_SB_list)):
        sum1 = sum1 + rho_SB_list[i] * Matrix_b[:,3][i]
    for i in range(len(rho_B_list)):
        sum2 = sum2 + rho_B_list[i] * Matrix_b[:,4][i]
    for i in range(len(rho_B_list)):
        sum3 = sum3 + rho_B_list[i] * Matrix_b[:,5][i] * Matrix_b[:,2][i]
    for i in range(len(rho_B_list)):
        sum4 = sum4 + rho_B_list[i]
    for i in range(len(q_list)):
        sum5 = sum5 + q_list[i]

    fvec[0] = (sigma*(m_sig**2)/g_sigma - sum1 +(kappa*(g_sigma*sigma)**2)/2.\
              + (lambda_0*(g_sigma*sigma)**3)/6.)**2
    fvec[1] = (omega*(m_w**2)/g_omega - sum2 + (zeta*(g_omega*omega)**3)/6.\
              + 2.*Lambda_w*g_omega*omega*(rho_03*g_rho)**2)**2
    fvec[2] = (rho_03*(m_rho**2)/g_rho - sum3 + 2.*Lambda_w*g_rho*rho_03*(omega*g_omega)**2)**2
    
    fvec[3] = (rho - sum4)**2
    fvec[4] = sum5**2

    
funcptr = myfunc.address # address in memory to myfunc

@njit
def Energy_density_Pressure(x, rho, theta, return_tag=False):
    """
    Computes the pressure and energy density for the equation of state (EOS) 
    based on the Relativistic Mean Field (RMF) model parameters,
    
    (Faster Version Using Numba)
    
    Args:
        x (array): An array containing the initial values for sigma, omega, rho, 
                   and chemical potential, obtained from the `initial_values` function.
        rho (float): The central density at which the EOS computation begins.
        theta (array): An array of 10 parameters that define the RMF model in the 
                       Lagrangian.
        return_tag (bool, optional): If False (default), returns only the energy 
                                     density and pressure. If True, returns additional 
                                     EOS components.

    Returns:
        tuple:
            If `return_tag` is False:
                energy_density (float): The energy density in natural units 
                                        (to convert to MeV/fm³, divide by MeV/fm³).
                pressure (float): The pressure in natural units.
            
            If `return_tag` is True:
                numpy array: A 1D array representing EOS components:
                    - EoS[0]: Number density in fm⁻³
                    - EoS[1]: Energy density in natural units
                    - EoS[2]: Pressure in natural units
                    - EoS[3]: Proton chemical potential in natural units
                    - EoS[4]: Neutron chemical potential in natural units
                    - EoS[5]: Electron chemical potential in natural units
                    - EoS[6]: Muon chemical potential in natural units
                    - EoS[7]: Proton fraction (dimensionless)
    """
    sigma, omega, rho_03, mu_n, mu_e = x

    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w = theta

    m_e = 2.5896041 * 10**-3
    m_mu = 0.5354479981
    m_n = 4.7583690772
    m_p = 4.7583690772

    J_B = 1 / 2.0
    b_B = 1

    m_l = np.array([m_e, m_mu])
    m_b = np.array([m_p, m_n])

    energy_b = 0
    energy_l = 0
    multi = 0

    EoS = np.empty(8, dtype=np.float64) #["rho", "energy density", "pressure", "mu_p", "mu_n", "mu_e", "mu_mu", "proton_fraction"]
    
    m_eff = m_n - (g_sigma * sigma)

    for i in range(len(Matrix_b)):
        mu_b = Matrix_b[i, 0] * mu_n - Matrix_b[i, 1] * mu_e

        EoS[3+i] = mu_b

        E_fb = mu_b - g_omega * omega - g_rho * rho_03 * Matrix_b[i, 2]

        k_fb_sq = E_fb**2 - m_eff**2
        if k_fb_sq < 0:
            k_fb_sq = 0.0
            E_fb = m_eff

        k_fb = math.sqrt(k_fb_sq)

        rho_B = ((2.0 * J_B) + 1.0) * b_B * (k_fb**3) / (6.0 * math.pi**2)

        if i == 0:
            #when b = proton,
            EoS[7] = rho_B / rho

        multi = multi + mu_b * rho_B
        energy_baryon = (1 / (8.0 * (math.pi**2))) * (
            k_fb * (E_fb**3)
            + (k_fb**3) * E_fb
            - np.log((k_fb + E_fb) / m_eff) * m_eff**4
        )

        energy_b = energy_b + energy_baryon

    for j in range(len(Matrix_l)):
        mu_l = Matrix_l[i, 0] * mu_n - Matrix_l[j, 1] * mu_e

        EoS[5+j] = mu_l

        k_fl_sq = mu_l**2 - m_l[j] ** 2
        if k_fl_sq < 0.0:
            k_fl_sq = 0.0
        k_fl = math.sqrt(k_fl_sq)

        rho_l = k_fl**3 / (3.0 * math.pi**2)

        multi = multi + mu_l * rho_l
        energy_lepton = (1 / (8.0 * (math.pi**2))) * (
            k_fl * (mu_l**3)
            + mu_l * (k_fl**3)
            - (m_l[j] ** 4) * np.log((k_fl + mu_l) / m_l[j])
        )

        energy_l = energy_l + energy_lepton

    sigma_terms = (
        0.5 * ((sigma * m_sig) ** 2)
        + (kappa * ((g_sigma * sigma) ** 3)) / 6.0
        + (lambda_0 * ((g_sigma * sigma) ** 4)) / 24.0
    )

    omega_terms = 0.5 * ((omega * m_w) ** 2) + (zeta * ((g_omega * omega) ** 4)) / 8.0

    rho_terms = 0.5 * ((rho_03 * m_rho) ** 2) + +3.0 * Lambda_w * (
        (g_rho * rho_03 * g_omega * omega) ** 2
    )

    energy_density = energy_b + energy_l + sigma_terms + omega_terms + rho_terms

    Pressure = multi - energy_density

    if return_tag:
        EoS[0] = rho
        EoS[1] = energy_density
        EoS[2] = Pressure
        return EoS
    else:
        EoS[1] = energy_density
        EoS[2] = Pressure
        return EoS

# define a function that computes the EoS 
def compute_EOS(eps_crust, pres_crust, theta, return_tag=False):
    """Generate core part equation of state, main function, from RMF model,

    (Faster Version Using Numba)
    
    Args:
        eps_crust (array): the energy density of crust EoS in g/cm³.
        pres_crust (array): the pressure from crust EoS model in dyn/cm².
        theta (array): An array representing the parameters used to determine a RMF model in the
        Lagrangian. In this case, the RMF model is defined by 10 parameters.

        return_tag (bool, optional): If False (default), returns only the energy 
                                     density and pressure. If True, returns additional 
                                     EOS components.
    Returns:
        If `return_tag` is False:
                energy_density (float): The energy density in natural units 
                                        (to convert to MeV/fm³, divide by MeV/fm³).
                pressure (float): The pressure in natural units.
            
        If `return_tag` is True:
                numpy array: A 1D array representing EOS components:
                    - EoS[0]: Number density in fm⁻³.
                    - EoS[1]: Energy density in natural units.
                    - EoS[2]: Pressure in natural units.
                    - EoS[3]: Proton chemical potential in natural units.
                    - EoS[4]: Neutron chemical potential in natural units.
                    - EoS[5]: Electron chemical potential in natural units.
                    - EoS[6]: Muon chemical potential in natural units.
                    - EoS[7]: Proton fraction (dimensionless).
    """
    dt    = 0.05
    rho_0 = 0.1505
    
    x_init   = np.array(initial_values(0.1 * rho_0, theta), dtype=np.float64)

    if return_tag:
        EoS = np.empty((124, 8), dtype=np.float64)
    else:
        Energy   = np.empty(124, dtype=np.float64)
        Pressure = np.empty(124, dtype=np.float64)
    for i in range(1, 125):

        rho = i * dt * rho_0
        args = np.append(theta, rho)
        xsol, fvec, success, info = lmdif(funcptr, x_init, 5, args)
        
        Re  = Energy_density_Pressure(x_init, rho, theta, return_tag)

        if return_tag:
            # Re = [ rho , energy_density , pressure, mu_n , mu_p , mu_e , mu_mu , proton_fraction ]
            Re[1] = Re[1] * oneoverfm_MeV / gcm3_to_MeVfm3
            Re[2] = Re[2] * oneoverfm_MeV / dyncm2_to_MeVfm3
            Re[3] = Re[3]
            Re[4] = Re[4]
            Re[5] = Re[5]
            Re[6] = Re[6]
            
            EoS[i-1] = Re
        else:
            Energy[i-1]   = Re[1]*oneoverfm_MeV / gcm3_to_MeVfm3
            Pressure[i-1] = Re[2]*oneoverfm_MeV / dyncm2_to_MeVfm3

        x_init = xsol
    
    if return_tag:
        end = 0
        for i in range(0, len(EoS) - 1):
            if EoS[i][1] > max(eps_crust / g_cm_3) and i > 18:
                end = i + 2
                break
            end += 1
        EoS = EoS[end::].T
        EoS[1] = EoS[1] * g_cm_3
        EoS[2] = EoS[2] * dyn_cm_2

        return EoS
    else:
        end = 0
        for i in range(0, len(Energy) - 1):
            if Energy[i] > max(eps_crust / g_cm_3) and i > 18:
                end = i + 2
                break
            end += 1
        ep = Energy[end::]
        pr = Pressure[end::]
    
        # tzzhou: migrating to new unit convention
        ep = ep * g_cm_3
        pr = pr * dyn_cm_2
    
        return ep, pr


##################### Date: 04 Nov 2024 #######################
###  João Cartaxo ### Tuhin Malik ### Constança Providência ###
def initial_guess_alpha(rho, theta):
    """ Outputs the sigma, omega, rho field value,
    
    (Faster Version Using Numba)
    
    Args:
        rho (float): given nuclear density
        theta (array): parameters to determine an RMF model in Lagrangian, here there are 11 parameters,
        where the last parameters are the proton fraction (alpha) and the number density rho.

    Returns:
        math.sqrt(sigma) (float): square root of the sigma term in the Lagrangian.
        math.sqrt(omega) (float): square root of the omega term in the Lagrangian.
        rho_03 (float): rho term in the Lagrangian.
    """
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w, alpha = theta
        
    sigma = g_sigma*rho/(m_sig**2)
    rho_03 = -g_rho*rho/(2.*(m_rho**2))
    omega = rho*((((m_w**2)/g_omega)+\
    (2.*Lambda_w*((g_rho*rho_03)**2)*g_omega))**(-1.))

    return math.sqrt(sigma), math.sqrt(omega), rho_03
    
@njit
def fields_alpha(x, args):
    """ Iterate the sigma, omega, and rho fields for a given proton fraction and density.
    
    (Faster Version Using Numba)
    
    Args:
        x (array): initial sqrt(sigma) sqrt(omega) and rho from initial_values function.
        args (array): parameters to determine a RMF model in Lagrangian; here, we have 12 parameters,
        where the last parameters are the proton fraction (alpha) and the density rho. 
        For pure neutron matter (PNM), alpha is 0, and for symmetric nuclear matter, alpha is 0.5.

    Returns:
        f (array): field equations which are then solved using the scipy root finding function.
    """
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w, alpha, rho = args

    m_n = 4.7583690772
    m_p = 4.7583690772

    sigma_sqrt, omega_sqrt, rho_03 = x

    sigma = sigma_sqrt**2
    omega = omega_sqrt**2

    m_eff = m_n - (g_sigma*sigma)

    # Proton
    rho_p    = alpha*rho
    kf_p     = (rho_p*(3*math.pi**2))**(1/3)
    E_fp     = (kf_p**2 + m_eff**2)**(1/2)
    rho_SB_p = (m_eff/(2.*math.pi**2))*(E_fp*kf_p - (m_eff**(2))*np.arctanh(kf_p/E_fp))

    #Neutron
    rho_n    = (1-alpha)*rho
    kf_n     = (rho_n*(3*math.pi**2))**(1/3)
    E_fn     = (kf_n**2 + m_eff**2)**(1/2)
    rho_SB_n = (m_eff/(2.*math.pi**2))*(E_fn*kf_n - (m_eff**(2))*np.arctanh(kf_n/E_fn))

    rho_b  = rho_p    + rho_n      
    rho_SB = rho_SB_p + rho_SB_n
    
    fvec_0 = (sigma*(m_sig**2)/g_sigma - rho_SB +(kappa*(g_sigma*sigma)**2)/2.\
              + (lambda_0*(g_sigma*sigma)**3)/6.)**2
    fvec_1 = (omega*(m_w**2)/g_omega - rho + (zeta*(g_omega*omega)**3)/6.\
              + 2.*Lambda_w*g_omega*omega*(rho_03*g_rho)**2)**2
    fvec_2 = (rho_03*(m_rho**2)/g_rho - (alpha-0.5)*rho + 2.*Lambda_w*g_rho*rho_03*(omega*g_omega)**2)**2

    f=[(fvec_0),(fvec_1),(fvec_2)]
    return f

@njit
def get_energy_pressure_alpha(x, rho, theta):
    """ Generate pressure and energy density at a given number density and proton fraction,
    
    (Faster Version Using Numba)
    
    Args:
        x (array): An array that consists of the initial values of sqrt(sigma), sqrt(omega), and rho 
        obtained from the initial_values function.
        rho (float): The central density from which the computation of the equation of state begins.
        theta (array): An array representing the parameters used to determine a RMF model in the
        Lagrangian. In this case, the RMF model is defined by 11 parameters, where the last parameters
        is the proton fraction (alpha).


    Returns:
        energy_density (float): EOS ingredient, energy density in natural units.
        pressure (float): EOS ingredient, pressure in natural units.

    """
    sigma_sqrt, omega_sqrt, rho_03 = x

    sigma = sigma_sqrt**2
    omega = omega_sqrt**2

    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w, alpha = theta
   
    
    m_n = 4.7583690772
    m_p = 4.7583690772

    m_eff    = m_n - (g_sigma*sigma)

   
    #Proton
    rho_p    = alpha*rho
    kf_p     = (rho_p*(3*math.pi**2))**(1/3)
    E_fp     = (kf_p**2 + m_eff**2)**(1/2)
    energy_p = (1/(8.*(math.pi**2)))*(kf_p*E_fp*(2*kf_p**2+m_eff**2) - np.log((kf_p + E_fp)/m_eff)*m_eff**4)
    integral_p = 1/4 * ( 1.5 * m_eff**4*np.arctanh(kf_p/E_fp) - 1.5*kf_p*m_eff**2*E_fp + kf_p**3*E_fp )
    Pressure_p = 1/3*(1/math.pi**2 * integral_p) 

    #Neutron
    rho_n    = (1-alpha)*rho
    kf_n     = (rho_n*(3*math.pi**2))**(1/3)
    E_fn     = (kf_n**2 + m_eff**2)**(1/2)
    energy_n = (1/(8.*(math.pi**2)))*(kf_n*E_fn*(2*kf_n**2+m_eff**2) - np.log((kf_n + E_fn)/m_eff)*m_eff**4)
    integral_n = 1/4 * ( 1.5 * m_eff**4*np.arctanh(kf_n/E_fn) - 1.5*kf_n*m_eff**2*E_fn + kf_n**3*E_fn )
    Pressure_n = 1/3*(1/math.pi**2 * integral_n) 

    #Total
    energy_b = energy_p + energy_n
    Pressure_b = Pressure_p + Pressure_n
    

    sigma_terms =  0.5*((sigma*m_sig)**2) + (kappa*((g_sigma*sigma)**3))/6.\
                    + (lambda_0*((g_sigma*sigma)**4))/24.

    # ε uses post-MF-substitution coefficients (ζ/8, 3·Λ_w); they absorb the matter-
    # meson coupling g_ω·ω·ρ_B via the ω-field equation of motion.
    omega_eps_terms = 0.5*((omega*m_w)**2) + (zeta*((g_omega*omega)**4))/8.
    rho_eps_terms   = 0.5*((rho_03*m_rho)**2) + 3.*Lambda_w*((g_rho*rho_03*g_omega*omega)**2)

    # P uses the raw Lagrangian coefficients (ζ/24, 1·Λ_w); P = T^{ii} carries no
    # matter-meson coupling to substitute.  The old code reused the ε coefficients,
    # which overstated P by (ζ/12)(g_ω ω)⁴ + 2·Λ_w (g_ρ ρ_03 g_ω ω)².
    omega_P_terms   = 0.5*((omega*m_w)**2) + (zeta*((g_omega*omega)**4))/24.
    rho_P_terms     = 0.5*((rho_03*m_rho)**2) + Lambda_w*((g_rho*rho_03*g_omega*omega)**2)

    energy_density = energy_b   + sigma_terms + omega_eps_terms + rho_eps_terms
    Pressure       = Pressure_p + Pressure_n  - sigma_terms + omega_P_terms   + rho_P_terms

    return energy_density, Pressure


def get_eos_alpha(theta, single_point = False):
    """ Generate EOS for a given alpha,
    
    (Faster Version Using Numba)
    
    Args:
        theta (array): An array representing the parameters used to determine a RMF model in the
        Lagrangian. In this case, the RMF model is defined by 11 parameters, where the last
        defined the proton fraction (alpha).
        single_point (boolean): Allows for the return of a single point of the EoS.

    Returns:
        rho (array): EOS ingredient, density in fm⁻³.
        energy_density (array): EOS ingredient, energy density in natural units.
        pressure (array): EOS ingredient, pressure in natural units.

    """
    if not single_point:
        x_init           = np.array(initial_guess_alpha(0.05, theta))
        dt               = 0.006
        N_points         = 125

        Density  = np.empty(N_points, dtype=float)
        Energy   = np.empty(N_points, dtype=float)
        Pressure = np.empty(N_points, dtype=float)
        for i in range(N_points):
            rho = 0.04 + i*dt

            arg    = np.append(theta, rho)
            sol    = optimize.root(fields_alpha, np.array(x_init).astype(np.float64) ,method='lm', args = arg.astype(np.float64))
            x_init = sol.x
            Re     = get_energy_pressure_alpha(x_init.astype(np.float64), rho, theta.astype(np.float64))

            Density[i]  = (rho)
            Energy[i]   = (Re[0]*oneoverfm_MeV/gcm3_to_MeVfm3)
            Pressure[i] = (Re[1]*oneoverfm_MeV/dyncm2_to_MeVfm3)

        rh = Density
        ep = Energy   * g_cm_3
        pr = Pressure * dyn_cm_2
        
        return rh, ep, pr
    else:
        rho    = single_point
        x_init = np.array(initial_guess_alpha(rho, theta))

        arg = np.append(theta, rho)
        sol = optimize.root(fields_alpha, np.array(x_init).astype(np.float64) ,method='lm', tol=1e-15, args = arg.astype(np.float64))
        Re  = get_energy_pressure_alpha(sol.x.astype(np.float64), rho, theta.astype(np.float64))

        return rho, Re[0] * oneoverfm_MeV / gcm3_to_MeVfm3 * g_cm_3 , Re[1] * oneoverfm_MeV / dyncm2_to_MeVfm3 *dyn_cm_2
