########################### Imports ############################
from TOVsolver.unit import g_cm_3, dyn_cm_2, km, Msun, MeV
from scipy import optimize
import pandas as pd
import sympy as sp
import numpy as np
import math


c = 3e10
G = 6.67428e-8
Msun = 1.989e33

dyncm2_to_MeVfm3 = 1.0 / (1.6022e33)
gcm3_to_MeVfm3 = 1.0 / (1.7827e12)
oneoverfm_MeV = 197.33

m_e = 2.5896 * 10**-3
m_mu = 0.53544
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



def Function(type='Typel99', couplings="Default"):
    """
    Defines the density-dependent couplings for the sigma, omega, and rho mesons 
    (g_sigma(rho), g_omega(rho), g_rho(rho)) based on various models. 

    The function returns symbolic expressions, allowing users to evaluate the couplings 
    and their derivatives at any given density rho. Coupling types and default values are 
    based on literature references, with options for user-defined functions.

    Parameters
    ----------
    type : str, optional
        Density-dependence model. Supported values are "Typel99", "Malik22",
        "Char23", and "UserDefined".
    couplings : str or list, optional
        Model-specific coupling constants. Use "Default" for the preset values
        of the selected model, or pass explicit coupling values.

    Returns
    -------
    tuple
        Six callable expressions: sigma, omega, and rho meson couplings,
        followed by their density derivatives.
    """

    if type == 'Typel99':
        """
        Ref. 2008.04491v1
        """
        if couplings == "Default":
            # DD-MEX model
            as_, av, ar, bs, bv, cs, cv, ds, dv, gs0, gv0, gr0, rho0 = [1.3970 , 1.3936 , 0.6202,
                                                                        1.3350 , 1.0191 , 
                                                                        2.0671 , 1.6060 , 
                                                                        0.4016 , 0.4556 , 
                                                                        10.7067, 13.3388, 7.2380, 0.153]
        else:
            as_, av, ar, bs, bv, br, cs, cv, cr, ds, dv, dr, gs0, gv0, gr0, rho0 = couplings

        x   = sp.symbols("x")
 
        gs  = sp.lambdify(x, sp.simplify(f"{gs0} * {as_} * (1+{bs}*(x/{rho0} + {ds})**(2))/(1+{cs}*(x/{rho0} + {ds})**(2))"))
        gw  = sp.lambdify(x, sp.simplify(f"{gv0} * {av}  * (1+{bv}*(x/{rho0} + {dv})**(2))/(1+{cv}*(x/{rho0} + {dv})**(2))"))
        gr  = sp.lambdify(x, sp.simplify(f"{gr0} * exp(-{ar}*(x/{rho0} - 1))"))

        dgs = sp.lambdify(x,sp.simplify(f"2*{gs0}*{as_}*(x/{rho0} + {ds})*({bs}-{cs})/(1+{cs}*(x/{rho0} + {ds})**(2))**(2) / {rho0}"))
        dgw = sp.lambdify(x,sp.simplify(f"2*{gv0}*{av} *(x/{rho0} + {dv})*({bv}-{cv})/(1+{cv}*(x/{rho0} + {dv})**(2))**(2) / {rho0}"))
        dgr = sp.lambdify(x,sp.simplify(f"- {ar} * ({gr0} * exp(-{ar}*(x/{rho0} - 1))) / {rho0} "))

        return gs, gw, gr, dgs, dgw, dgr

    elif type == 'Malik22':
        """
        https://doi.org/10.3847/1538-4357/ac5d3c
        """
        if couplings == "Default":
            # DDBm model
            as_, av, ar, gs0, gv0, grho0, rho0 = [0.086372, 0.054065, 0.509147, 9.180364, 10.981329, 3.826364*2, 0.150]
        else:
            as_, av, ar, gs0, gv0, grho0, rho0 = couplings

        x0   = sp.symbols("x0")

        gs  = sp.simplify(f"{gs0} * exp(-((x0/{rho0})**{as_} - 1.0))")
        gw  = sp.simplify(f"{gv0} * exp(-((x0/{rho0})**{av} - 1.0))")
        gr  = sp.simplify(f"{grho0} * exp(-{ar} * ((x0 / {rho0}) - 1.0))")

        dgs = sp.lambdify(x0,sp.simplify(f"-({gs} * {as_} / {rho0}) * (x0 / {rho0})**({as_} - 1.0)"))
        dgw = sp.lambdify(x0,sp.simplify(f"-({gw} * {av}  / {rho0}) * (x0 / {rho0})**({av} - 1.0)"))
        dgr = sp.lambdify(x0,sp.simplify(f"-({gr} * {ar}  / {rho0})"))

        gs = sp.lambdify(x0, gs)
        gw = sp.lambdify(x0, gw)
        gr = sp.lambdify(x0, gr)

        return gs, gw, gr, dgs, dgw, dgr

    elif type == 'Char23':
        """
        Functions used in: PhysRevD.108.103045 e-Print: 2307.12364 [nucl-th]
        """
        if couplings == "Default":
            # Model I
            as_, av, ar, bs, bv, br, cs, cv, cr, ds, dv, dr, rho0 = [8.225494 , 10.426752, 0.64584657,
                                                                     2.7079569, 1.6468675, 5.2033131 ,
                                                                     2.4776689, 6.8349408, 0.4262597 ,
                                                                     3.8630221, 1.4458185, -0.1824181, 0.16194209]
        else:
            as_, av, ar, bs, bv, br, cs, cv, cr, ds, dv, dr, rho0 = couplings

        x   = sp.symbols("x")
        n0  = 0.16   ## In calculation of Char23 the n0 plays as normalization factor and fixed to 0.16 fm-3
        
        gs_expr = sp.simplify(f"{as_} + ({bs} + {ds}*(x/{n0})**(3))*exp(-{cs}*x/{n0})")
        gw_expr = sp.simplify(f"{av}  + ({bv} + {dv}*(x/{n0})**(3))*exp(-{cv}*x/{n0})")
        gr_expr = sp.simplify(f"({ar}  + ({br} + {dr}*(x/{n0})**(3))*exp(-{cr}*x/{n0}))*2")   # Multiplied by 2 ( Lagrange density defintion is different )

        dgs = sp.lambdify(x, sp.diff(gs_expr, x))
        dgw = sp.lambdify(x, sp.diff(gw_expr, x))
        dgr = sp.lambdify(x, sp.diff(gr_expr, x))
        gs  = sp.lambdify(x, gs_expr)
        gw  = sp.lambdify(x, gw_expr)
        gr  = sp.lambdify(x, gr_expr)

        return gs, gw, gr, dgs, dgw, dgr
    
    elif type == "UserDefined":
        if couplings[-1] == "latex":
            from sympy.parsing.latex import parse_latex
            couplings = [parse_latex(couplings[i]) for i in range(3)]

        x   = sp.symbols("x")

        gs  = sp.simplify(couplings[0])
        dgs = sp.lambdify(x,sp.diff(gs, x))
        gs  = sp.lambdify(x, gs)

        gw  = sp.simplify(couplings[1])
        dgw = sp.lambdify(x,sp.diff(gw, x))
        gw  = sp.lambdify(x, gw)

        gr  = sp.simplify(couplings[2])
        dgr = sp.lambdify(x,sp.diff(gr, x))
        gr  = sp.lambdify(x, gr)

        return gs, gw, gr, dgs, dgw, dgr
    



############################### Beta Equilibrium ###############################
def initial_guess(rho, theta):
    """
    Computes initial values for the sigma, omega, and rho meson fields, 
    as well as the chemical potentials for neutrons and electrons at a given density.

    This method provides an initial approximation of these fields and potentials 
    based on the specified nuclear density and parameters of the chosen RMF model.

    Parameters:
        rho (float): Nuclear density at which the initial guess is calculated.
        theta (array): Array of parameters specific to the chosen RMF model. 
                       The number and meaning of these parameters vary by model.

    Returns:
        tuple:
            sigma (float): Initial value of the sigma meson field.
            omega (float): Initial value of the omega meson field.
            rho_03 (float): Initial value of the rho meson field.
            mu_n (float): Initial neutron chemical potential.
            mu_e (float): Initial electron chemical potential.
    """

    m_sig, m_w, m_rho, gsf, gwf, grf, dgsf, dgwf, dgrf, rho0 = theta

    g_sigma  = gsf(rho)
    g_omega  = gwf(rho)
    g_rho    = grf(rho)

    sigma  = g_sigma*rho/(m_sig**2)
    omega  = rho*((m_w**2)/g_omega)
    rho_03 = -g_rho*rho/(2.*(m_rho**2))

    m_eff_n = m_b[1]-(g_sigma*sigma) # Thesis Tiago 2.28  m_eff_Neutron
    mu_n    = m_eff_n + g_omega*omega + g_rho*rho_03*Matrix_b[1, 2]
    mu_e    = 0.12*m_l[0]*(rho/rho0)**(2/3.)

    return math.sqrt(sigma), math.sqrt(omega), rho_03, math.sqrt(mu_n), math.sqrt(mu_e)



def beta_equilibrium_function(x, args):
    """
    Iteratively adjusts the sigma, omega, and rho meson fields, as well as the chemical potentials 
    of neutrons and electrons, to achieve beta equilibrium at a specified density.

    Parameters:
        x (array): Initial values for sigma, omega, rho meson fields, and chemical potentials 
                   obtained from an initial guess function.
        args (array): Model parameters defining the RMF Lagrangian for the chosen RMF model.

    Returns:
        tuple:
            sigma (float): Adjusted value of the sigma meson field.
            omega (float): The adjusted value of the omega meson field.
            rho_03 (float): Adjusted value of the rho meson field.
            mu_n (float): Adjusted neutron chemical potential.
            mu_e (float): Adjusted electron chemical potential.
    """
    
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, dg_sigma, dg_omega, dg_rho, rho0, rho = args
    
    sigma_sqrt, omega_sqrt, rho_03, mu_n_sqrt, mu_e_sqrt = x
    
    sigma = sigma_sqrt**2
    omega = omega_sqrt**2
    mu_n  = mu_n_sqrt**2
    mu_e  = mu_e_sqrt**2
    
    rho_B_list  = []
    rho_l_list  = []
    rho_SB_list = []
    q_list      = []
    
    m_eff    = [m - (g_sigma*sigma) for m in m_b]
    
    Sigma_0R = dg_omega*omega*rho - dg_sigma*(sigma**2)*(m_sig**2)/g_sigma + dg_rho*(rho_03**2)*(m_rho**2)/g_rho
    
    for i in range(len(Matrix_b)):

        mu_b = Matrix_b[i,0]*mu_n - Matrix_b[i, 1]*mu_e
        
        E_fb = mu_b - g_omega*omega - g_rho*rho_03*Matrix_b[i,2] - Sigma_0R
        
        k_fb_sq = E_fb**2 - m_eff[i]**2
        if k_fb_sq <= 0:
            k_fb_sq = 0
            E_fb    = m_eff[i]
        
        k_fb = math.sqrt(k_fb_sq)
        
        rho_B  = k_fb**3 / (3.*math.pi**2)
        rho_SB = (m_eff[i]/(2.*math.pi**2))*(E_fb*k_fb - (m_eff[i]**2)*np.log((E_fb + k_fb )/m_eff[i]))
        
        rho_B_list.append(rho_B)
        rho_SB_list.append(rho_SB)
        
        Q_B = Matrix_b[i,1]*rho_B
        q_list.append(Q_B)
        
    for j in range(len(Matrix_l)):
        
        mu_l = Matrix_l[j,0]*mu_n - Matrix_l[j,1]*mu_e
        
        k_fl_sq = mu_l**2 - m_l[j]**2
        if k_fl_sq  < 0.0:
            k_fl_sq = 0.0
            mu_l    = m_l[j]
        k_fl = math.sqrt(k_fl_sq)
        
        rho_l = k_fl**3 / (3.*math.pi**2)
        Q_L   = Matrix_l[j,1]*rho_l
        q_list.append(Q_L)
        rho_l_list.append(rho_B)
    
    f = [(sigma*(m_sig**2)/g_sigma  - sum(rho_SB_list)),
        ( omega*(m_w**2)/g_omega    - sum(rho_B_list) ),
        ( rho_03*(m_rho**2)/g_rho   - sum([rho_B_list[k]*Matrix_b[k][2] for k in range(len(Matrix_b))])), 
        ( rho                       - sum(rho_B_list)),
        ( sum(q_list) )]

    return f

def get_energy_pressure(x, rho, theta):
    """
    Computes the energy density and pressure at a specified density, using given parameters 
    for an RMF model.

    Parameters:
        x (array): Array containing initial values of the sigma, omega, and rho meson fields, 
                   as well as the chemical potentials, obtained from an initial guess function.
        rho (float): Nuclear density at which to compute energy density and pressure.
        theta (array): Array of parameters defining the RMF model in the Lagrangian.

    Returns:
        tuple:
            energy_density (float): The energy density in natural units, essential for the EOS.
            pressure (float): The pressure in natural units, also essential for the EOS.
    """
    
    sigma_sqrt, omega_sqrt, rho_03, mu_n_sqrt, mu_e_sqrt = x
    sigma = sigma_sqrt**2
    omega = omega_sqrt**2
    mu_n  = mu_n_sqrt**2
    mu_e  = mu_e_sqrt**2
    
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, dg_sigma, dg_omega, dg_rho, rho0 = theta
    
    energy_b = 0
    energy_l = 0
    multi    = 0
    m_eff    = [m - (g_sigma*sigma) for m in m_b]

    rho_S       = sigma * m_sig**2 / g_sigma
    Sigma_0R    = dg_omega*omega*rho - dg_sigma*sigma*rho_S + dg_rho*(rho_03**2)*(m_rho**2)/g_rho
    Pressure_bl = 0
    q_list      = []
    for i in range(len(Matrix_b)):     
        
        mu_b = Matrix_b[i,0]*mu_n - Matrix_b[i, 1]*mu_e
        
        E_fb = mu_b - g_omega*omega - g_rho*rho_03*Matrix_b[i,2] - Sigma_0R
        
        k_fb_sq = E_fb**2 - m_eff[i]**2

        if k_fb_sq <= 0:
            k_fb_sq = 0.0
            E_fb = m_eff[i]
        
        k_fb = math.sqrt(k_fb_sq)

        rho_B = (k_fb**3) / (3.*math.pi**2)
        
        energy_baryon = (1/(8.*(math.pi**2)))*(k_fb*E_fb*(2*k_fb**2+m_eff[i]**2) - np.log((k_fb + E_fb)/m_eff[i])*m_eff[i]**4)

        energy_b +=  energy_baryon
        
        integral_b   = 1/4 * ( 1.5 * m_eff[i]**4*np.arctanh(k_fb/E_fb) - 1.5*k_fb*m_eff[i]**2*E_fb + k_fb**3*E_fb )
        Pressure_bl += 1/3*(1/math.pi**2 * integral_b)

        if i == 0:
            alpha = rho_B/rho
        
    for j in range(len(Matrix_l)):
        
        mu_l = Matrix_l[j, 0]*mu_n - Matrix_l[j, 1]*mu_e
        
        k_fl_sq = mu_l**2 - m_l[j]**2
        if k_fl_sq  < 0.0:
            k_fl_sq = 0.0
            mu_l    = m_l[j]
        k_fl = math.sqrt(k_fl_sq)
        
        energy_lepton = (1/(8.*(math.pi**2)))*(k_fl*mu_l*(2*k_fl**2+m_l[j]**2)-(m_l[j]**4)*np.log((k_fl+mu_l)/m_l[j]))
        energy_l     += energy_lepton
        integral_l    = 1/4 * ( 1.5 * m_l[j]**4*np.arctanh(k_fl/mu_l) - 1.5*k_fl*m_l[j]**2*mu_l + k_fl**3*mu_l )
        Pressure_bl  += 1/3*(1/math.pi**2) * integral_l
        
    sigma_terms = 0.5*((sigma*m_sig)**2)
    omega_terms = 0.5*((omega*m_w)**2)
    rho_terms   = 0.5*((rho_03*m_rho)**2)
    
    energy_density = energy_b + energy_l + sigma_terms + omega_terms + rho_terms
        
    Pressure = Pressure_bl - sigma_terms + omega_terms + rho_terms + Sigma_0R*rho
    return energy_density, Pressure, alpha

def compute_eos(eps_crust, pres_crust, theta):
    """
    Computes the core part of the equation of state (EOS) using a chosen RMF model, 
    complementing it with crust EOS data.

    Parameters:
        eps_crust (array): Energy density values for the crust EOS
        pres_crust (array): Pressure values for the crust EOS
        theta (array): Array of parameters defining the RMF model in the Lagrangian 
                       for the core EOS.

    Returns:
        tuple:
            energy_density (float): Energy density of the core in natural units, an essential EOS ingredient.
            pressure (float): Pressure of the core in natural units, is also essential for the EOS.
    """

    rho0 = theta[9]
    x_init           = np.array(initial_guess(0.04, theta))
    dt               = 0.006
    N_points         = 125

    Density  = np.empty(N_points, dtype=float)
    Energy   = np.empty(N_points, dtype=float)
    Pressure = np.empty(N_points, dtype=float)
    Alpha    = np.empty(N_points, dtype=float)
    for i in range(N_points):
        rho = 0.04 + i*dt

        theta_in    = theta.copy()
        theta_in[3] = theta_in[3](rho)
        theta_in[4] = theta_in[4](rho)
        theta_in[5] = theta_in[5](rho)
        theta_in[6] = theta_in[6](rho)
        theta_in[7] = theta_in[7](rho)
        theta_in[8] = theta_in[8](rho)
        theta_in    = np.array([float(t) for t in theta_in], dtype=np.float64)

        arg    = np.append(theta_in, float(rho))                
        sol    = optimize.root(beta_equilibrium_function, x_init ,method='lm', args = arg)
        if not sol.success:
            raise ValueError("Did not converge")
        x_init = sol.x
        Re     = get_energy_pressure(x_init, rho, theta_in)
    
        Energy[i]   = Re[0]
        Pressure[i] = Re[1]
        Alpha[i]    = Re[2]
        Density[i]  = rho

    end = 0
    for i in range(0, len(Energy) - 1):
        if Energy[i] > max(eps_crust) and i > 18:
            end = i + 2
            break
        end += 1
    
    Density  = Density[end::]
    Energy   = Energy[end::]
    Pressure = Pressure[end::]
    Alpha    = Alpha[end::]

    return Density, Energy, Pressure, Alpha





######################################  Alpha Depedent #######################3##########
def initial_guess_alpha(rho, theta):
    """
    Provides initial estimates for the sigma, omega, and rho meson fields at a specified 
    nuclear density, based on the parameters of a selected RMF model.

    Parameters:
        rho (float): Nuclear density at which to calculate initial field values.
        theta (array): Array of parameters defining the RMF model in the Lagrangian.

    Returns:
        tuple:
            sigma (float): Initial value of the sigma meson field.
            omega (float): Initial value of the omega meson field.
            rho_03 (float): Initial value of the rho meson field.
    """
    m_sig, m_w, m_rho, gsf, gwf, grf, dgsf, dgwf, dgrf, rho0, alpha  = theta
    g_sigma  = gsf(rho)
    g_omega  = gwf(rho)
    g_rho    = grf(rho)

    sigma  = g_sigma*rho/(m_sig**2)
    omega  = rho*((m_w**2)/g_omega)
    rho_03 = -g_rho*rho/(2.*(m_rho**2))
    
    return math.sqrt(sigma), math.sqrt(omega), rho_03

def fields_with_alpha(x, args):
    """
    Iteratively adjusts the sigma, omega, and rho meson fields at a given nuclear density, 
    based on a specified proton fraction, using the parameters of a selected RMF model.

    Parameters:
        x (array): Initial values of sigma, omega, rho meson fields, and chemical potentials 
                   from an initial guess function.
        args (array): Parameters defining the RMF model in the Lagrangian. The last element 
                      in this array specifies the proton fraction (alpha):
                      - alpha = 0 for pure neutron matter
                      - alpha = 0.5 for symmetric nuclear matter
                      - alpha = rho_p / (rho_p + rho_n)

    Returns:
        tuple:
            sigma (float): Adjusted value of the sigma meson field.
            omega (float): Adjusted value of the omega meson field.
            rho_03 (float): Adjusted value of the rho meson field.
    """

    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, dg_sigma, dg_omega, dg_rho, rho0, alpha, rho = args

    sigma_sqrt, omega_sqrt, rho_03 = x

    sigma = sigma_sqrt**2
    omega = omega_sqrt**2
    
    m_eff = m_b[1] - (g_sigma*sigma)

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

    f =[( rho                       - rho_b           ),
        ( sigma *(m_sig**2)/g_sigma - rho_SB          ),
        ( omega *(m_w**2)  /g_omega - rho             ),
        ( rho_03*(m_rho**2)/g_rho   - (alpha-0.5)*rho )]
    return f

def get_energy_pressure_alpha(x, rho, theta):
    """
    Computes the energy density and pressure for a specified nuclear density and proton fraction, 
    using parameters from a chosen RMF model.

    Parameters:
        x (array): Array containing initial values of sigma, omega, rho meson fields, and 
                   chemical potentials from an initial guess function.
        rho (float): Central density at which the equation of state (EOS) calculation begins.
        theta (array): Parameters defining the RMF model in the Lagrangian. The last element in 
                       this array specifies the proton fraction (alpha):
                       - alpha = 0 for pure neutron matter
                       - alpha = 0.5 for symmetric nuclear matter
                       - alpha = rho_p / (rho_p + rho_n)

    Returns:
        tuple:
            energy_density (float): Energy density in natural units, essential for the EOS.
            pressure (float): Pressure in natural units, is essential for the EOS.
    """
    
    sigma_sqrt, omega_sqrt, rho_03 = x

    sigma = sigma_sqrt**2
    omega = omega_sqrt**2
    
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, dg_sigma, dg_omega, dg_rho, rho0, alpha = theta
    
    m_eff    = m_b[1] - (g_sigma*sigma)

    rho_S    = sigma * m_sig**2 / g_sigma
    Sigma_0R = dg_omega*omega*rho - dg_sigma*sigma*rho_S + dg_rho*(rho_03**2)*(m_rho**2)/g_rho

    #Proton
    rho_p    = alpha*rho
    kf_p     = (rho_p*(3*math.pi**2))**(1/3)
    E_fp     = (kf_p**2 + m_eff**2)**(1/2)
    energy_p = (1/(8.*(math.pi**2)))*(kf_p*E_fp*(2*kf_p**2+m_eff**2) - np.log((kf_p + E_fp)/m_eff)*m_eff**4)
    integral_p = 1/4 * ( 1.5 * m_eff**4*np.arctanh(kf_p/E_fp) - 1.5*kf_p*m_eff**2*E_fp + kf_p**3*E_fp )
    Pressure_p = 1/3*(1/math.pi**2 * integral_p) +  Sigma_0R*rho_p
    
    #Neutron
    rho_n    = (1-alpha)*rho
    kf_n     = (rho_n*(3*math.pi**2))**(1/3)
    E_fn     = (kf_n**2 + m_eff**2)**(1/2)
    energy_n = (1/(8.*(math.pi**2)))*(kf_n*E_fn*(2*kf_n**2+m_eff**2) - np.log((kf_n + E_fn)/m_eff)*m_eff**4)
    integral_n = 1/4 * ( 1.5 * m_eff**4*np.arctanh(kf_n/E_fn) - 1.5*kf_n*m_eff**2*E_fn + kf_n**3*E_fn )
    Pressure_n = 1/3*(1/math.pi**2 * integral_n) +  Sigma_0R*rho_n

    #Total
    energy_b = energy_p + energy_n

    sigma_terms = 0.5*((sigma*m_sig)**2)
    omega_terms = 0.5*((omega*m_w)**2)
    rho_terms   = 0.5*((rho_03*m_rho)**2)
    
    energy_density = energy_b   + sigma_terms + omega_terms + rho_terms
    Pressure       = Pressure_p + Pressure_n  - sigma_terms + omega_terms + rho_terms

    return energy_density, Pressure

def compute_eos_alpha(theta):
    """
    Generates the equation of state (EOS) table for a sequence of densities, considering a specified 
    proton fraction, using parameters from the chosen RMF model.

    Parameters:
        eps_crust (array): Energy density values for the crust EOS, in MeV/fm³ (including a G/c² factor).
        pres_crust (array): Pressure values for the crust EOS, in MeV/fm³ (including a G/c⁴ factor).
        theta (array): Parameters defining the RMF model in the Lagrangian. The last element specifies 
                       the proton fraction (alpha):
                       - alpha = 0 for pure neutron matter
                       - alpha = 0.5 for symmetric nuclear matter
                       - alpha = rho_p / (rho_p + rho_n)

    Returns:
        tuple:
            energy_density (float): Energy density in natural units, an essential EOS ingredient.
            pressure (float): Pressure in natural units, also essential for the EOS.
    """
    rho0 = theta[9]
    x_init           = np.array(initial_guess_alpha(0.04, theta))
    dt               = 0.006
    N_points         = 125

    Density  = np.empty(N_points, dtype=float)
    Energy   = np.empty(N_points, dtype=float)
    Pressure = np.empty(N_points, dtype=float)
    
    for i in range(N_points):
        rho = 0.04 + i*dt

        theta_in    = theta.copy()
        theta_in[3] = theta_in[3](rho)
        theta_in[4] = theta_in[4](rho)
        theta_in[5] = theta_in[5](rho)
        theta_in[6] = theta_in[6](rho)
        theta_in[7] = theta_in[7](rho)
        theta_in[8] = theta_in[8](rho)
        theta_in    = np.array([float(t) for t in theta_in], dtype=np.float64)

        arg    = np.append(theta_in, float(rho)) 
        sol    = optimize.root(fields_with_alpha, x_init ,method='lm', args = arg)
        x_init = sol.x
        Re     = get_energy_pressure_alpha(x_init, rho, theta_in)

        Density[i]  = rho
        Energy[i]   = Re[0]
        Pressure[i] = Re[1]
    
    return Density, Energy, Pressure
