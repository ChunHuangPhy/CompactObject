import TOVsolver.constant as constant
import TOVsolver.solver_code as TOV_solver
import TOVsolver.EoS_import as EoS_import
import EOSgenerators.crust_EOS as crust

import TOVsolver.main as main
from scipy import interpolate
from TOVsolver.unit import  km, Msun, MeV,fm,g_cm_3,dyn_cm_2
import numpy as np
import math

oneoverfm_MeV = constant.oneoverfm_MeV
c = constant.c
G = constant.G

def MRlikihood_kernel(eps_total,pres_total,x,d1):
    """Computing likelihood from a distribution of MR measurement
    
    Args:
        eps_total (array): the energy density of full EoS in MeV/fm3, times a G/c**2 factor
        pres_total (array): the pressure from full EoS model in MeV/fm3, times a G/c**4 factor
        x (kde.kernel): the distribution kernel of MR measurement.
        d1 (float): the sampled density of this measurement

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """
    kernel = x
    MR = []
    if d1 ==0 :
        likelihood = -1e101
    else:
        d1 = 10**(d1)
        if   all(x<y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<y for x, y in zip(pres_total[:], pres_total[1:])):
            MR = main.OutputMR("",eps_total,pres_total,[d1*g_cm_3])[0]
        if len(MR) == False:
            likelihood = -1e101
        else:
            likelihood = np.log(kernel.evaluate((MR[1]/km, MR[0]/Msun)))
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood
    
def TidalLikihood_kernel(eps_total,pres_total,x,d1):
    """Computing likelihood from a distribution of Gravitational wave measurement
    
    Args:
        eps_total (array): the energy density of full EoS in MeV/fm3, times a G/c**2 factor
        pres_total (array): the pressure from full EoS model in MeV/fm3, times a G/c**4 factor
        x (kde.kernel): containing kernelGW and chirp, kernelGW is the distribution kde.kernel 
        sampled from full GW measurement, in [chrip mass, M2/M1, tidal of M1, tidal of M2] sequence.
        chrip mass is the sampling from chrip mass term in GW events solely.
        d1 (float): the sampled density of this measurement

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """
    kernelGW,chrip = x
    MRT = []
    chrip_mass = chrip.resample(1)
    if d1 ==0 :
        likelihood = -1e101
    else:
        d1 = 10**(d1)
        if   all(x<y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<y for x, y in zip(pres_total[:], pres_total[1:])):
            MRT = main.OutputMRTpoint(d1,eps_total,pres_total).T
            M1 = TOV_solver.m1_from_mc_m2(chrip_mass, MRT[1][0])
            Tidal_line = main.OutputMRT('',eps_total,pres_total).T
        if len(MRT[0]) == False or len(Tidal_line[0]) == False:
            likelihood = -1e101
        else:
            MTspline = interpolate.interp1d(Tidal_line[1]/Msun,Tidal_line[2])
            point = np.array([[chrip_mass[0][0]], [MRT[1][0] / M1[0][0]],[ MTspline(M1)[0][0]], [MRT[2][0]]])
            likelihood = np.log(kernelGW.evaluate(point))
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood

def MRlikihood_Gaussian(eps_total,pres_total,x,d1):
    """Computing likelihood from a simulation gaussian distribution of MR measurement
    
    Args:
        eps_total (array): the energy density of full EoS in MeV/fm3, times a G/c**2 factor
        pres_total (array): the pressure from full EoS model in MeV/fm3, times a G/c**4 factor
        x (float array): [Mvalue, Rvalue, Mwidth, Rwidth], Mvalue is the Mass center value of this 
        simulated measurement, Rvalue is the Radius center of it, Mwidth is the 1-sigma width of
        this Mass measurement, Rwidth is the 1-sigma width of this radius measurement.
        d1 (float): the sampled density of this measurement

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """

    Mvalue, Rvalue, Mwidth,Rwidth = x
    
    sigma_x = Rwidth
    sigma_y = Mwidth
    MR = []
    if d1 ==0:
        likelihood = -1e101
    else:
        d1 = 10**(d1)
        # if all(x<=y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<=y for x, y in zip(pres_total[:], pres_total[1:])):
        MR = main.OutputMR("",eps_total,pres_total,[d1*g_cm_3])[0]
        if len(MR) == False:
            likelihood = -1e101
        else:
            
            fx = 1/(sigma_x*sigma_y*(np.sqrt(2*np.pi))**2)*np.exp(-np.power(MR[1]/km-Rvalue, 2.)/(2*np.power(sigma_x,2.))-np.power(MR[0]/Msun-Mvalue, 2.)/(2*np.power(sigma_y,2.)))
            likelihood = np.log(fx)
            
            
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood
    
def Masslikihood_Gaussian(eps_total,pres_total,x,d1):
    """Computing likelihood from a simulation gaussian distribution of Mass measurement
    
    Args:
        eps_total (array): the energy density of full EoS in MeV/fm3, times a G/c**2 factor
        pres_total (array): the pressure from full EoS model in MeV/fm3, times a G/c**4 factor
        x (float array): [Mvalue, Mwidth], Mvalue is the Mass center value of this 
        simulated measurement, Mwidth is the 1-sigma width of this Mass measurement. 
        d1 (float): the sampled density of this measurement

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """
    Mvalue, Mwidth = x
    
    sigma_y = Mwidth
    MR = []
    if d1 ==0 :
        likelihood = -1e101
    else:
        d1 = 10**(d1)
        if   all(x<y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<y for x, y in zip(pres_total[:], pres_total[1:])):
            MR = main.OutputMR("",eps_total,pres_total,[d1*g_cm_3])[0]
        if len(MR) == False:
            likelihood = -1e101
        else:
            if MR[0]>= Mvalue:
                likelihood = 0
            else:    
                fx = 1/(sigma_y*(np.sqrt(2*np.pi))**2)*np.exp(-np.power(MR[0]/Msun-Mvalue, 2.)/(2*np.power(sigma_y,2.)))
                likelihood = np.log(fx)
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood


def Kliklihood(theta,K_low,K_up):
    """Computing likelihood from a hard cut constraint of K.
    
    Args:
    
        theta (array): An array representing the parameters used to determine a RMF model in the 
        Lagrangian. In this case, the RMF model is defined by 7 parameters. 
        K_low (float): lower bound of this K constraint.
        K_up (float): upper bound of this K constraint.

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, lambda_w = theta
    
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    
    m_b = 939
    E_per =-16.28
    m_sig =495/oneoverfm_MeV
    m_eff = (0.55+np.random.random_sample()*0.1)*m_b
    kappa_one = kappa*oneoverfm_MeV
    m_sig_one = m_sig*oneoverfm_MeV
    m_w_one = m_w*oneoverfm_MeV
    rho = 0.1505* (197.33**3 )
    k_f = (0.5*3.*rho*math.pi**2)**(1/3.)
    E_f = (k_f**2 + m_eff**2)**(1/2)

    g_w_omega = E_per - E_f + m_b
    sigma = m_b - m_eff
    m_sig_star_sqr_g_sig_sqr = (m_sig_one**2/g_sigma**2)+ (kappa_one*sigma) + 0.5*lambda_0*sigma**2
    m_w_star_sqr = m_w_one**2 + 0.5*zeta*(g_omega**2)*(g_w_omega**2)
    # Components o f rho_s_prime
    term1 = k_f/E_f
    term2 = E_f**2 + (2*m_eff**2)
    term3 = 3*m_eff**2
    term4 = np.log((k_f+E_f)/m_eff)
    rho_s_prime = (1/math.pi**2)*(term1*term2 - term3*term4)
    dM_drho = -(m_eff/E_f)*(1./(m_sig_star_sqr_g_sig_sqr+rho_s_prime))
    M_Ef = m_eff/E_f
    dEf_drho = math.pi**2/(2*k_f*E_f)
    dW0_drho = g_omega**2/m_w_star_sqr
    
    K = 9.*rho*(dEf_drho+dW0_drho+(M_Ef*dM_drho))
    center = (K_low + K_up)/2
    width = (K_up - K_low)/2
    p_K = -0.5*abs( center - K)**10./width**10.
    return p_K

def Jliklihood(theta,J_low,J_up):
    """Computing likelihood from a hard cut constraint of J.
    
    Args:
    
        theta (array): An array representing the parameters used to determine a RMF model in the 
        Lagrangian. In this case, the RMF model is defined by 7 parameters.
        K_low (float): lower bound of this J constraint.
        K_up (float): upper bound of this J constraint.

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w = theta
    
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    
    m_b = 939
    E_per =-16.28
    m_sig =495/oneoverfm_MeV
    m_eff = (0.55+np.random.random_sample()*0.1)*m_b
    kappa_one = kappa*oneoverfm_MeV
    m_sig_one = m_sig*oneoverfm_MeV
    m_w_one = m_w*oneoverfm_MeV
    rho = 0.1505* (197.33**3 )
    k_f = (0.5*3.*rho*math.pi**2)**(1/3.)
    E_f = (k_f**2 + m_eff**2)**(1/2)
    g_w_omega = E_per - E_f + m_b
    J_0 = k_f**2/(6*E_f)
    J_1 = (1/8.)*(rho*g_rho**2)/(m_rho**2+(2*Lambda_w*(g_w_omega**2)*g_rho**2))
    J = J_0 + J_1
    
    center = (J_low + J_up)/2
    width = (J_up - J_low)/2
    p_J = -0.5*abs( center - J)**10./width**10.
    return p_J

def Lliklihood(theta,L_low,L_up):
    """Computing likelihood from a hard cut constraint of L.
    
    Args:
    
        theta (array): An array representing the parameters used to determine a RMF model in the 
        Lagrangian. In this case, the RMF model is defined by 7 parameters. 
        K_low (float): lower bound of this L constraint.
        K_up (float): upper bound of this L constraint.

    Returns:
        likelihood (float): likelihood feed back for this given paramter set-up.
        
    """
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w = theta
    
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    
    m_b = 939
    E_per =-16.28
    m_sig =495/oneoverfm_MeV
    m_eff = (0.55+np.random.random_sample()*0.1)*m_b
    kappa_one = kappa*oneoverfm_MeV
    m_sig_one = m_sig*oneoverfm_MeV
    m_w_one = m_w*oneoverfm_MeV
    rho = 0.1505* (197.33**3 )
    k_f = (0.5*3.*rho*math.pi**2)**(1/3.)
    E_f = (k_f**2 + m_eff**2)**(1/2)
    J_0 = k_f**2/(6*E_f)
    J_1 = (1/8.)*(rho*g_rho**2)/(m_rho**2+(2*Lambda_w*(g_w_omega**2)*g_rho**2))
    J= J_0 + J_1 
    g_w_omega = E_per - E_f + m_b
    m_star_sqr_E_f_sqr = m_eff**2/E_f**2 
    rho_m_star = 3*rho/m_eff 
    sigma = m_b - m_eff
    m_sig_star_sqr_g_sig_sqr = (m_sig**2/g_sigma**2)+\
    (kappa*sigma) + 0.5*lambda_0*sigma**2
    
    # Components o f rho_s_prime
    term1 = k_f/E_f
    term2 = E_f**2 + (2*m_eff**2)
    term3 = 3*m_eff**2
    term4 = np.log((k_f+E_f)/m_eff)
    rho_s_prime = (1/math.pi**2)*(term1*term2 - term3*term4)
    dM_drho = -(m_eff/E_f)*(1./(m_sig_star_sqr_g_sig_sqr+rho_s_prime))
    m_w_star_sqr = m_w**2 + 0.5*zeta*(g_omega**2)*(g_w_omega**2)
    L_0 = J_0*(1+(m_star_sqr_E_f_sqr*(1. -(rho_m_star*dM_drho))))
    L_1 = 3.*J_1*(1.- 32.*(g_omega**2/m_w_star_sqr)*g_w_omega*Lambda_w*J_1)
    L= L_0 + L_1
    
    center = (L_low + L_up)/2
    width = (L_up - L_low)/2
    p_L = -0.5*abs( center - J)**10./width**10.
    return p_L



############################################ Date: 04 Nov 2024 #########################################################
def chiEFT_PNM( EoS_PNM, type="Gaussian", contraint_quantity="e", enlargement=0):
    """
    Authors: João Cartaxo, Tuhin Malik, Constança Providência
    
    Calculate the log-likelihood for the equation of state (EoS) of pure neutron matter using 
    chiEFT data extracted from chiral effective field theory. This can be achieved by assigning 
    constants to the energy per neutron ( E/N ) or the pressure ( p ), utilizing either a 
    Gaussian or Super-Gaussian likelihood model.
    
    Parameters
    ----------
    EoS_PNM : np.ndarray
        Array with PNM equation of state data, where:
        - EoS_PNM[0] is the density in fm^-3,
        - EoS_PNM[1] is the energy density MeV.fm^-3.
        - EoS_PNM[2] is the pressure in MeV.fm^-3.
        
    type : str, optional
        Specifies the type of distribution to use for the likelihood function.
        - "Gaussian" (default) or "Super Gaussian".
        
    contraint_quantity : str, optional
        Specifies which quantity (energy or pressure) to use for the log-likelihood calculation.
        - "e" for energy, "p" for pressure. The default is "e".
        
    enlargement : float, optional
        Enlargement factor (as a percentage in decimal form) for the Super-Gaussian distribution, e.g., 0.05 for 5%.
        Only applicable if type="Super Gaussian".
        
    Returns
    -------
    log_likelihood : float
    The sum of log-likelihoods over constraint on number density 0.08, 0.12 and 0.16 fm^-3.
    
    Explanation
    -----------
    - The likelihood calculation:
      - "Gaussian": Uses a standard Gaussian log-likelihood based on the discrepancy between EoS data and constraints.
      - "Super Gaussian": Employs an adjusted likelihood featuring a flattened peak of the Gaussian distribution.  
      
    Data Sources
    ------------
    - Energy per neutron constraints are taken from: Huth et al., Nature, vol 606, pp 276–280 (2022).
    - Pressure constraints are taken from: K. Hebeler et al., ApJ 773, 11 (2013).
    """
    def Gaussian(x, mu, sigma):
        """
        Compute the log of a Gaussian (Normal) probability density function (PDF).
        
        Parameters
        ----------
        x : float or np.ndarray
            The data point(s) at which to evaluate the Gaussian PDF.
        mu : float
            The mean (center) of the Gaussian distribution.
        sigma : float
            The standard deviation (spread) of the Gaussian distribution.
        
        Returns
        -------
        log_pdf : np.ndarray
            The log of the Gaussian PDF at each value of `x`.
        
        Explanation
        -----------
        - The Gaussian PDF is defined as:
          f(x) = (1 / sqrt(2 * pi * sigma^2)) * exp(-0.5 * ((x - mu) / sigma)^2)
        - We compute the log of this PDF for numerical stability:
          log(f(x)) = -0.5 * log(2 * pi * sigma^2) - 0.5 * ((x - mu) / sigma)^2
        """
        
        # Calculate the normalization term: -0.5 * log(2 * pi * sigma^2)
        log_of_norm = -0.5 * np.log(2 * np.pi * sigma**2)
        
        # Calculate the exponent term: -0.5 * ((x - mu) / sigma)^2
        log_of_exp = -0.5 * ((x - mu) / sigma) ** 2
        
        # Return the sum of both terms, representing the log of the Gaussian PDF
        return log_of_norm + log_of_exp

    def Super_Gaussian(x, mu, sigma, enlargement):
        """
        Compute the log of a Super-Gaussian probability density function (PDF).
        
        Parameters
        ----------
        x : float or np.ndarray
            The data point(s) at which to evaluate the Super-Gaussian PDF.
        mu : float
            The mean (center) of the Super-Gaussian distribution.
        sigma : float
            The standard deviation (spread) of the original Gaussian distribution.
        enlargement : float
            The percentage enlargement factor for `sigma` to increase the width of the distribution.
            For example, 0.05 for 5% enlargement, or 0.10 for 10% enlargement.
    
        Returns
        -------
        log_pdf : float or np.ndarray
            The log of the Super-Gaussian PDF at each value of `x`.
        
        Explanation
        -----------
        - The Super-Gaussian function extends the Gaussian by applying an enlargement factor to `sigma`,
          effectively "widening" the distribution.
        - Based on the reference: https://arxiv.org/pdf/2407.18452
        
        - The PDF calculation involves:
          - Adjusted variance term, `denom_1 = 2 * sigma_enlarged^2`
          - A damping factor, `denom_2 = 1 + exp((|x - mu| - sigma_enlarged) / 0.015)`, which modulates the tail behavior.
        """
        
        # Calculate the enlarged standard deviation by applying the percentage increase to sigma
        sigma_enlarged = sigma * (1 + enlargement)
    
        # Compute the first denominator term, incorporating the enlarged standard deviation
        denom_1 = 2 * sigma_enlarged ** 2
        
        # Compute the second denominator term, which dampens large deviations from the mean
        denom_2 = 1 + np.exp((np.abs(x - mu) - sigma_enlarged) / 0.015)
        
        # Calculate the Super-Gaussian PDF
        fx = 1 / (denom_1 * denom_2)
    
        log_fx = np.log(fx)
    
        return np.clip(log_fx, -1e20, np.inf)

    
    EoS_PNM = (EoS_PNM.T[EoS_PNM.T[:, 0] < 0.2]).T
    
    chidata_rho=[0.08,0.12,0.16]

    if contraint_quantity=="e":
        #### https://www.nature.com/articles/s41586-022-04750-w
        ### Nature volume 606, pages276–280 (2022), Huth et al
        ### Combined band including all the chiEFT sources in Fig. 4
        chidata_e     = np.array([8.90713476783698, 12.3640996602492, 16.1183465458664])  ## E/N : energy per neutron in MeV scale
        chidata_e_err = np.array([1.64779161947911, 2.36976217440553, 3.98357870894679])  ## error +- sigma
        real_e        = np.interp(chidata_rho, EoS_PNM[0], EoS_PNM[1]) / chidata_rho - 939.0
    
        if type=="Gaussian":
            return sum(Gaussian(real_e, chidata_e, chidata_e_err))
        elif type=="Super Gaussian":
            return sum(Super_Gaussian(real_e, chidata_e, chidata_e_err, enlargement))
    
    elif contraint_quantity=="p":
        ### K. Hebeler et al 2013 ApJ 773 11
        ## NN+3N data (Figure 2) 
        chidata_p     = np.array([0.505714285714279,1.24142857142857, 2.4857142857143])    ## in MeV.fm-3 units
        chidata_p_err = np.array([0.097142857142855,0.304285714285714, 0.691428571428572]) ## error +-sigma
        real_p = np.interp(chidata_rho, EoS_PNM[0], EoS_PNM[2])

        if type=="Gaussian":
            return sum(Gaussian(real_p, chidata_p, chidata_p_err))
        elif type=="Super Gaussian":
            return sum(Super_Gaussian(real_p, chidata_p, chidata_p_err, enlargement))
            
########################################################################################################################


############################################ Date: 05 Nov 2024 #########################################################
from InferenceWorkflow.pQCD import constraints
def ln_pQCD(EOS, rho_list=[0.92], points=1000):
    """
    Calculates the log-likelihood for the beta-equilibrium equation of state (EoS) using pQCD data. 
    
    Parameters
    ----------
    EoS : np.ndarray
        Array with beta-equilibrium equation of state data, where:
        - EoS[0] is the density in fm^-3,
        - EoS[1] is the energy density GeV.fm^-3.
        - EoS[2] is the pressure in GeV.fm^-3.
        
    rho_list : list, optional
       Specifies at which number densities (in fm^-3) to compute the pQCD contraint at.
        
    points : int, optional
        Number of points used to compute the weight, a higher number allows for greater precisions.
        
    Returns
    -------
    log_mean : float
        The averaged sum of log-likelihoods over the number densities specified
        by `rho_list`.

    References
    ----------
    Komoltsev and Kurkela (2021), Gorda et al. (2022), and Gorda et al.
    Phys. Rev. Lett. 127, 162003 (2021).
    """
    
    log_mean = 0
    for rho in rho_list:
        energy   = np.interp(rho, EOS[0], EOS[1])
        pressure = np.interp(rho, EOS[0], EOS[2])
        weight   = np.empty(points)

        for i in range(points):    
            X = np.exp(np.random.uniform( np.log(1), np.log(4) )) # Exp of the Log-linear distribution
            # To apply random weightage of the renormalization scale between X=1 to X=4            
            weight[i] = int(constraints(X, energy, pressure, rho)) 
        log_mean += np.log(weight.mean())
    log_mean = log_mean/len(rho_list)
                            
    return log_mean


########################################################################################################################











