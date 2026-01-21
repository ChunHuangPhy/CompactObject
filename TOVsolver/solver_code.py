# Python packages
import numpy as np
from scipy.constants import pi
from scipy.integrate import ode
from scipy.interpolate import interp1d
import TOVsolver.unit as unit


def m1_from_mc_m2(mc, m2):
    """a function that feed back the companion star mass from GW event measurement.

    Args:
        mc (float): chrip mass of a GW event, unit in solar mass.
        m2 (float or numpy array): the determined mass for one of the star, this is computed from sampling of EoS.

    Returns:
        m1 (array): the companion star mass in solar mass.
    """
    
    m2 = np.array(m2)
    num1 = (2.0 / 3.0) ** (1.0 / 3.0) * mc**5.0
    denom1 = (
        9 * m2**7.0 * mc**5.0
        + np.sqrt(3.0)
        * np.sqrt(abs(27 * m2**14.0 * mc**10.0 - 4.0 * m2**9.0 * mc**15.0))
    ) ** (1.0 / 3.0)
    denom2 = 2.0 ** (1.0 / 3.0) * 3.0 ** (2.0 / 3.0) * m2**3.0
    return num1 / denom1 + denom1 / denom2


def pressure_adind(P, epsgrid, presgrid):
    idx = np.searchsorted(presgrid, P)
    if idx == 0:
        eds = epsgrid[0] * np.power(P / presgrid[0], 3.0 / 5.0)
        adind = (
            5.0
            / 3.0
            * presgrid[0]
            * np.power(eds / epsgrid[0], 5.0 / 3.0)
            * 1.0
            / eds
            * (eds + P)
            / P
        )
    if idx == len(presgrid):
        eds = epsgrid[-1] * np.power(P / presgrid[-1], 3.0 / 5.0)
        adind = (
            5.0
            / 3.0
            * presgrid[-1]
            * np.power(eds / epsgrid[-1], 5.0 / 3.0)
            * 1.0
            / eds
            * (eds + P)
            / P
        )
    else:
        ci = np.log(presgrid[idx] / presgrid[idx - 1]) / np.log(
            epsgrid[idx] / epsgrid[idx - 1]
        )
        eds = epsgrid[idx - 1] * np.power(P / presgrid[idx - 1], 1.0 / ci)
        adind = (
            ci
            * presgrid[idx - 1]
            * np.power(eds / epsgrid[idx - 1], ci)
            * 1.0
            / eds
            * (eds + P)
            / P
        )
    return adind


def TOV(r, y, inveos):
    pres, m = y

    # eps = 10**inveos(np.log10(pres))
    eps = inveos(pres)
    dpdr = -(eps + pres) * (m + 4.0 * pi * r**3.0 * pres)
    dpdr = dpdr / (r * (r - 2.0 * m))
    dmdr = 4.0 * pi * r**2.0 * eps

    return np.array([dpdr, dmdr])


def TOV_def(r, y, inveos, ad_index):
    """Packs the complete set of TOV equations for neutron stars with tidal deformability

    Args:
        r (float): Radius as integration variable
        y (array): Array of integration variables containing [pressure, mass, h, b]
            where h and b are variables used for tidal deformability calculation
        inveos (function): Inverse equation of state function that maps pressure to energy density
        ad_index (float): Adiabatic index used in the tidal deformability equations

    Returns:
        array: Array containing derivatives [dpdr, dmdr, dhdr, dfdr] for the integration
    """
    
    pres, m, h, b = y

    # energy_density = 10**inveos(np.log10(pres))
    eps = inveos(pres)
    dpdr = -(eps + pres) * (m + 4.0 * pi * r**3.0 * pres)
    dpdr = dpdr / (r * (r - 2.0 * m))
    dmdr = 4.0 * pi * r**2.0 * eps
    dhdr = b

    dfdr = 2.0 * np.power(1.0 - 2.0 * m / r, -1) * h * (
        -2.0
        * np.pi
        * (5.0 * eps + 9.0 * pres + (eps + pres) ** 2.0 / (pres * ad_index))
        + 3.0 / np.power(r, 2)
        + 2.0
        * np.power(1.0 - 2.0 * m / r, -1)
        * np.power(m / np.power(r, 2) + 4.0 * np.pi * r * pres, 2)
    ) + 2.0 * b / r * np.power(1.0 - 2.0 * y[1] / r, -1) * (
        -1.0 + m / r + 2.0 * np.pi * np.power(r, 2) * (eps - pres)
    )

    return np.array([dpdr, dmdr, dhdr, dfdr])


def tidal_deformability(y2, Mns, Rns):
    """Compute Tidal deformability from y2, neutron star mass and raius

    Args:
        y2 (array): midiate varrible that computing tidal
        Mns (array): neutron star mass in g/cm3
        Rns (array): neutron star radius in cm.

    Returns:
        tidal_def (array): neutron star tidal deformability with unit-less.
    """
    C = Mns / Rns
    Eps = (
        4.0
        * C**3.0
        * (13.0 - 11.0 * y2 + C * (3.0 * y2 - 2.0) + 2.0 * C**2.0 * (1.0 + y2))
        + 3.0
        * (1.0 - 2.0 * C) ** 2.0
        * (2.0 - y2 + 2.0 * C * (y2 - 1.0))
        * np.log(1.0 - 2.0 * C)
        + 2.0 * C * (6.0 - 3.0 * y2 + 3.0 * C * (5.0 * y2 - 8.0))
    )
    tidal_def = (
        16.0 / (15.0 * Eps) * (1.0 - 2.0 * C) ** 2.0 * (2.0 + 2.0 * C * (y2 - 1.0) - y2)
    )

    return tidal_def


def solveTOV_tidal(center_rho, energy_density, pressure, max_radius=30e5*unit.cm):
    """Solve TOV equation from given Equation of state in the neutron star core density range

    Args:
        center_rho(array): This is the energy density here is fixed in main that is np.logspace(14.3, 15.6, 50)
        energy_density (array): Desity array of the neutron star EoS, in MeV/fm^{-3}. Notice here for simiplicity, we omitted G/c**4 magnitude, so (value in MeV/fm^{-3})*G/c**4, could convert to the energy density we are using, please check the Test_EOS.csv to double check the order of magnitude.
        pressure (array): Pressure array of neutron star EoS, also in nautral unit with MeV/fm^{-3}, still please check the Test_EOS.csv, the conversion is (value in dyn/cm3)*G/c**4.

    Returns:
        tuple[float, float, float]: (mass_Msun, radius_km, tidal_deformability)
            Tidal deformability is dimensionless.
    """
    # Unit conversions
    c = 3e10
    G = 6.67428e-8

    # tzzhou: migrating to new units
    center_rho = center_rho / unit.g_cm_3
    energy_density = energy_density * G / c**2 / unit.g_cm_3
    pressure = pressure * G / c**4 / unit.dyn_cm_2

    unique_pressure_indices = np.unique(pressure, return_index=True)[1]
    unique_pressure = pressure[np.sort(unique_pressure_indices)]

    if np.any(np.diff(unique_pressure) == 0):
        raise ValueError("Pressure values are not unique")

    if np.any(np.diff(energy_density) == 0):
        raise ValueError("energy_density values are not unique")

    # Interpolate pressure vs. energy density
    eos = interp1d(energy_density, pressure, kind="cubic", fill_value="extrapolate")

    # Interpolate energy density vs. pressure
    inveos = interp1d(
        unique_pressure,
        energy_density[unique_pressure_indices],
        kind="cubic",
        fill_value="extrapolate",
    )

    Pmin = pressure[20]

    r = 1e-18
    dr = 10.0
    rhocent = center_rho * G / c**2.0
    # pcent = 10**eos(np.log10(rhocent))
    pcent = eos(rhocent)
    P0 = pcent - (2.0 * pi / 3.0) * (pcent + rhocent) * (3.0 * pcent + rhocent) * r**2.0
    m0 = 4.0 / 3.0 * pi * rhocent * r**3.0
    h0 = r**2.0
    b0 = 2.0 * r
    stateTOV = np.array([P0, m0, h0, b0])
    ad_index = pressure_adind(P0, energy_density, pressure)
    
    # 1. Use a more stable integrator for stiff ODEs
    sy = ode(TOV_def, None).set_integrator("dopri5", atol=1e-8, rtol=1e-8)

    # 2. Set initial values
    sy.set_initial_value(stateTOV, r).set_f_params(inveos, ad_index)

    prev_mass = m0
    rho_current = rhocent
    
    # 3. Add more stopping conditions
    while (sy.successful() and 
           stateTOV[0] > Pmin and 
           sy.t < max_radius and
           rho_current > 1e-13 * rhocent):  # Stop if density drops too much
        
        stateTOV = sy.integrate(sy.t + dr)
        
        # TODO: This judgment related to rho_current seems unnecessary in most cases. 
        # Perhaps we can simplify it further.
        # 4. Calculate current density
        # rho_current = inveos(stateTOV[0])
        
        # 5. Sanity check for mass (it should never decrease)
        if stateTOV[1] < prev_mass:
            break
        prev_mass = stateTOV[1]
        
        # 6. Adaptive step size with safety factor and restrictions
        dpdr, dmdr, dhdr, dfdr = TOV_def(sy.t + dr, stateTOV, inveos, ad_index)
        
        delta= 0.23 ## optimal delta https://articles.adsabs.harvard.edu/pdf/1971ApJ...170..299B
        dr = delta / ((1.0 / stateTOV[1]) * dmdr - (1.0 / stateTOV[0]) * dpdr)    
        
        # 8. Reduce step size near surface where gradients are steep
        if stateTOV[0] < Pmin * 10:
            dr = min(dr, 50.0)

    # Clean up results - if we broke early for some reason
    if not sy.successful() or stateTOV[0] <= 0:
        # Handle the case where integration failed but we have previous valid values
        Mb = prev_mass
        Rns = sy.t
        # For tidal deformability, we need y, which depends on h and b 
        # If the last values of h and b are not valid, we can't compute y accurately
        # In this case, return a default tidal value (e.g., 0)
        tidal = 0  
    else:
        # Normal case - integration completed successfully
        Mb = stateTOV[1]
        Rns = sy.t
        y = Rns * stateTOV[3] / stateTOV[2]
        tidal = tidal_deformability(y, Mb, Rns)

    return Mb * c**2.0 / G * unit.g, Rns * unit.cm, tidal


def solveTOV(center_rho, Pmin, eos, inveos, max_radius=30e5*unit.cm):
    """Solve the Tolman-Oppenheimer-Volkoff (TOV) equation to determine the structure of a neutron star
    
    This function numerically integrates the TOV equations from the center outward to find the 
    mass and radius of a neutron star with a given central density and equation of state.
    The integration uses an adaptive step size method to ensure numerical stability,
    especially near the surface where pressure gradients become steep.
    
    Args:
        center_rho (float): The central energy density of the neutron star in unit.g_cm_3
        Pmin (float): Minimum pressure threshold that defines the star's surface, in unit.G / unit.c**4
        eos (function): Function mapping energy density to pressure (ρ → P)
                        Takes energy density in unit.G / unit.c**2, returns pressure in unit.G / unit.c**4
        inveos (function): Inverse equation of state, function mapping pressure to energy density (P → ρ)
        max_radius (float, optional): Safety parameter to prevent runaway integration in case
                                     the star's surface is not properly detected. Defaults to 30e5*unit.cm
    
    Returns:
        tuple[float, float]: (mass_Msun, radius_km)
            Mass is in solar masses; radius is in kilometers.
    
    Notes:
        - Integration begins at a small non-zero radius (r = 1e-18 * unit.cm) to avoid
          the coordinate singularity at r=0
        - Initial pressure is calculated using a Taylor expansion around r=0 since
          the TOV equation has an indeterminate form (0/0) at the origin
        - We use adaptive step sizing for numerical stability
        - Multiple stopping conditions ensure robustness against EoS-related numerical issues
    """
    
    # Initialize at a small radius to avoid singularity at r=0
    # Small enough for valid Taylor expansion but large enough to avoid floating-point issues
    r = 1e-18 * unit.cm  # Defined small length scale for computation
    dr = 10.0 * unit.cm  # Initial step size

    # Central pressure from equation of state
    pcent = eos(center_rho)
    
    # Calculate initial pressure using Taylor expansion around r=0
    # This approximation is derived from the TOV equation as r→0, avoiding the 0/0 indeterminate form
    # P(r) ≈ P(0) - (4π/3)(P+ρ)(3P+ρ)r² for small r
    P0 = (
        pcent
        - (4.0 * pi / 3.0) * (pcent + center_rho) * (3.0 * pcent + center_rho) * r**2.0
    )
    
    # Initial mass assuming uniform density within tiny sphere
    m0 = 4.0 / 3.0 * pi * center_rho * r**3.0
    
    # Initialize state vector [pressure, mass]
    stateTOV = np.array([P0, m0])
    
    # Use a more stable integrator for stiff ODEs with appropriate error tolerances
    sy = ode(TOV, None).set_integrator("dopri5", atol=1e-8, rtol=1e-8)
    
    # Set initial values and pass the inverse EoS function to the integrator
    sy.set_initial_value(stateTOV, r).set_f_params(inveos)
    
    # Variables to track integration progress and ensure stability
    prev_mass = m0
    rho_current = center_rho
    
    # Integration loop with multiple stopping conditions:
    # 1. Integration must be successful
    # 2. Pressure must remain above minimum threshold (defines stellar surface)
    # 3. Radius must not exceed maximum value (prevents runaway)
    # 4. Density must remain reasonably high (additional surface detection)
    
    while (sy.successful() and 
           stateTOV[0] > Pmin and 
           sy.t < max_radius and
           rho_current > 1e-13 * center_rho):  # Stop if density drops too much
        
        stateTOV = sy.integrate(sy.t + dr)
        
        # TODO: This judgment related to rho_current seems unnecessary in most cases. 
        # Perhaps we can simplify it further.
        # 4. Calculate current density
        # rho_current = inveos(stateTOV[0])
        
       
        # 5. Sanity check for mass (it should never decrease)
        if stateTOV[1] < prev_mass:
            break
        prev_mass = stateTOV[1]
        
        # 6. Adaptive step size with safety factor and restrictions
        dpdr, dmdr = TOV(sy.t + dr, stateTOV, inveos)
        delta= 0.23 ## optimal delta https://articles.adsabs.harvard.edu/pdf/1971ApJ...170..299B
        dr = delta / ((1.0 / stateTOV[1]) * dmdr - (1.0 / stateTOV[0]) * dpdr)         
        
        # 8. Reduce step size near surface where gradients are steep
        if stateTOV[0] < Pmin * 10:
            dr = min(dr, 50.0 * unit.cm)
    
    # Clean up results - if we broke early for some reason
    if not sy.successful() or stateTOV[0] <= 0:
        # Return the last valid values
        return prev_mass * unit.c**2 / unit.G, sy.t
    
    # at the end of this function, we rescale the quantities back
    return stateTOV[1] * unit.c**2 / unit.G, sy.t
