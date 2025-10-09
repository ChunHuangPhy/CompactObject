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
    """Right-hand side of the TOV + tidal ODE system.

    Args:
        r (float): raius as integrate varible
        y (psudo-varible): containing pressure, mass, h and b as intergarte varibles
        to solve out the TOV equation
        inveos: the invert of the eos, pressure and energy density relation to integrate
        and interpolate.

    Returns:
        Mass (array): The array that contains all the Stars' masses, in M_sun as a 
        Units.
        Radius (array): The array that contains all the Stars's radius, in km.
        Tidal Deformability (array): The array that contains correpsonding Tidal property, 
        These are dimension-less.
    """
    pres, m,h,b = y
    
    #energy_density = 10**inveos(np.log10(pres))
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


# Function solves the TOV equation, returning mass and radius
def solveTOV_tidal(center_rho, energy_density, pressure):
    """Solve TOV equation from given Equation of state in the neutron star core density range

    Args:
        center_rho(array): This is the energy density here is fixed in main that is np.logspace(14.3, 15.6, 50)
        energy_density (array): Desity array of the neutron star EoS, in MeV/fm^{-3}. Notice here for simiplicity, we omitted G/c**4 magnitude, so (value in MeV/fm^{-3})*G/c**4, could convert to the energy density we are using, please check the Test_EOS.csv to double check the order of magnitude.
        pressure (array): Pressure array of neutron star EoS, also in nautral unit with MeV/fm^{-3}, still please check the Test_EOS.csv, the conversion is (value in dyn/cm3)*G/c**4.

    Returns:
        Mass (array): The array that contains all the Stars' masses, in M_sun as a 
        Units.
        Radius (array): The array that contains all the Stars's radius, in km.
        Tidal Deformability (array): The array that contains correpsonding Tidal property, 
        These are dimension-less.
    """
    # eos = UnivariateSpline(np.log10(eps), np.log10(pres), k=1, s=0)
    # inveos = UnivariateSpline(np.log10(pressure), np.log10(energy_density), k=1, s=0)
    # We could change this to Double Log Interpolation。

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
        print(energy_density)
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

    r = 4.441e-16
    dr = 10.0
    rmax = 50 * 1e5
    rhocent = center_rho * G / c**2.0
    # pcent = 10**eos(np.log10(rhocent))
    pcent = eos(rhocent)
    P0 = pcent - (2.0 * pi / 3.0) * (pcent + rhocent) * (3.0 * pcent + rhocent) * r**2.0
    m0 = 4.0 / 3.0 * pi * rhocent * r**3.0
    h0 = r**2.0
    b0 = 2.0 * r
    stateTOV = np.array([P0, m0, h0, b0])
    ad_index = pressure_adind(P0, energy_density, pressure)
    sy = ode(TOV_def, None).set_integrator("dopri5")

    # have been modified from Irida to this integrator
    sy.set_initial_value(stateTOV, r).set_f_params(inveos, ad_index)

    while sy.successful() and stateTOV[0] > Pmin:
        stateTOV = sy.integrate(sy.t + dr)
        dpdr, dmdr, dhdr, dfdr = TOV_def(sy.t + dr, stateTOV, inveos, ad_index)
        dr = 0.46 * (1.0 / stateTOV[1] * dmdr - 1.0 / stateTOV[0] * dpdr) ** (-1.0)
    Mb = stateTOV[1]
    Rns = sy.t
    y = Rns * stateTOV[3] / stateTOV[2]
    tidal = tidal_deformability(y, Mb, Rns)

    return Mb * c**2.0 / G * unit.g, Rns * unit.cm, tidal


def solveTOV(center_rho, Pmin, eos, inveos):
    """Solve TOV equation from given Equation of state in the neutron star core density range

    Args:
        center_rho(float): The energy density in unit unit.g_cm_3
        Pmin (float): In unit.G / unit.c**4
        eos (function): pressure vs. energy density, energy density in unit.G / unit.c**2, pressure in unit.G / unit.c**4
        inveos (function): energy density vs. pressure

    Returns:
        Mass (array): The array that contains all the Stars' masses, in M_sun as a 
        Units.
        Radius (array): The array that contains all the Stars's radius, in km.
    """

    #eos = UnivariateSpline(np.log10(energy_density), np.log10(pres), k=1, s=0)
    #inveos = UnivariateSpline(np.log10(pres), np.log10(energy_density), k=1, s=0)
    #We could change this to Double Log Interpolation。
    c = 3e10
    G = 6.67428e-8
    Msun = 1.989e33

    eos = UnivariateSpline(energy_density, pressure, k=3, s=0)
    inveos = UnivariateSpline(pressure, energy_density, k=3, s=0)
    Pmin = pressure[20]
    r = 4.441e-16
    dr = 10.
    center_rho = center_rho * G/c**2.
    
    #pcent = 10**eos(np.log10(rhocent))
    pcent = eos(center_rho)
    P0 = (
        pcent
        - (4.0 * pi / 3.0) * (pcent + center_rho) * (3.0 * pcent + center_rho) * r**2.0
    )
    m0 = 4.0 / 3.0 * pi * center_rho * r**3.0
    stateTOV = np.array([P0, m0])

    sy = ode(TOV, None).set_integrator("dopri5")

    # have been modified from Irida to this integrator
    sy.set_initial_value(stateTOV, r).set_f_params(inveos)

    while sy.successful() and stateTOV[0] > Pmin:
        stateTOV = sy.integrate(sy.t + dr)
        dpdr, dmdr = TOV(sy.t + dr, stateTOV, inveos)
        dr = 0.46 * (1.0 / stateTOV[1] * dmdr - 1.0 / stateTOV[0] * dpdr) ** (-1.0)

    # at the end of this function, we rescale the quantities back
    return stateTOV[1] * unit.c**2 / unit.G, sy.t