# Python packages
import numpy as np
from scipy.interpolate import interp1d

# Import files
import TOVsolver.solver_code as TOV_solver
import TOVsolver.EoS_import as EoS_import
import TOVsolver.speed_of_sound as speed_of_sound
from TOVsolver.unit import g_cm_3, G, c, km, Msun


def _has_physical_mr_points(mass, radius):
    """Return True when an MR sequence contains enough plausible NS points."""
    if len(mass) < 3 or len(radius) < 3:
        return False

    mass = np.asarray(mass, dtype=float)
    radius = np.asarray(radius, dtype=float)
    radius_km = radius / km
    mass_msun = mass / Msun

    physical = (
        np.isfinite(radius_km)
        & np.isfinite(mass_msun)
        & (radius_km > 4.0)
        & (radius_km < 40.0)
        & (mass_msun > 0.02)
        & (mass_msun < 4.0)
    )
    return np.count_nonzero(physical) >= 3

# Global Variables
def OutputMR(input_file="", density=[], pressure=[], central_density_range=np.logspace(14.3, 15.6, 50) * g_cm_3):
    """Outputs the mass, radius, and tidal deformability
    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        MR (tuple): tuple with mass, radius.
    """

    #############This is something we need to change, like the input for this EOS import should
    ############# be one file contatining Whole EOS. that first column is density and second is pressure
    energy_density, pressure = EoS_import.EOS_import(input_file, density, pressure)
    original_energy_density = energy_density.copy()
    original_pressure = pressure.copy()
    ############# Lets the user only input the EOS file path, then this EOS_import should have file
    ############# as input. and the outputMR should have a file as input too?

    # Notice that we only rescale quantities inside this function
    energy_density = energy_density * G / c**2
    pressure = pressure * G / c**4

    # eos = UnivariateSpline(np.log10(energy_density), np.log10(pres), k=1, s=0)
    # inveos = UnivariateSpline(np.log10(pres), np.log10(energy_density), k=1, s=0)
    # We could change this to Double Log Interpolation。

    unique_pressure_indices = np.unique(pressure, return_index=True)[1]
    unique_pressure = pressure[np.sort(unique_pressure_indices)]

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

    Radius = []
    Mass = []

    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    for i in range(len(central_density_range)):
        try:
            center_rho = central_density_range[i] * G / c**2
            M, R = TOV_solver.solveTOV(center_rho, Pmin, eos, inveos)
            Mass.append(M)
            Radius.append(R)
        # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
        except OverflowError as e:
            # print("This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated.")
            break

    if not _has_physical_mr_points(Mass, Radius):
        Radius = []
        Mass = []
        for center_density in central_density_range:
            try:
                M, R, _ = TOV_solver.solveTOV_tidal(
                    center_density,
                    original_energy_density,
                    original_pressure,
                )
                Mass.append(M)
                Radius.append(R)
            except OverflowError:
                break

    MR = np.vstack((Mass, Radius)).T
    # print("Mass Radius file will be generated and stored as  2-d array. The first column is mass, second one is Radius")

    return MR


def OutputMRT(input_file="", density=[], pressure=[]):
    """Outputs the mass, radius, and tidal deformability
    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        MRT (tuple): tuple with mass, radius, and tidal deformability.
    """

    #############This is something we need to change, like the input for this EOS import should
    ############# be one file contatining Whole EOS. that first column is density and second is pressure
    energy_density, pressure = EoS_import.EOS_import(input_file, density, pressure)
    ############# Lets the user only input the EOS file path, then this EOS_import should have file
    ############# as input. and the outputMR should have a file as input too?

    Radius = []
    Mass = []
    tidal = []
    density = np.logspace(14.3, 15.6, 50) * g_cm_3
    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    for i in range(len(density)):
        try:
            M, R, T = TOV_solver.solveTOV_tidal(density[i], energy_density, pressure)
            Radius.append(R)
            Mass.append(M)
            tidal.append(T)
        # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
        except OverflowError as e:
            # print("This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated.")
            break
    MRT = np.vstack((Radius, Mass, tidal)).T
    # print("Mass Radius and tidal will be generated as the 3-d array. The first column is Radius, second one is mass,last is tidal")

    return MRT


def OutputC_s(input_file="", density=[], pressure=[]):
    """Calls function to open csv (if needed) and check equation of state validity.
        Then calls function to calculate speed of sound.

    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        C_s (array): numpy 1D array. List of speeds of sound.
    """

    energy_density, pressure = EoS_import.EOS_import(input_file, density, pressure)
    C_s = speed_of_sound.speed_of_sound_calc(energy_density, pressure)
    return C_s


def OutputMRTpoint(central_density, energy_density, pressure):
    """Outputs the mass, radius, and tidal deformability (single point)
    Args:
        central_density (float): central density that we want to compute
        density (array, optional): numpy 1Darray. Density of EoS
        pressure (array, optional): numpy 1Darray. pressure of EoS

    Returns:
        MRT (tuple): tuple with mass, radius and tidal.
    """
    Radius = []
    Mass = []
    tidal = []
    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    try:
        M, R, T = TOV_solver.solveTOV_tidal(central_density* g_cm_3, energy_density, pressure)
        Mass.append(M)
        Radius.append(R)
        tidal.append(T)
    # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
    except OverflowError as e:
        print(
            "This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated."
        )
    MRT = np.vstack((Radius, Mass, tidal)).T

    return MRT
