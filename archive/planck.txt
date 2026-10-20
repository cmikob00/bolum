import math
import numpy as np
from scipy.integrate import quad
from scipy import constants

# Get fundamental physical constants
h = constants.h  # Planck's constant (J s)
c = constants.speed_of_light  # Speed of light (m/s)
k = constants.k  # Boltzmann constant (J/K)
sigma = constants.sigma # Stefan-Boltzmann constant (W/m^2/K^4)

def planck_wavelength_radiance(wavelength, T):
    """
    Calculates the spectral radiance according to Planck's Law (wavelength form).

    Parameters:
    wavelength (float or array-like): Wavelength in meters.
    T (float): Temperature in Kelvin.

    Returns:
    float or array-like: Spectral radiance in W/(m^2 sr m).
    """
    if wavelength is None:
        return 0
    # The form of Planck's law used here is for spectral radiance in W/(m^2 sr m)
    radiance = (2.0 * h * c**2) / (wavelength**5 * (np.exp(h * c / (wavelength * k * T)) - 1.0))
    return radiance

def integrate_planck_range_nm(T, min_wavelength_nm, max_wavelength_nm):
    """
    Integrates Planck's law over a specific wavelength range in nanometers.

    Parameters:
    T (float): Temperature in Kelvin.
    min_wavelength_nm (float): Minimum wavelength in nanometers.
    max_wavelength_nm (float): Maximum wavelength in nanometers.

    Returns:
    float: Total radiance in W/(m^2 sr) within the specified range.
    """
    # Convert nanometers to meters for the Planck function
    min_wavelength_m = min_wavelength_nm * 1e-9
    max_wavelength_m = max_wavelength_nm * 1e-9

    # Use scipy.integrate.quad for numerical integration
    # The result is the total radiance (W/(m^2 sr)) within the range
    # args=(T,)) is a tuple containing a single element
    result, error = quad(planck_wavelength_radiance, min_wavelength_m, max_wavelength_m, args=(T,))

    return result

# --- Example Usage ---
temperature_k = 5800  # Sun's surface temperature in Kelvin
min_range_nm = 400   # Start of visible light range (blue)
max_range_nm = 700   # End of visible light range (red)

radiance_in_range = integrate_planck_range_nm(temperature_k, min_range_nm, max_range_nm)

print(f"Temperature: {temperature_k} K")
print(f"Wavelength Range: {min_range_nm} nm to {max_range_nm} nm")
print(f"Total Radiance in range: {radiance_in_range:.4e} W/(m^2 sr)")

# convert to watts
# P (W) = radiance (W/m^2*sr) * area (m^2) * 4*pi (sr)
area = 10.  # m^2, example
P = radiance_in_range * area * (4. * math.pi)
