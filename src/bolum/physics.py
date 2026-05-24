'''
bolum physics model file
'''

from bolum.constants import *

# estimate compressive strength based on density and porosity
def compute_strength(init_strength, porosity, density):

    if init_strength > 0:
        # use input material strength in Pa
        strength = init_strength
    else:
        # compute mat strength normalized to average mat strength
        # of ordinary chondrite and primitive chondrite (Brown et al 2002)
        strength = 1.e5 * (1. - porosity) * (density / 3000.)  # density/3000 term normalizes to silicate

    print(f"strength of bolide is {strength} Pa")

    return strength

# calculate area (assuming spherical)
def area_init_calc(diameter):

    # print out initial area
    area = 4. * pi * (diameter / 2.)**2.
    print(f"area of bolide = {area} m^2")

    return area

# calculate mass, given diameter and density
def mass_init_calc(diameter, density):

    # print out initial mass
    mass = (4. / 3.) * pi * (diameter / 2.)**3. * density
    print(f"mass of bolide = {mass} kg")

    return mass

# calculate temperature at altitude
def temp_altitude_calc(h):

    # altitude-dependent atmospheric model from NASA GRC, added 16 Oct 2025
    # converted Celsius temperature units to Kelvin by adding 273.15
    if(h >= 25000.):
        T = -131.21 + 0.00299 * h + 273.15
    if(h >= 11000. and h < 25000.):
        T = -56.46 + 273.15
    if(h < 11000.):
        T = 15.04 - 0.00649 * h + 273.15

    # Linear lapse rate model: T = T0 - L * h
    # T = T0 - L * h

    # T = max(T, 100)

    return T

# calculate atmospheric pressure
def atmos_pressure_calc(h, T):

    # altitude-dependent atmospheric model from NASA GRC
    # pressure here is in kPa
    if(h >= 25000.):
        p = 2.488 * (T / 216.6)**(-11.388)
    if(h >= 11000. and h < 25000.):
        p = 22.65 * math.exp(1.73 - 0.000157 * h)
    if(h < 11000.):
        p = 101.29 * (T / 288.08)**(5.256)
    
    # convert pressure to Pa 
    p = p * 1000.

    # standard exponential atmospheric model
    # p = p_air0 * np.exp(-h / H)

    return p

# calculate atmospheric density
def atmos_density_calc(h, p, T):

    # temperature- and pressure-dependent atmospheric model from NASA GRC
    # convert pressure bback to kPa first
    rho = (p / 1000.) / (0.2869 * T)

    # standard exponential atmospheric model
    # rho = rho_air0 * np.exp(-h / H)

    return rho

# calculate speed of sound
def sound_speed_calc(T, p, rho):

    # can use "right-hand side" of ideal gas equation
    # a = np.sqrt(gamma * R * T)  # R here is the specific gas constant (~287 J/(kg*K))

    # can also use "left-hand side" of ideal gas equation
    a = np.sqrt((gamma * p) / rho)

    return a

# calculate Mach number
def mach_num_calc(v, a):

    M = v / a

    return M

# calculate dynamic pressure
def dynamic_pressure_calc(rho, v):

    q = 0.5 * rho * v**2.

    return q

# calculate drag force
def Fdrag_calc(rho, v, area, Cd):

    # calculate cross-sectional area of a sphere
    # added 27 Oct 2025 CJM
    # area_cs = pi * r^2
    # if area = 4 * pi * r^2, then
    # area_cs = area / 4. is cross-sectional area
    area_cs = area / 4.

    # drag is defined as
    Fdrag = 0.5 * Cd * rho * v**2. * area_cs
    # units are kg/m^3 * m^2/s^2 * m^2 = kg * m/s^2 = N

    return Fdrag

# calculate kinetic energy
def KE_calc(mass, v):

    KE = 0.5 * mass * v**2.

    return KE

# calculate stagnation temperature
def T_stag_calc(M, T):

    T_stag = T * (1. + 0.5 * (gamma - 1.) * M**2.)

    return T_stag

# calculate stagnation pressure
def p_stag_calc(q, p):

    p_stag = q + p

    return p_stag

# calculate heat flux
def heatflux_calc(Fdrag, v, area):

    # cross-sectional area of a sphere
    # added 27 Oct 2025 CJM
    area_cs = area / 4.

    q_h = (Fdrag * v) / area_cs

    return q_h

# calculate surface temperature
def T_surf_calc(epsilon, q_h):

    T_surf = (q_h / (epsilon * sigma_sb))**0.25

    return T_surf

# calculate mass loss rate
def dmdt_calc(q_h, area, L_ablation):

    # units: q_h in W/m^2, area in m^2, and L_ablation in J/kg
    #        W/m^2 * m^2 / J/kg
    #        (J/s)/m^2 * m^2 / J/kg
    #        = kg/s
    # assumes spherical bolide in thermal equilibrium,
    # that is, mass is lost equally across surface area
    dmdt = (q_h * area) / L_ablation

    return dmdt

# function to update mass of bolide
def mass_update_calc(mass, dmdt, dt):

    mass = mass - (dmdt * dt)
    if mass < 0.:
        mass = 0.
    print(f"updated mass = {mass:.4e} kg")

    return mass

# function to calculate updated diameter
def diameter_update_calc(mass, density, n_fragments, fragmented):

    if fragmented == False:
        diameter = 2. * (((3. * mass) / (4. * pi * density))**(1. / 3.))
        print(f"updated diameter   = {diameter:.4e} meters")

    if fragmented == True:
        mass_fragment = mass / n_fragments
        diameter = 2. * (((3. * mass_fragment) / (4. * pi * density))**(1. / 3.))
        print(f"updated diameter   = {diameter:.4e} meters")

    return diameter

def area_update_calc(diameter, n_fragments, fragmented):

    if fragmented == False:
        area = 4. * pi * (diameter / 2.)**2.
        print(f"updated area       = {area:.4e} m^2")

    if fragmented == True:
        area_fragment = 4. * pi * (diameter / 2.)**2.
        area = area_fragment * n_fragments
        print(f"updated area       = {area:.4e} m^2")

    return area

# find intermediate total radiated energy
def E_rad_total_calc(E_rad_total, luminosity, dt):

    E_rad_total = E_rad_total + luminosity * dt
    print(f"total radiated energy is {E_rad_total:.3e} Joules")
    print(f"total radiated energy is {E_rad_total / jkt:.3e} kilotons")

    return E_rad_total

# find total deposited energy
def E_deposited_calc(E_deposited, Fdrag, v, dmdt, dt):

    E_dot_mech = Fdrag * v + 0.5 * v**2. * abs(dmdt)
    E_deposited = E_deposited + E_dot_mech * dt
    print(f"total deposited energy is {E_deposited:.3e} Joules")
    print(f"total deposited energy is {E_deposited / jkt:.3e} kilotons")

    return E_deposited

# find peak luminosity and peak radiant intensity
def peak_find(bolide_outputs, lum_tot_array):

    # find peak luminosity in W
    peak_power = np.max(lum_tot_array)
    print(f"bolide peak luminosity  = {peak_power:.3e} W")
    bolide_outputs.write(f"bolide peak luminosity = {peak_power:.3e} W\n")

    # find peak radiant intensity in W/sr
    peak_rad_i = peak_power / (4. * math.pi)
    print(f"bolide peak radiant intensity = {peak_rad_i:.3e} W/sr")
    bolide_outputs.write(f"bolide peak radiant intensity = {peak_rad_i:.3e} W\n")

    return peak_power, peak_rad_i

def E_event_calc(bolide_outputs, E_rad_tot, E_deposited):

    # from Brown et al, 2002, the totl event energy can be calculated by
    # E_event = 8.2508 * E_visible^0.885
    # this assumes a 6,000 K blackbody as the radiater (bolide)
    # for simplicity, since we only calculate the visible radiated energy
    # this can be updated in the future
    # # first, must convert total radiated energy into kt from J
    # E_rad_tot = E_rad_tot / jkt
    # # then calculate event energy
    # E_event = 8.2508 * (E_rad_tot)**0.885

    # # then convert back to joules
    # E_event = E_event * jkt

    # new, more direct way: calculate event energy from total E deposited
    E_event = E_deposited

    print(f"approximate total radiated energy of event = {E_rad_tot / jkt:.3e} kilotons")
    bolide_outputs.write(f"approximate total energy of event = {E_rad_tot / jkt:.3e} kilotons\n")
    
    print(f"approximate total energy of event = {E_event / jkt:.3e} kilotons")
    bolide_outputs.write(f"approximate total energy of event = {E_event / jkt:.3e} kilotons\n")

    return E_event