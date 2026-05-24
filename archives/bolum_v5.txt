# Bolide Luminosity Model (bolum) Version 5
# C.J. Miko
# Instructions: Modify inputs in the input deck (bolum_in.txt) and run to see if/when your bolide fragments
#               and the resultant light curve.
# Current line count: 1446 lines
#
# Recent Additions:
#   October 2025 - v5:
#   improved plotting and visualization
#   added altitude-dependent atmospheric temperature, pressure, and density model from NASA GRC
#   added writing out atmospheric temperature, pressure, and density to output file
#   fixed Mach scaling term in shock luminosity calculation (lum_shock_calc)
#   combined Mach scaling and density scaling term into one term called Shock Intensity Scaling Factor (SISF)
#   added total radiated energy in kilotons to output file
#   output quantities now written directly to numpy arrays instead of lists: started with np.zeros arrays,
#   then append with a new row after each timestep
#   added timestep counter (helps with array indexing)
#   added total event energy calculation from Brown et al, 2002
#   added max number of cycles parameter (100,000)
#   improved physics and logic for fragmentation: now accounting for total surface area of n_fragments
#   and no artificially-modeled "flare"
#   updates to diameter and area update functions to account for tracking multiple fragments
#   updates to drag force calculations and heat flux to account for cross-sectional area
#   consolidated fragmentation logic and math into one self-contained function
#   added ballistic trajectory functionality - seems to work fine
#
#   July 2025 - v4:
#   added scaled luminous efficiency function based on velocity and density (Ceplecha et al 1998)
#   added print statement to tell when the bolide fragments
#   added plotting output vs time in W/sr and output vs altitude in W/sr
#   changed angle variable in input deck to be angle_deg for consistency
#   added in x and y components for position, velocity, and acceleration due to gravity
#   added ballistic coefficient (default 1.5) for drag calculcation
#   added shock front luminosity calculations - includes scaling for density and stagnation pressure
#   made reading in input deck function
#   added peak power, radiant intensity, and total energy radiated finding functions
#
#   Room for Improvement:
#   need to consider adding in true fragment tracking (separate module currently)
#   include radiation transport to a fixed sensor in the far field to get power/irrad at a distance

# Library Imports
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties

# Constants
pi                = math.pi             # numerical value of pi
g                 = 9.80665             # gravity (m/s^2)
R                 = 287.05              # J/kg·K, specific gas constant for dry air
sigma_sb          = 5.670374419e-8      # Stefan–Boltzmann constant
rad2deg           = 180. / math.pi      # conversion factor from radians to degrees
deg2rad           = math.pi / 180.      # conversion factor from degrees to radians
jkt               = 4.185e12            # conversion factor from joules to kilotons
R_Earth           = 6378.14e3           # Earth's radius in meters
C_Earth           = 2. * pi * R_Earth   # Earth's circumference in meters

# air constants
gamma             = 1.4             # adiabatic index for air
rho_air0          = 1.225           # sea-level atmospheric density (kg/m^3)
p_air0            = 101325.         # sea-level air pressure in Pa
H                 = 7160.0          # scale height of atmosphere (m)
T0                = 288.15          # temperature at sea-level
L                 = 0.0065          # temperature lapse rate in K / m
rho_cutoff        = 0.01            # atmospheric density cutoff for strong shock formation

# bolide constants
lum_efficiency    = 0.03            # baseline luminous efficiency
strength_baseline = 1.e6            # baseline strength in Pascals
epsilon           = 0.9             # surface emissivity
L_ablation        = 5.e6            # heat of ablation in J/kg
eta_0             = 0.03            # luminous efficiency scaling factor
Cd                = 1.5             # ballistic coefficient for hypersonic flows
v_ref             = 20000.          # m/s, reference velocity for luminous efficiency scaling
n                 = 2.              # velocity exponent for luminous efficiency scaling (Ceplecha et al 1998)
m                 = 0.5             # density exponent for luminous efficiency scaling (Ceplecha et al 1998)
k_shock           = 2.              # bolide effective area factor for shock luminosity calculations
alpha             = 0.85            # smoothing weight for mediating flares
luminosity_prev   = 0.              # dummy variable to store luminosity for smooth scaling

# other simulation parameters
mxcycl            = 100000          # max number of cycles
new_row           = np.zeros(1)     # new row addition to numpy arrays at end of timestep
new_col           = np.zeros(1)     # new column addition to numpy arrays at end of timestep


# I/O functions

# reading in input deck
def rdinput(inputpath):
    params = {}
    with open(inputpath, "r") as f:
        lines = f.readlines()[2:]  # skip first two lines

    for line in lines:
        if "=" not in line:
            continue
        key, rest = line.split("=", 1)
        key = key.strip()
        # remove comments
        value_str = rest.split("#")[0].strip()

        # cast to int if whole number, otherwise float
        try:
            value = float(value_str)
        except ValueError:
            continue

        params[key] = value

    # assign to variables
    diameter       = params["diameter"]
    velocity       = params["velocity"]
    theta_deg      = params["theta_deg"]
    init_strength  = params["init_strength"]
    density        = params["density"]
    porosity       = params["porosity"]
    altstart       = params["altstart"]
    xstart         = params["xstart"]
    tstart         = params["tstart"]
    dt             = params["dt"]
    tstop          = params["tstop"]
    n_fragments    = params["n_fragments"]
    SISF_flag      = params["SISF_flag"]

    return diameter, velocity, theta_deg, density, porosity, init_strength, altstart, xstart, tstart, dt, tstop, n_fragments, SISF_flag

# writing initial parameters to output file
def write_init_params(diameter, velocity, theta_deg, density, porosity, init_strength, altstart, xstart, tstart, dt, tstop, n_fragments):

    # open output file
    init_params = open('init_params.txt', 'w')

    init_params.write(f"Initialization Parameters for Bolide Luminosity and Fragmentation Simulation\n")
    init_params.write(f"Diameter (m)                    {diameter:.3f}\n")
    init_params.write(f"velocity (m/s)                  {velocity:.3f}\n")
    init_params.write(f"Entry Angle to Horizon (deg)    {theta_deg:.3f}\n")
    init_params.write(f"Density (kg/m^3)                {density:.3f}\n")
    init_params.write(f"Porosity                        {porosity:.3f}\n")
    init_params.write(f"Material Strength (Pa)          {init_strength:.3f}\n")
    init_params.write(f"Initial Altitude (m)            {altstart:.3f}\n")
    init_params.write(f"Initial Position (m)            {xstart:.3f}\n")
    init_params.write(f"Starting Time (s)               {tstart:.3f}\n")
    init_params.write(f"Timestep (s)                    {dt:.3f}\n")
    init_params.write(f"Maximum Time (s)                {tstop:.3f}\n")
    init_params.write(f"Number of Fragments             {n_fragments}\n")

    init_params.close()

def wr_out(bolide_outputs, luminosity, y, t, mass, diameter, area, KE, v, acc, q_h, Fdrag, E_rad_total, q, M, T_stag, p_stag, T_surf, T, p, rho, a, SISF, theta_deg):
    
    # writing outputs after timestep
    bolide_outputs.write(f"time = {t:.4e}   lum = {luminosity:.4e}   alt = {y:.4e}     mass = {mass:.4e}      dia = {diameter:.4e}     area = {area:.4e}\n")
    bolide_outputs.write(f"KE = {KE:.4e}     v = {v:.4e}     acc = {acc:.4e}     E_rad_ttl = {E_rad_total:.4e} Fdrag = {Fdrag:.4e}   theta_deg = {theta_deg:.3f}\n")
    bolide_outputs.write(f"q_h = {q_h:.4e}    q = {q:.4e}     M = {M:.4e}       T_stag = {T_stag:.4e}    p_stag = {p_stag:.4e}  T_surf = {T_surf:.4e}\n")
    bolide_outputs.write(f"T = {T:.4e}      p = {p:.4e}     rho = {rho:.4e}     a = {a:.4e}         SISF = {SISF:.4e}\n")
    bolide_outputs.write(f"\n")

# plotting outputs function
def plot_outputs(basename, header, times, lums, alts, frag_t, frag_h, velocities, accels, qs, strength, T_stags, p_stags, T_surfs,\
                 E_rad_totals, peak_power, peak_rad_i, E_rad_tot, E_event):

    owtname = basename + f'/bolide_plots.pdf'
    pdf_pages = PdfPages(owtname)

    # Plot outputs
    # luminosity in Watts vs time in linear scale
    fig1 = plt.figure(figsize=(9.0, 6.5))
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * lums[0]
    ymax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, lums)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminosity (Watts)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.gcf().text(.14, .77, f"Frag Alt = {frag_h / 1000.:.1f} km", fontsize=10, color='black')
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Approx Event E = {E_event / jkt:.3e} kt", fontsize=10, color='black')
    plt.axhline(peak_power, linestyle='--', color='red', label=f"Peak Power = {peak_power:.3e} W")
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Light Curve: Luminosity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    pdf_pages.savefig(fig1)

    # luminosity in Watts vs time in log scale
    fig2 = plt.figure(figsize=(9.0, 6.5))
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * lums[0]
    ymax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, lums)
    plt.xscale('linear')
    plt.yscale('log')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminosity (Watts)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.gcf().text(.14, .77, f"Frag Alt = {frag_h / 1000.:.1f} km", fontsize=10, color='black')
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Approx Event E = {E_event / jkt:.3e} kt", fontsize=10, color='black')
    plt.axhline(peak_power, linestyle='--', color='red', label=f"Peak Power = {peak_power:.3e} W")
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Light Curve: Luminosity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    pdf_pages.savefig(fig2)

    # radiant intensity in W/sr vs time in linear scale
    fig3 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * radiant_intensities[0]
    ymax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, radiant_intensities)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radiant Intensity (W/sr)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.gcf().text(.14, .77, f"Frag Alt = {frag_h / 1000.:.1f} km", fontsize=10, color='black')
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Approx Event E = {E_event / jkt:.3e} kt", fontsize=10, color='black')
    plt.axhline(peak_rad_i, linestyle='--', color='red', label=f"Peak Rad I = {peak_rad_i:.3e} W/sr")
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Light Curve: Radiant Intensity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    pdf_pages.savefig(fig3)

    # radiant intensity in W/sr vs time in log scale
    fig4 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * radiant_intensities[0]
    ymax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, radiant_intensities)
    plt.xscale('linear')
    plt.yscale('log')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radiant Intensity (W/sr)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.gcf().text(.14, .77, f"Frag Alt = {frag_h / 1000.:.1f} km", fontsize=10, color='black')
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Approx Event E = {E_event / jkt:.3e} kt", fontsize=10, color='black')
    plt.axhline(peak_rad_i, linestyle='--', color='red', label=f"Peak Rad I = {peak_rad_i:.3e} W/sr")
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Light Curve: Radiant Intensity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    pdf_pages.savefig(fig4)

    # luminosity vs altitude
    fig5 = plt.figure(figsize=(9.0, 6.5))
    ymin = 0.0
    ymax = 100.
    xmin = 0.8 * lums[0]
    xmax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(lums, alts / 1000.)
    plt.xscale('linear')
    plt.xlabel('Luminosity (Watts)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.gcf().text(.73, .77, f"Frag Time = {frag_t:.3f} s", fontsize=10, color='black')
    plt.gcf().text(.69, .74, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.65, .71, f"Approx Event E = {E_event / jkt:.3e} kt", fontsize=10, color='black')
    plt.axvline(peak_power, linestyle='--', color='red', label=f'Peak Power = {peak_power:.3e} W')
    if frag_h is not None: 
        plt.axhline(frag_h / 1000., linestyle='--', color='orange', label=f'Frag Alt = {frag_h / 1000.:.1f} km')
    plt.title(f'Bolide Light Curve: Altitude vs Luminosity\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    pdf_pages.savefig(fig5)

    # radiant intensity in W/sr vs altitude
    fig6 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    ymin = 0.0
    ymax = 100.
    xmin = 0.8 * radiant_intensities[0]
    xmax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(radiant_intensities, alts / 1000.)
    plt.xscale('linear')
    plt.xlabel('Radiant Intensity (W/sr)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.gcf().text(.73, .77, f"Frag Time = {frag_t:.3f} s", fontsize=10, color='black')
    plt.gcf().text(.69, .74, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.65, .71, f"Approx Event E = {E_event / jkt:.3e} kt", fontsize=10, color='black')
    plt.axvline(peak_rad_i, linestyle='--', color='red', label=f'Peak Rad I = {peak_rad_i:.3e} W/sr')
    if frag_h is not None: 
        plt.axhline(frag_h / 1000., linestyle='--', color='orange', label=f'Frag Alt = {frag_h / 1000.:.1f} km')
    plt.title(f'Bolide Light Curve: Altitude vs Radiant Intensity\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    pdf_pages.savefig(fig6)

    # altitude vs time
    fig7 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, alts / 1000.)
    plt.xscale('log')
    plt.xlabel('Time (seconds))')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_h is not None:
        plt.axhline(frag_h / 1000., linestyle='--', color='orange', label=f"Frag Alt = {frag_h / 1000.:.1f} km")
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Altitude vs Time\n')
    plt.grid(True)
    if frag_h or frag_t is not None:
        plt.legend()
    pdf_pages.savefig(fig7)

    # velocity vs time
    fig8 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, velocities/1000.)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Velocity (km/s)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Velocity vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig8)

    # acceleration vs time
    fig9 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, accels)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Acceleration (m/s^2)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Acceleration vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig9)

    # dynamic pressure vs time
    fig10 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, qs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Dynamic Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axhline(strength, linestyle='--', color='magenta', label=f"Mat Strength = {strength:.1e} Pa")
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Dynamic Pressure vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig10)

    # dynamic pressure vs altitude
    fig11 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(alts/1000., qs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Dynamic Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axhline(strength, linestyle='--', color='magenta', label=f"Mat Strength = {strength:.1e} Pa")
        plt.axvline(frag_h/1000., linestyle='--', color='orange', label=f"Frag Alt = {frag_h / 1000.:.1f} km")
    plt.title(f'Bolide Dynamic Pressure vs Altitude\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig11)

    # stagnation temperature vs time
    fig12 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, T_stags)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Stagnation Temperature (K)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Stagnation Temperature vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig12)

    # stagnation pressure vs time
    fig13 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, p_stags)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Stagnation Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Stagnation Pressure vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig13)

    # surface temperature vs time
    fig14 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, T_surfs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Surface Temperature (K)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Surface Temperature vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig14)

    # total radiated energy vs time
    fig15 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, E_rad_totals)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Total Radiated Energy (J)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_t is not None:
        plt.axvline(frag_t, linestyle='--', color='green', label=f"Frag Time = {frag_t:.3f} s")
    plt.title(f'Bolide Total Radiated Energy vs Time\n')
    if frag_t is not None:
        plt.legend()
    plt.grid(True)
    pdf_pages.savefig(fig15)

    pdf_pages.close()

    print("All post-processing complete")


# functions for bolide luminosity and fragmentation calculations

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

# function to update mass of bolide
def mass_update_calc(mass, dmdt, dt):

    mass = mass - (dmdt * dt)
    if mass < 0.:
        mass = 0.
    print(f"updated mass       = {mass:.4e} kg")

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

# function to calculate initial kinematic quantities
def trajectory_init(bolide_outputs, xstart, altstart, velocity, tstart, theta_deg):

    # initial altitude
    z   = altstart + R_Earth
    R   = z
    alt = altstart
    print(f"initial altitude = {z:.4e} meters")
    bolide_outputs.write(f"Initial altitude = {z:.4e} meters\n")

    # initial position
    x = xstart
    print(f"initial range    = {x:.4e} meters")
    bolide_outputs.write(f"Initial position = {x:.4e} meters\n")


    # initial velocity
    v = velocity
    print(f"initial velocity = {v:.4e} m/s")
    bolide_outputs.write(f"Initial velocity = {v:.4e} m/s\n")

    # initial time
    t = tstart
    print(f"initial time = {t:.4f} seconds")
    bolide_outputs.write(f"Initial time     = {t:.4e} s\n\n")

    # printing out initial entry angle
    theta = np.radians(theta_deg)
    print(f"initial angle of entry in radians     = {theta:.3f} radians\n")
    bolide_outputs.write(f"Initial angle of entry = {theta_deg:.3f} deg\n")
    bolide_outputs.write(f"Initial angle of entry = {theta:.3f} rad\n\n")

    # calculate initial x and y velocity components of bolide
    vx = v * math.cos(theta)         # horizontal velocity
    vz = -v * math.sin(theta)        # vertical velocity (negative since downward)
    theta = np.arctan(abs(vz) / vx)  # angle of motion from horizontal (used for drag direction)
    bolide_outputs.write(f"Initial x-component velocity = {vx:.4e} m/s\n")
    bolide_outputs.write(f"Initial z-component velocity = {vz:.4e} m/s\n")

    # initialize initial alpha angle
    alpha = 0.  # radians
    bolide_outputs.write(f"Initial alpha angle          = {alpha * rad2deg:.4e} deg\n")

    # initialize initial delta angle
    delta = (pi / 2.) - theta
    bolide_outputs.write(f"Initial delta angle          = {delta * rad2deg:.4e} deg\n")

    # initialize initial gamma angle
    gamma = delta + alpha
    bolide_outputs.write(f"Initial gamma angle          = {gamma * rad2deg:.4e} deg\n")

    # initialize initial flight path angle
    phi = (pi / 2.) - (alpha + delta)
    bolide_outputs.write(f"Initial flight path angle    = {phi * rad2deg:.4e} deg\n")

    # initialize initial downrange distance
    s = C_Earth * (alpha / (2. * pi))  # meters
    bolide_outputs.write(f"Initial downrange distance   = {s:.4e} meters\n")

    # Earth surface calcs
    surfx = R_Earth * math.cos((pi / 2.) - alpha)
    surfz = R_Earth * math.sin((pi / 2.) - alpha)
    bolide_outputs.write(f"Initial Earth surface x-comp = {surfz:.4e} meters\n")
    bolide_outputs.write(f"Initial Earth surface z-comp = {surfz:.4e} meters\n\n")

    return x, z, R, alt, v, t, theta, vx, vz, alpha, delta, gamma, phi, s, surfx, surfz

# function to update bolide trajectory quantities
def trajectory_update(x, z, R, v, theta, vx, vz, alpha, delta, gamma, phi, s, Fdrag, dt, m):

    g = 9.80665  # acceleration due to gravity at sea-level (m/s^2)

    # x, z, R, alt, v, t, theta, vx, vz, alpha, delta, gamma, phi, s are given by trajectory_init

    # store given quantities for the timestep
    x0     = x
    z0     = z
    R0     = R
    vx0    = vx
    vz0    = vz
    v0     = v

    theta0 = theta
    alpha0 = alpha
    delta0 = delta
    gamma0 = gamma
    phi0   = phi

    s0     = s

    # calculate sum of forces in x and z directions
    g = g * (R_Earth / (R_Earth + z))**2.  # Earth radius = 6378.14 km

    # calc accels first
    # gravitational forces
    Fgx = -m * g * math.sin(alpha0)
    Fgz = -m * g * math.cos(alpha0)
    gx = Fgx / m
    gz = Fgz / m

    # and drag forces
    Fdragx = -Fdrag * math.cos(theta0)
    Fdragz = Fdrag * math.sin(theta0)
    dragx = Fdragx / m
    dragz = Fdragz / m
    
    # then sum up the accelerations
    ax = gx + dragx
    az = gz + dragz
    acc = (ax**2. + az**2.)**0.5

    print(f"ax  = {ax} m/s^2")
    print(f"az  = {az} m/s^2")
    print(f"acc = {acc} m/s^2")

    # next calculate final velocities
    vx = vx0 + ax * dt
    vz = vz0 + az * dt
    v = (vx**2. + vz**2.)**0.5

    print(f"vx = {vx} m/s")
    print(f"vz = {vz} m/s")
    print(f"v  = {v} m/s")

    # then calculate final positions
    x = x0 + vx0 * dt + 0.5 * ax * dt**2.
    z = z0 + vz0 * dt + 0.5 * az * dt**2.

    dx = x - x0
    dz = z - z0

    d = (dx**2. + dz**2.)**0.5

    print(f"x = {x} m")
    print(f"z = {z} m")
    print(f"d = {d} m")

    # now find new altitude using law of cosines
    R = (R0**2. + d**2. - 2. * R0 * d * math.cos(gamma0))**0.5
    alt = R - R_Earth
    print(f"alt = {alt} m")

    # find delta_alpha and new alpha angle with law of sines
    sindalpha = math.sin(gamma0) * (d / R)
    dalpha = math.asin(sindalpha)
    alpha = alpha0 + dalpha
    print(f"dalpha = {dalpha * rad2deg} deg")
    print(f"alpha  = {alpha * rad2deg} deg")

    # find new theta angle
    tantheta = abs(vz) / vx
    theta = math.atan(tantheta)
    dtheta = theta - theta0
    print(f"theta = {theta * rad2deg} deg")
    
    # find new delta angle
    delta = (pi / 2.) - theta
    print(f"delta = {delta * rad2deg} deg")

    # find new downrange distance s
    s = C_Earth * (alpha / (2. * pi))
    ds = s - s0
    print(f"s  = {s} m")
    print(f"ds = {ds} m")

    # find new flight path angle phi
    phi = (pi / 2.) - (alpha + delta)
    print(f"phi = {phi * rad2deg} deg")

    # find new gamma angle
    gamma = delta + alpha
    print(f"gamma = {gamma * rad2deg} deg")

    # Earth surface calcs
    surfx = R_Earth * math.cos((pi / 2.) - alpha)
    surfz = R_Earth * math.sin((pi / 2.) - alpha)
    R_E = (surfx**2. + surfz**2.)**0.5
    print(f"surfx = {surfx} m")
    print(f"surfz = {surfz} m")
    print(f"Radius of Earth = {R_E} m")

    return x, z, R, alt, v, theta, vx, vz, acc, alpha, delta, gamma, phi, s, surfx, surfz

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
def dynamic_pressure(rho, v):

    q = 0.5 * rho * v**2.

    return q

# calculate initial velocities
def vel_comp_calc(theta, velocity):

    # forward is positive x, downward is negative y
    vx = velocity * np.cos(theta)         # horizontal velocity
    vz = -velocity * np.sin(theta)        # vertical velocity (negative since downward)
    v = np.sqrt(vx**2 + vz**2)            # total speed
    theta = np.arctan(abs(vz) / vx)           # angle of motion from horizontal (used for drag direction)

    return vx, vz, v

def angle_calc(vx, vy):

    angle = np.arctan2(-vy, vx)           # angle of motion from horizontal (used for drag direction)

    return angle

# calculate kinetic energy
def KE_calc(mass, v):

    KE = 0.5 * mass * v**2.

    return KE

# calculate drag
def Fdrag_calc(rho, v, area, Cd):

    # calculate cross-sectional area of a sphere
    # added 27 Oct 2025 CJM
    # area_cs = pi * r^2
    # if area = 4 * pi * r^2, then
    # area_cs = area / 4.
    area_cs = area / 4.

    # drag is defined as
    Fdrag = 0.5 * Cd * rho * v**2. * area_cs
    # units are kg/m^3 * m^2/s^2 * m^2 = kg * m/s^2 = N

    return Fdrag

# calculate acceleration
def accel_calc(Fdrag, mass, y, angle):

    # first, calculate acceleration due to gravity
    g = g * (R_Earth / (R_Earth + y))**2.  # Earth radius = 6378.14 km
    acc_grav = -g * np.sin(angle)

    # acc = Fdrag / m --> 0.5 * Cd * rho * v^2 * area / mass
    acc_drag = Fdrag / mass

    # total acceleration
    acc = acc_drag + acc_grav

    return acc

# calculate heat flux
def heatflux_calc(Fdrag, v, area):

    # cross-sectional area of a sphere
    # added 27 Oct 2025 CJM
    area_cs = area / 4.

    q_h = (Fdrag * v) / area_cs

    return q_h

# calculate luminous efficiency scaling factor eta
def luminous_efficiency_scaled(v, rho, eta_0, v_ref, rho_air0, n, m):

    # m is 0.5 and n is 2 (Ceplecha et al 1998)
    eta = eta_0 * (v / v_ref)**n * (rho / rho_air0)**m

    return eta

# calculate luminosity due to drag
def lum_drag_calc(eta, Fdrag, v, dt):
        
    # acc = Fdrag / mass
    # Work = Force * Distance
    # v * dt = distance, so
    # dE = Force * Distance = Fdrag * (v * dt), units of energy (J)
    # dE / dt --> units of energy / time (J/s = W)
    dE = Fdrag * (v * dt)
    luminosity_drag = eta * dE / dt

    return luminosity_drag

# dynamic emissivity function
def epsilon_dynamic_calc(rho, p_stag, p_air0):

    # Scaled between 0.01 and 0.1, based on atmospheric density and stagnation pressure
    # can refine with look-up tables or empirical data
    epsilon_dynamic = 0.01 + 0.02 * (p_stag / p_air0)

    # clip from 0.01 at a minimum to 0.1 at a maximum
    epsilon_dynamic = max(epsilon_dynamic, 0.01)
    epsilon_dynamic = min(epsilon_dynamic, 0.1)

    return epsilon_dynamic

# calculate effective area
def A_eff_calc(k_shock, area, fragmented, n_fragments):

    A_eff = k_shock * area

    # do not artificially increase the effective area because now we are tracking multiple
    # fragments CJM 25 Oct 2025
    # if fragmented:
       # A_eff = (1.5 + 0.5 * n_fragments) * area

    return A_eff

# calculate luminosity due to shock front
def lum_shock_calc(A_eff, T_stag, epsilon_dynamic, M, rho, rho_air0, SISF_flag):

    # basic luminosity equation
    # luminosity_shock = epsilon_dynamic * sigma_sb * A_eff * T_stag**4.

    # Shock Intensity Scaling calculations, CJM 17 Oct 2025
    if(SISF_flag == 1):

        # Smooth Mach-dependent scaling 
        # see Revelle (1997) and Ceplecha et al (1998)
        # luminous efficiency scales with Mach number, and shock-produced light only appears when
        # the flow is strongly hypersonic
        mach_scale = ((M - 3.) / 5.0)**0.5  # grows from 0 to 1 as M increases from 3 to 8, then increases

        if M < 3:
            mach_scale = 0.0 # no significant shock radiation below Mach 3

        # Atmospheric density scaling
        # motivation: accounts for increased stagnation heating and radiative efficiency in denser air,
        # but suppresses unrealistic high-altitude luminosity
        rho_scale = (rho / rho_air0)**0.8

        # Shock Intensity Scaling Factor is mach_scale * rho_scale
        SISF = mach_scale * rho_scale

    elif(SISF_flag == 0):
        
        # No shock intensity scaling, CJM 17 Oct 2025
        SISF = 1.

    elif(SISF_flag != 1 and SISF_flag != 0):

        # print error message and end simulation
        print(f"Shock Intensity Scaling Factor Flag incorrectly set: ending simulation")

        return

    # Now calculate luminosity of the shock
    luminosity_shock = SISF * epsilon_dynamic * sigma_sb * A_eff * T_stag**4.

    return luminosity_shock, SISF

# calculate change in kinetic energy
def E_dot_mech_calc(Fdrag, v):

    # calculating change in mechanical energy per timestep
    E_dot_mech = Fdrag * v
    # units are kg * m/s^2 * m/s = kg * m^2/s^3 = J/s = W

    return E_dot_mech

# calculate stagnation temperature
def T_stag_calc(M, T):

    T_stag = T * (1. + 0.5 * (gamma - 1.) * M**2.)

    return T_stag

# calculate stagnation pressure
def p_stag_calc(q, p):

    p_stag = q + p

    return p_stag

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

# estimate compressive strength based on density and porosity
def compute_strength(init_strength):

    if init_strength > 0:
        # use input material strength in Pa
        strength = init_strength
    else:
        # compute mat strength normalized to average mat strength
        # of ordinary chondrite and primitive chondrite (Brown et al 2002)
        strength = 1.e5 * (1. - porosity) * (density / 3000.)  # density/3000 term normalizes to silicate

    print(f"strength of bolide is {strength} Pa")

    return strength

def frag_init_calc(bolide_outputs, alt, t, fragmented, frag_altitude, frag_time, n_fragments, density, mass):
        
        # change fragmentation logic to True
        fragmented      = True

        # record fragmentation altitude and time
        frag_altitude   = alt
        frag_time       = t

        # divide mass up equally amongst fragments and recompute total surface area
        mass_fragment   = mass / n_fragments
        r_fragment      = ((3. * mass_fragment) / (4. * pi * density))**(1. / 3.)
        dia_fragment    = 2. * r_fragment
        area_fragment   = 4. * pi * (dia_fragment / 2.)**2.
        area            = n_fragments * area_fragment
        mass            = mass_fragment * n_fragments

        # write outputs
        print(f"Bolide fragmented at {frag_altitude} meters and {frag_time} seconds!\n")
        bolide_outputs.write(f"Bolide fragmented at {frag_altitude:.4e} m and {frag_time:.4e} s!\n\n")

        return fragmented, frag_altitude, frag_time, area, mass


# find intermediate total radiated energy
def E_rad_total_calc(E_rad_total, luminosity, dt):

    E_rad_total = E_rad_total + luminosity * dt
    print(f"total radiated energy is {E_rad_total:.3e} Joules")
    print(f"total radiated energy is {E_rad_total / jkt:.3e} kilotons")

    return E_rad_total

# find peak luminosity and peak radiant intensity
def peak_find(bolide_outputs, lum_array):

    # find peak luminosity in W
    peak_power = np.max(lum_array)
    print(f"bolide peak luminosity            = {peak_power:.3e} W")
    bolide_outputs.write(f"bolide peak luminosity            = {peak_power:.3e} W\n")

    # find peak radiant intensity in W/sr
    peak_rad_i = peak_power / (4. * math.pi)
    print(f"bolide peak radiant intensity     = {peak_rad_i:.3e} W/sr")
    bolide_outputs.write(f"bolide peak radiant intensity     = {peak_rad_i:.3e} W\n")

    return peak_power, peak_rad_i

# calculate total radiated energy
def final_rad_E_calc(bolide_outputs, E_rad_ttl_array):

    E_rad_tot = E_rad_ttl_array[-1, -1]

    print(f"total radiated energy of bolide   = {E_rad_tot:.3e} J")
    bolide_outputs.write(f"total radiated energy of bolide   = {E_rad_tot:.3e} J\n")

    print(f"total radiated energy of bolide   = {E_rad_tot / jkt:.3e} kilotons")
    bolide_outputs.write(f"total radiated energy of bolide   = {E_rad_tot / jkt:.3e} kilotons\n")

    return E_rad_tot

def E_event_calc(bolide_outputs, E_rad_tot):

    # from Brown et al, 2002, the totl event energy can be calculated by
    # E_event = 8.2508 * E_optical^0.885
    # this assumes a 6,000 K blackbody as the radiater (bolide)
    # for simplicity, since we only calculate the total radiated energy,
    # not the optical radiated energy, we will divide E_rad_tot by 2
    # this can be updated in the future

    # first, must convert total radiated energy into kt from J
    E_rad_tot = E_rad_tot / jkt

    # then calculate event energy using the 1/2 relationship stated above
    E_event = 8.2508 * (E_rad_tot / 2.)**0.885

    # then convert back to joules
    E_event = E_event * jkt

    print(f"approximate total energy of event = {E_event:.3e} J")
    bolide_outputs.write(f"approximate total energy of event = {E_event:.3e} J\n")
    
    print(f"approximate total energy of event = {E_event / jkt:.3e} kilotons")
    bolide_outputs.write(f"approximate total energy of event = {E_event / jkt:.3e} kilotons\n")

    return E_event
    
# main loop function
def bolide_luminosity_model(diameter, velocity, theta_deg, zstart, xstart, density, porosity, tstart, dt, tstop, n_fragments, SISF_flag):

    # initializing counter
    timestep = int(1)
    print(f"first timestep = {timestep}")

    # compute initial area
    area = area_init_calc(diameter)

    # compute initial mass
    mass = mass_init_calc(diameter, density)

    # compute material strength of bolide
    strength = compute_strength(init_strength)
    
    # define header for plotting
    header = (f"Dia: {diameter:.1f} m  Vel: {velocity/1000.:.1f} km/s  theta: {theta_deg:.1f} deg from horiz  Mat Strength: {strength:.1e} Pa")

    # open output file and write first line and header
    bolide_outputs = open('bolide_out.txt', 'w')
    bolide_outputs.write(f"Bolide Luminosity and Fragmentation Simulation\n")
    bolide_outputs.write(f"{header}\n\n")

    # write out initial material properties of bolide
    bolide_outputs.write(f"Initial area of bolide      = {area:.4e} m^2\n")
    bolide_outputs.write(f"Initial mass of bolide      = {mass:.4e} kg\n")
    bolide_outputs.write(f"Material strength of bolide = {strength:.4e} Pa\n\n")

    # set up initial quantities for trjectory
    x, z, R, alt, v, t, theta, vx, vz, alpha, delta, gamma, phi, s,\
    surfx, surfz = trajectory_init(bolide_outputs, xstart, zstart, velocity, tstart, theta_deg)

    # initialize total energy radiated
    E_rad_total = 0.0
    print(f"setting total radiated energy to {E_rad_total:.4e} J\n")
    bolide_outputs.write(f"Setting total radiated energy to {E_rad_total:.4e} J\n\n")

    # initializing arrays
    # atmospheric temperature, pressure, density, and sound speed
    T_array         = np.zeros(1)
    p_array         = np.zeros(1)
    rho_array       = np.zeros(1)
    a_array         = np.zeros(1)

    # luminosity, altitudes, times
    lum_array       = np.zeros(1)
    SISF_array      = np.zeros(1)
    alt_array       = np.zeros(1)
    time_array      = np.zeros(1)
    mass_array      = np.zeros(1)
    dia_array       = np.zeros(1)
    area_array      = np.zeros(1)

    # energies, velocities, and accelerations
    KE_array        = np.zeros(1)
    v_array         = np.zeros(1)
    accel_array     = np.zeros(1)
    heatflux_array  = np.zeros(1)
    Fdrag_array     = np.zeros(1)
    E_rad_ttl_array = np.zeros(1)

    # shock physics/flow quantities
    q_array         = np.zeros(1)
    Mach_array      = np.zeros(1)
    T_stag_array    = np.zeros(1)
    p_stag_array    = np.zeros(1)
    T_surf_array    = np.zeros(1)

    # initial conditions for bolide fragmentation
    fragmented      = False
    frag_altitude   = None
    frag_time       = None

    # main loop
    while timestep < mxcycl:   

        # 1. calculate atmospheric quantities and thermal/flow variables:
        # compute temperature at altitude
        T = temp_altitude_calc(alt)

        # compute atmospheric pressure
        p = atmos_pressure_calc(alt, T)

        # compute air density
        rho = atmos_density_calc(alt, p, T)

        # compute speed of sound
        a = sound_speed_calc(T, p, rho)

        # calculate Mach number of bolide
        M = mach_num_calc(v, a)

        # calculate dynamic pressure at bolide
        q = dynamic_pressure(rho, v)

        # compute drag of bolide
        Fdrag = Fdrag_calc(rho, v, area, Cd)

        # calculate kinetic energy of bolide
        KE = KE_calc(mass, v)

        # compute stagnation temperature
        T_stag = T_stag_calc(M, T)

        # compute stagnation pressure
        p_stag = p_stag_calc(q, p)

        # compute heat flux at bolide
        q_h = heatflux_calc(Fdrag, v, area)

        # compute surface temperature
        T_surf = T_surf_calc(epsilon, q_h)

        # 2. perform fragmentation if fragmentation conditions are met
        if not fragmented and q > strength:
            fragmented, frag_altitude, frag_time, area, mass = frag_init_calc(bolide_outputs, alt, t, fragmented,\
                                                                              frag_altitude, frag_time, n_fragments,\
                                                                              density, mass)

        # 3. calculate luminous output quantities
        # calculate velocity- and density-dependent luminous efficiency
        eta = luminous_efficiency_scaled(v, rho, eta_0, v_ref, rho_air0, n, m)

        # calculate change in KE
        E_dot_mech = E_dot_mech_calc(Fdrag, v)

        # compute change in energy and luminosity due to drag
        # uses scaled luminous efficiency eta calculated in previous step (Ceplecha et al 1998)
        luminosity_drag = lum_drag_calc(eta, Fdrag, v, dt)

        # calculate dynamic emissivity for shock luminosity calculation
        epsilon_dynamic = epsilon_dynamic_calc(rho, p_stag, p_air0)

        # calculate effective area for shock luminosity calculation
        A_eff = A_eff_calc(k_shock, area, fragmented, n_fragments)

        # compute luminosity due to luminous shock front
        luminosity_shock, SISF = lum_shock_calc(A_eff, T_stag, epsilon_dynamic, M, rho, rho_air0, SISF_flag)

        # compute total luminosity
        luminosity = luminosity_drag + luminosity_shock

        # --- smoothing step: only apply at the timestep immediately after fragmentation ---
        # if just_fragmented:
            # luminosity = alpha * luminosity_prev + (1 - alpha) * luminosity
            # just_fragmented = False  # Reset flag after smoothing

        # update previous luminosity for next timestep
        luminosity_prev = luminosity

        # 4. update dynamical variables: acceleration, velocity, positions, and angle of entry
        x, z, R, alt, v, theta, vx, vz, acc, alpha, delta, gamma, phi, s,\
        surfx, surfz = trajectory_update(x, z, R, v, theta, vx, vz, alpha, delta, gamma, phi, s, Fdrag, dt, mass)

        # 5. update physical properties: mass, diameter, area
        # mass
        dmdt = dmdt_calc(q_h, area, L_ablation)  # compute mass loss rate due to ablation
        mass  = mass_update_calc(mass, dmdt, dt)
        
        # diameter
        diameter = diameter_update_calc(mass, density, n_fragments, fragmented)

        # area
        area = area_update_calc(diameter, n_fragments, fragmented)

        # 6. update energy and time quantities
        # calculate total radiated energy
        E_rad_total = E_rad_total_calc(E_rad_total, luminosity, dt)

        # 7. assign quantities to array indices
        # temperature, pressure, density, and sound speed
        T_array[timestep-1]         = T
        p_array[timestep-1]         = p
        rho_array[timestep-1]       = rho
        a_array[timestep-1]         = a

        # luminosity, altitudes, and times
        lum_array[timestep-1]       = luminosity 
        SISF_array[timestep-1]      = SISF
        alt_array[timestep-1]       = alt 
        time_array[timestep-1]      = t 
        mass_array[timestep-1]      = mass
        dia_array[timestep-1]       = diameter
        area_array[timestep-1]      = area

        # energies and velocities
        KE_array[timestep-1]        = KE
        v_array[timestep-1]         = v
        accel_array[timestep-1]     = acc
        heatflux_array[timestep-1]  = q_h
        Fdrag_array[timestep-1]     = Fdrag
        E_rad_ttl_array[timestep-1] = E_rad_total

        # shock physics/flow quantities
        q_array[timestep-1]         = q
        Mach_array[timestep-1]      = M
        T_stag_array[timestep-1]    = T_stag
        p_stag_array[timestep-1]    = p_stag
        T_surf_array[timestep-1]    = T_surf

        # finally, update timestep
        t = t + dt
        print(f"updated time       = {t:.4f} seconds\n")
        
        # write to output file
        wr_out(bolide_outputs, luminosity, z, t, mass, diameter, area, KE, v, acc, q_h, Fdrag,\
               E_rad_total, q, M, T_stag, p_stag, T_surf, T, p, rho, a, SISF, theta_deg)

        # 8. simulation end conditions
        # stop simulation if time exceeds tstop
        if t > tstop:
            print(f"time exceeded tstop: end simulation\n")
            break

        # stop simulation if mass drops to zero
        if mass <= 0.:
            print(f"mass is <= 0 kg: end simulation\n")
            break

        # stop simulation if altitude is <= zero
        if alt <= 0.:
            print(f"altitude is <= 0 m: end simulation\n")
            break
        
        # stop simulation if velocity is <= zero
        if v <= 0.:
            print(f"velocity is <= 0 m/s: end simulation")
            break

        # 9. if simulation not ended, add another row to arrays, then update counter
        # temperature, pressure, density, and sound speed
        T_array         = np.vstack([T_array, new_row])
        p_array         = np.vstack([p_array, new_row])
        rho_array       = np.vstack([rho_array, new_row])
        a_array         = np.vstack([a_array, new_row])

        # luminosities, altitudes, and times
        lum_array       = np.vstack([lum_array, new_row])
        SISF_array      = np.vstack([SISF_array, new_row])
        alt_array       = np.vstack([alt_array, new_row])
        time_array      = np.vstack([time_array, new_row]) 
        mass_array      = np.vstack([mass_array, new_row])
        dia_array       = np.vstack([dia_array, new_row])
        area_array      = np.vstack([area_array, new_row])

        # energies and velocities
        KE_array        = np.vstack([KE_array, new_row])
        v_array         = np.vstack([v_array, new_row])
        accel_array     = np.vstack([accel_array, new_row])
        heatflux_array  = np.vstack([heatflux_array, new_row])
        Fdrag_array     = np.vstack([Fdrag_array, new_row])
        E_rad_ttl_array = np.vstack([E_rad_ttl_array, new_row])

        # shock physics/flow quantities
        q_array         = np.vstack([q_array, new_row])
        Mach_array      = np.vstack([Mach_array, new_row])
        T_stag_array    = np.vstack([T_stag_array, new_row])
        p_stag_array    = np.vstack([p_stag_array, new_row])
        T_surf_array    = np.vstack([T_surf_array, new_row])

        # updating to next timestep counter
        timestep = timestep + 1
        print(f"next timestep = {timestep}\n")

    # 10. calculate final outputs once simulation ends and close output file
    # find peak luminosity and total radiated energy
    peak_power, peak_rad_i = peak_find(bolide_outputs, lum_array)

    # find final radiated energy
    E_rad_tot = final_rad_E_calc(bolide_outputs, E_rad_ttl_array)

    # find total energy deposited into environment
    E_event = E_event_calc(bolide_outputs, E_rad_tot)

    # print out fragmentation data if fragmented
    if fragmented:
        print(f"fragmentation time                = {frag_time:.3f} seconds")
        print(f"fragmentation altitude            = {frag_altitude / 1000.:.3f} km")

    # close output file
    bolide_outputs.close()

    return (
        header,
        strength,
        T_array,
        p_array,
        rho_array,
        a_array,
        lum_array,
        SISF_array,
        alt_array,
        time_array,
        mass_array,
        dia_array,
        area_array,
        frag_altitude,
        frag_time,
        KE_array,
        v_array,
        accel_array,
        heatflux_array,
        Fdrag_array,
        E_rad_ttl_array,
        q_array,
        Mach_array,
        T_stag_array,
        p_stag_array,
        T_surf_array,
        peak_power,
        peak_rad_i,
        E_rad_tot,
        E_event
    )

# Main Program
basename  = f"/Users/julie/Desktop/Projects/bolides/"
inputpath = f"{basename}bolum_in.txt"

# read in and assign input quantities
diameter, velocity, theta_deg, density, porosity, init_strength,\
zstart, xstart, tstart, dt, tstop,\
n_fragments, SISF_flag = rdinput(inputpath)

# write out initial parameters
write_init_params(diameter, velocity, theta_deg, density, porosity, init_strength,\
                  zstart, xstart, tstart, dt, tstop, n_fragments)

# call main loop
header, strength, T_array, p_array, rho_array, a_array, lum_array, SISFs_array,\
alt_array, time_array, mass_array, dia_array, area_array, frag_altitude, frag_time,\
KE_array, v_array, accel_array, heatflux_array, Fdrag_array,\
E_rad_ttl_array, q_array, Mach_array, T_stag_array, p_stag_array, T_surf_array,\
peak_power, peak_rad_i, E_rad_tot, E_event = bolide_luminosity_model(diameter, velocity, theta_deg,\
                                                                     zstart, xstart, density, porosity,\
                                                                     tstart, dt, tstop, n_fragments,\
                                                                     SISF_flag)

# error checking
print(f"total cycles: {lum_array.size}")
print(f"total cycles: {time_array.size}")
print(f"lum array shape:  {lum_array.shape}")
print(f"time array shape: {time_array.shape}")

print(f"STOP all done bolide simulation complete")

# plot outputs
plot_outputs(basename, header, time_array, lum_array, alt_array, frag_time, frag_altitude,\
             v_array, accel_array, q_array, strength, T_stag_array, p_stag_array,\
             T_surf_array, E_rad_ttl_array, peak_power, peak_rad_i, E_rad_tot, E_event)