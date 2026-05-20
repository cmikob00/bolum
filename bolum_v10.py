#######################################################################
# Bolide Luminosity Model (bolum) Version 10                           #
# Author: C.J. Miko                                                   #
# Instructions: Modify inputs in the input deck (bolum_in.) and run   #
# Please see available readme for additional information              #
# Total line count: 1358 lines                                        #
# This code is UNCLASSIFIED                                           #
#######################################################################

import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#### File Path ####
basename = f"/Users/julie/Desktop/Projects/bolides/Events/"
####################

############ Variable Declarations ####################

# Constants
pi       = math.pi
g        = 9.80665            # gravitational acceleration for Earth
R        = 287.05             # dry air specific gas constant (J / mol*kg)
sigma_sb = 5.670374419e-8     # Stefan-Boltzmann constant
rad2deg  = 180. / math.pi     # radians to degrees
deg2rad  = math.pi / 180      # degrees to radians
jkt      = 4.185e12           # conversion factor joules/kt
R_Earth  = 6378.14e3          # Earth radius
C_Earth  = 2. * pi * R_Earth  # circumference of Earth

# air constants
gamma      = 1.2      # adiabatic constant of air (1.2 = shocked air)
rho_air0   = 1.225    # sea-level air density
p_air0     = 101325.  # sea-level air pressure
H          = 7160.0   # scale height of atmosphere
T0         = 288.15   # ambient air temperature at sea-level
L          = 0.0065   # temperature lapse rate of atmosphere

# bolide constants
epsilon           = 0.9    # emissivity of bolide
flare_duration    = 0.1   # flare duration after fragmentation in sec

# other simulation parameters
mxcycl  = 100000           # maximum number of cycles
new_row = np.zeros(1)

#####################################################

############# I/O functions #######################
# reading in input deck
def rdinput(inputpath):

    params = {}

    with open(inputpath, "r") as f:
        lines = f.readlines()[2:]

    for line in lines:

        # skip blank/comment lines
        if "=" not in line:
            continue

        key, rest = line.split("=", 1)

        key = key.strip()
        value_str = rest.split("#")[0].strip()

        try:
            value = float(value_str)
        except ValueError:
            continue

        params[key] = value

    #
    # REQUIRED PARAMETERS
    #
    try:
        diameter      = params["diameter"]
        velocity      = params["velocity"]
        theta_deg     = params["theta_deg"]
        init_strength = params["init_strength"]
        density       = params["density"]

    except KeyError as e:
        raise ValueError(f"Missing required input parameter: {e}")

    #
    # OPTIONAL PARAMETERS WITH DEFAULTS
    #
    strength_secondary  = params.get("strength_secondary", 1.e15)
    strength_tertiary   = params.get("strength_tertiary", 1.e15)
    strength_quaternary = params.get("strength_quaternary", 1.e15)
    porosity            = params.get("porosity", 0.1)
    Cd                  = params.get("Cd", 1.0)
    L_ablation          = params.get("L_ablation", 5.e6)
    n_fragments         = params.get("n_fragments", 50)
    flare_duration      = params.get("flare_duration", 0.1)
    rho_tau_scale       = params.get("rho_tau_scale", 0.0)
    zstart              = params.get("zstart", 100000.0)
    xstart              = params.get("xstart", 0.0)
    tstart              = params.get("tstart", 0.0)
    dt                  = params.get("dt", 0.003)
    tstop               = params.get("tstop", 1.0)

    n_frag_init = n_fragments

    return (
        diameter,
        velocity,
        theta_deg,
        init_strength,
        strength_secondary,
        strength_tertiary,
        strength_quaternary,
        density,
        porosity,
        Cd,
        L_ablation,
        n_fragments,
        n_frag_init,
        flare_duration,
        rho_tau_scale,
        zstart,
        xstart,
        tstart,
        dt,
        tstop
    )

# writing out initial parameters to output file
def write_init_params(outputpath, diameter, velocity, theta_deg,
                      init_strength, strength_secondary, strength_tertiary, strength_quaternary,
                      density, porosity, Cd, L_ablation,
                      n_fragments, n_frag_init, flare_duration, rho_tau_scale,
                      zstart, xstart, tstart, dt, tstop):

    init_params = open(f'{outputpath}/init_params.dat', 'w')

    init_params.write(f"Initialization Parameters for Bolide Luminosity and Fragmentation Simulation\n")
    init_params.write(f"Diameter (m):                      {diameter:.3f}\n")
    init_params.write(f"velocity (m/s):                    {velocity:.3f}\n")
    init_params.write(f"Entry Angle to Horizon (deg):      {theta_deg:.3f}\n")
    init_params.write(f"Initial Strength (Pa):             {init_strength:.3e}\n")
    init_params.write(f"Secondary Strength (Pa):           {strength_secondary:.3e}\n")
    init_params.write(f"Tertiary Strength (Pa):            {strength_tertiary:.3e}\n")
    init_params.write(f"Quaternary Strength (Pa):          {strength_quaternary:.3e}\n")
    init_params.write(f"Density (kg/m^3):                  {density:.3f}\n")
    init_params.write(f"Porosity:                          {porosity:.3f}\n")
    init_params.write(f"Drag Coefficient:                  {Cd:.3f}\n")
    init_params.write(f"Heat of Ablation (J/kg):           {L_ablation:.3e}\n")
    init_params.write(f"Number of Fragments per stage:     {n_fragments}\n")
    init_params.write(f"Initially set number of fragments: {n_frag_init}\n")
    init_params.write(f"Flare Duration (sec):              {flare_duration:.3f}\n")
    init_params.write(f"Luminous Efficiency scaling term:  {rho_tau_scale:.1f}\n")
    init_params.write(f"Initial Altitude (m):              {zstart:.3f}\n")
    init_params.write(f"Initial Position (m):              {xstart:.3f}\n")
    init_params.write(f"Starting Time (s):                 {tstart:.3f}\n")
    init_params.write(f"Timestep (s):                      {dt:.3f}\n")
    init_params.write(f"Maximum Time (s):                  {tstop:.3f}\n")
    
    init_params.close()

# define general output file
def wr_out(bolide_outputs, luminosity, y, t, mass, diameter, area, KE, v, acc, q_h, Fdrag,
           E_rad_total, E_deposited, q, M, T_stag, p_stag, T_surf, T, p, rho, a, theta_deg, tau,
           frag_stage, lum_per_fragment):

    bolide_outputs.write(f"time = {t:.4e}       lum = {luminosity:.4e}     tau = {tau:.4e}   dia = {diameter:.4e}    alt = {y:.4e}\n")
    bolide_outputs.write(f"E_rad_ttl = {E_rad_total:.4e}  E_dep = {E_deposited:.4e}   KE = {KE:.4e}    Fdrag = {Fdrag:.4e}  frag_stage = {frag_stage}\n")
    bolide_outputs.write(f"v = {v:.4e}          acc = {acc:.4e}     theta_deg = {theta_deg:.3f}  mass = {mass:.4e}   area = {area:.4e}\n")
    bolide_outputs.write(f"T_stag = {T_stag:.4e}     T_surf = {T_surf:.4e}  q_h = {q_h:.4e}   q = {q:.4e}      p_stag = {p_stag:.4e}\n")
    bolide_outputs.write(f"T = {T:.4e}          p = {p:.4e}       rho = {rho:.4e}   a = {a:.4e}      M = {M:.4e}\n\n")
    if lum_per_fragment is not None:
        bolide_outputs.write(f"lum_per_fragment = {lum_per_fragment:.4e}\n")
        bolide_outputs.write(f"\n")

# define power vs time output file
def wr_power_vs_time(bolide_pvt, luminosity, t):

    bolide_pvt.write(f'{luminosity:.4e}   {t:.4e}\n')

def wr_power_vs_alt(bolide_pva, luminosity, alt):

    bolide_pva.write(f'{luminosity:.4e}   {alt:.4e}\n')

# plotting function — simplified to single luminosity
def plot_outputs(basename, header, times, lums, alts, frag_times, frag_alts, velocities, accels, qs,
                 T_stags, p_stags, T_surfs, E_rad_totals, peak_power, peak_rad_i, E_rad_tot, E_event):

    owtname = basename + f'/bolide_plots.pdf'
    pdf_pages = PdfPages(owtname)

    # 1. luminosity vs time (linear/linear)
    fig1 = plt.figure(figsize=(9.0, 6.5))
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * lums[0]
    ymax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, lums, color='blue', label='Total Luminosity')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminosity (Watts)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_power, linestyle='--', color='red', label=f"Peak = {peak_power:.3e} W")
    plt.title(f'Bolide Light Curve: Luminosity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig1)

    # 2. luminosity in Watts vs time in log/log scale
    fig2 = plt.figure(figsize=(9.0, 6.5))
    xmin = 0.003
    xmax = 1.05 * times[-1]
    ymin = 0.8 * lums[0]
    ymax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, lums, color='blue', label='Total Luminosity')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Luminosity (Watts)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
         plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_power, linestyle='--', color='red', label=f"Peak Power = {peak_power:.3e} W")
    plt.title(f'Bolide Light Curve: Luminosity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig2)

    # 3. radiant intensity in W/sr vs time in linear/linear scale
    fig3 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    xmin = -0.1
    xmax = 1.05 * times[-1]
    ymin = 0.8 * radiant_intensities[0]
    ymax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, radiant_intensities, color='blue', label='Total Rad I')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radiant Intensity (W/sr)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
         plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_rad_i, linestyle='--', color='red', label=f"Peak Rad I = {peak_rad_i:.3e} W/sr")
    plt.title(f'Bolide Light Curve: Radiant Intensity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig3)

    # 4. radiant intensity in W/sr vs time in log/log scale
    fig4 = plt.figure(figsize=(9.0, 6.5))
    # convert luminosity in watts to W/sr
    radiant_intensities = lums / (4. * pi)
    xmin = 0.003
    xmax = 1.05 * times[-1]
    ymin = 0.8 * radiant_intensities[0]
    ymax = 2. * peak_rad_i
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(times, radiant_intensities, color='blue', label='Total Rad I')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radiant Intensity (W/sr)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.14, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.14, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.14, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axhline(peak_rad_i, linestyle='--', color='red', label=f"Peak Rad I = {peak_rad_i:.3e} W/sr")
    plt.title(f'Bolide Light Curve: Radiant Intensity vs Time\n')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig4)

    # luminosity vs altitude
    fig5 = plt.figure(figsize=(9.0, 6.5))
    ymin = 0.0
    ymax = 100.
    xmin = 0.8 * lums[0]
    xmax = 2. * peak_power
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.plot(lums, alts / 1000., color='blue', label='Luminosity per Altitude')
    plt.xscale('linear')
    plt.xlabel('Luminosity (Watts)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0] / 1000., linestyle='--', color='green', alpha=0.6, label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.65, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.65, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axvline(peak_power, linestyle='--', color='red', label=f'Peak Power = {peak_power:.3e} W')
    plt.title(f'Bolide Light Curve: Altitude vs Luminosity\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
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
    plt.plot(radiant_intensities, alts / 1000., color='blue', label='Total Rad I')
    plt.xscale('linear')
    plt.xlabel('Radiant Intensity (W/sr)')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0] / 1000., linestyle='--', color='green', label=f"Frag at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec")
    plt.gcf().text(.65, .74, f"Total Rad E = {E_rad_tot / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .71, f"Total Rad E = {E_rad_tot:.3e} J", fontsize=10, color='black')
    plt.gcf().text(.65, .68, f"Approx Event E = {E_event / jkt:.3f} kt", fontsize=10, color='black')
    plt.gcf().text(.65, .65, f"Approx Event E = {E_event:.3e} J", fontsize=10, color='black')
    plt.axvline(peak_rad_i, linestyle='--', color='red', label=f'Peak Rad I = {peak_rad_i:.3e} W/sr')
    plt.title(f'Bolide Light Curve: Altitude vs Radiant Intensity\n')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig6)

    # altitude vs time
    fig7 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, alts / 1000.)
    plt.xscale('log')
    plt.xlabel('Time (seconds))')
    plt.ylabel('Altitude (km)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Altitude vs Time\n')
    plt.grid(True)
    # plt.legend(loc='upper right')
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig7)

    # velocity vs time
    fig8 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, velocities/1000.)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Velocity (km/s)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    # plt.legend()
    plt.title(f'Bolide Velocity vs Time\n')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig8)

    # acceleration vs time
    fig9 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, accels)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Acceleration (m/s^2)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Acceleration vs Time\n')
    plt.grid(True)
    # plt.legend()
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig9)

    # dynamic pressure vs time
    fig10 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, qs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Dynamic Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
        plt.axhline(init_strength, linestyle='--', color='orange', alpha=0.6)
        plt.gcf().text(.14, .84, f"Init Strength =  {init_strength:.1e} Pa", fontsize=9, color='orange')
    plt.title(f'Bolide Dynamic Pressure vs Time\n')
    # plt.legend()
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig10)

    # dynamic pressure vs altitude
    fig11 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(alts/1000., qs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Dynamic Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_alts:
        plt.axhline(frag_alts[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Dynamic Pressure vs Altitude\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig11)

    # stagnation temperature vs time
    fig12 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, T_stags)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Stagnation Temperature (K)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Stagnation Temperature vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig12)

    # stagnation pressure vs time
    fig13 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, p_stags)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Stagnation Pressure (Pa)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Stagnation Pressure vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig13)

    # surface temperature vs time
    fig14 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, T_surfs)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Surface Temperature (K)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Surface Temperature vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig14)

    # total radiated energy vs time
    fig15 = plt.figure(figsize=(9.0, 6.5))
    plt.plot(times, E_rad_totals)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Total Radiated Energy (J)')
    plt.gcf().text(.20, .89, header, fontsize=10, color='black')
    if frag_times:
        plt.axvline(frag_times[0], linestyle='--', color='green', alpha=0.6)
        plt.gcf().text(.14, .86, f"Frag at {frag_alts[0]/1000.:.1f} km", fontsize=9, color='green')
    plt.title(f'Bolide Total Radiated Energy vs Time\n')
    plt.grid(True)
    plt.gcf().text(.10, .95, 'UNCLASSIFIED', fontsize=10, color='green')
    plt.gcf().text(.80, .05, 'UNCLASSIFIED', fontsize=10, color='green')
    pdf_pages.savefig(fig15)

    pdf_pages.close()

    print("All post-processing complete")

############ New luminosity functions (v7+) ############
def tau_calc(v, mass, rho, rho_tau_scale, n_fragments):

    # This luminous efficiency function is designed to output
    # bolide luminous efficiency based on optical (broadband Si)
    # sensor response - not the entire spectrum

    v_kms = v / 1000.0
    if n_fragments > 1:
        m_eff = mass / n_fragments
    else:
        m_eff = mass

    ln_v = np.log(v_kms) if v_kms > 1e-5 else -10.0
    ln_m = np.log(m_eff) if m_eff > 1e-10 else -20.0

    if v_kms < 25.372:
        ln_tau = (0.567 - 10.307 * ln_v + 9.781 * ln_v**2. - 3.0414 * ln_v**3. + 0.3213 * ln_v**4. + 0.347 * np.tanh(0.38 * ln_m))
    else:
        ln_tau = -1.4286 + ln_v + 0.347 * np.tanh(0.38 * ln_m)

    # exponentiate
    tau = np.exp(ln_tau)

    # apply air density scaling for large bolides (small bolides don't need it)
    # this is not bulletproof
    # updated to 1.e-3 and 0.4 for v9 on 7 April 2026 cjm
    if rho_tau_scale == 1:
        tau = tau * (rho / 1.e-3)**0.45  # empirical term set to match Chelyabinsk luminosity

    # tau is now in percent, so divide by 100 to get true value
    tau = tau / 100.

    # Extra Stuff:
    # Bolometric τ (Borovička base) + airburst scaling
    # Small bolides: pure Borovička (realistic 2-7%)
    # Large (>~few m, airburst regime): empirical boost to match Chelyabinsk (12-17%)
    # if mass > 5.e5:  # Threshold: ~5m stony (~1e6 kg); below = fireballs
    #     tau_scale = 1.0 + 2.8 * (mass / 1.e7)**0.35  # ~1x small, 3.5x Chelyabinsk
    #     tau = tau * tau_scale
    # else: tau remains ~6% for v~19 km/s, m~1e3 kg (perfect for Benešov/Košice)

    return tau

def lum_calc(tau, Fdrag, v, dmdt, n_frag_init, fragmented, t, flare_end_time):

    # calculate luminosity (recall tau was set to match optical outputs)
    # based on ablation only - not modeling a shock front
    E_dot_total = (Fdrag * v) + (0.5 * v**2. * abs(dmdt))
    luminosity = tau * E_dot_total

    # apply "pancake" model to flare if fragmented
    if fragmented and t <= flare_end_time:
            luminosity = luminosity * n_frag_init**0.333

    return luminosity

############ Fragmentation logic (now multi-stage) ############
def frag_attempt(bolide_outputs, alt, t, frag_stage, frag_times, frag_alts,
                 current_strength, n_fragments, density, mass, area, q,
                 flare_end_time, flare_duration):

    # Try to fragment if dynamic pressure > current strength
    if q > current_strength:
        frag_stage = frag_stage + 1
        frag_times.append(t)
        frag_alts.append(alt)

        n_fragments = n_fragments**frag_stage

        mass_fragment  = mass / n_fragments
        r_fragment     = ((3. * mass_fragment) / (4. * pi * density))**(1. / 3.)
        dia_fragment   = 2. * r_fragment
        area_fragment  = 4. * pi * (dia_fragment / 2.)**2.
        area_new       = n_fragments * area_fragment
        mass_new       = mass_fragment * n_fragments # conserved
        flare_end_time = t + flare_duration

        print(f"Fragmentation stage {frag_stage} at {alt/1000.:.1f} km, t={t:.3f} s (q={q:.2e} Pa)")
        print(f'flare end time = {flare_end_time:.4e} sec')
        # print(type(flare_end_time))
        bolide_outputs.write(f"Fragmentation stage {frag_stage} at {alt:.4e} m / {t:.4e} s\n\n")

        return True, area_new, mass_new, frag_stage, frag_times, frag_alts, n_fragments, flare_end_time
    else:   
        return False, area, mass, frag_stage, frag_times, frag_alts, n_fragments, flare_end_time

################ Other Functions ###################

# function to calculate initial kinematic quantities
def trajectory_init(bolide_outputs, xstart, zstart, velocity, tstart, theta_deg):

    # initial altitude
    z   = zstart + R_Earth
    R   = z
    alt = zstart
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

    # print(f"ax  = {ax} m/s^2")
    # print(f"az  = {az} m/s^2")
    # print(f"acc = {acc} m/s^2")

    # next calculate final velocities
    vx = vx0 + ax * dt
    vz = vz0 + az * dt
    v = (vx**2. + vz**2.)**0.5
    print(f"updated velocity = {v:.3f} m/s")

    # print(f"vx = {vx} m/s")
    # print(f"vz = {vz} m/s")
    # print(f"v  = {v} m/s")

    # then calculate final positions
    x = x0 + vx0 * dt + 0.5 * ax * dt**2.
    z = z0 + vz0 * dt + 0.5 * az * dt**2.

    dx = x - x0
    dz = z - z0

    d = (dx**2. + dz**2.)**0.5

    # print(f"x = {x} m")
    # print(f"z = {z} m")
    # print(f"d = {d} m")

    # now find new altitude using law of cosines
    R = (R0**2. + d**2. - 2. * R0 * d * math.cos(gamma0))**0.5
    alt = R - R_Earth
    print(f"updated alt = {alt:.3f} m")

    # find delta_alpha and new alpha angle with law of sines
    sindalpha = math.sin(gamma0) * (d / R)
    dalpha = math.asin(sindalpha)
    alpha = alpha0 + dalpha
    # print(f"dalpha = {dalpha * rad2deg} deg")
    # print(f"alpha  = {alpha * rad2deg} deg")

    # find new theta angle
    tantheta = abs(vz) / vx
    theta = math.atan(tantheta)
    dtheta = theta - theta0
    # print(f"theta = {theta * rad2deg} deg")
    
    # find new delta angle
    delta = (pi / 2.) - theta
    # print(f"delta = {delta * rad2deg} deg")

    # find new downrange distance s
    s = C_Earth * (alpha / (2. * pi))
    ds = s - s0
    # print(f"s  = {s} m")
    # print(f"ds = {ds} m")

    # find new flight path angle phi
    phi = (pi / 2.) - (alpha + delta)
    # print(f"phi = {phi * rad2deg} deg")

    # find new gamma angle
    gamma = delta + alpha
    # print(f"gamma = {gamma * rad2deg} deg")

    # Earth surface calcs
    surfx = R_Earth * math.cos((pi / 2.) - alpha)
    surfz = R_Earth * math.sin((pi / 2.) - alpha)
    R_E = (surfx**2. + surfz**2.)**0.5
    # print(f"surfx = {surfx} m")
    # print(f"surfz = {surfz} m")
    # print(f"Radius of Earth = {R_E} m")

    return x, z, R, alt, v, theta, vx, vz, acc, alpha, delta, gamma, phi, s, surfx, surfz

# estimate compressive strength based on density and porosity
def compute_strength(init_strength, porosity):

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

############ Main simulation function ############
def bolide_luminosity_model(outputpath, diameter, velocity, theta_deg,
                            init_strength, strength_secondary, strength_tertiary, strength_quaternary,
                            density, porosity, Cd, L_ablation,
                            n_fragments, n_frag_init, flare_duration, rho_tau_scale,
                            zstart, xstart, tstart, dt, tstop):

    # initialize timestep
    timestep = int(1)
    print(f"first timestep = {timestep}")

    # compute initial area and mass
    area = area_init_calc(diameter)
    mass = mass_init_calc(diameter, density)

    # compute initial material strength of bolide
    strength_initial = compute_strength(init_strength, porosity)
    current_strength = strength_initial

    # initialize fragmentation model
    fragmented     = False  # still useful as flag for first frag
    flare_end_time = None
    frag_stage     = 0
    frag_times     = []
    frag_alts      = []

    # define header
    header = (f"Dia: {diameter:.1f} m | Vel: {velocity/1000.:.1f} km/s | theta: {theta_deg:.1f} deg | "
    f"Strength: {strength_initial:.1e} Pa")

    # open output files
    # general output file
    bolide_outputs = open(f'{outputpath}/bolide_out.dat', 'w')
    bolide_outputs.write(f"Bolide Luminosity and Fragmentation Simulation — 8\n")
    bolide_outputs.write(f"{header}\n\n")
    bolide_outputs.write(f"Initial area     = {area:.4e} m²\n")
    bolide_outputs.write(f"Initial mass     = {mass:.4e} kg\n")
    bolide_outputs.write(f"Initial strength = {strength_initial:.4e} Pa\n\n")

    # power vs time output file
    bolide_pvt = open(f'{outputpath}/power_vs_time.dat', 'w')
    bolide_pvt.write(f'{header}\n')
    bolide_pvt.write(f'Power (W)    Time (sec)\n')

    # power vs altitude output file
    bolide_pva = open(f'{outputpath}/power_vs_alt.dat', 'w')
    bolide_pva.write(f'{header}\n')
    bolide_pva.write(f'Power (W)    Altitude (m)\n')

    # trajectory init (your original function — assuming it exists)
    (x, z, R, alt, v, t, theta, vx, vz, alpha, delta, gamma, phi, s,
    surfx, surfz) = trajectory_init(bolide_outputs, xstart, zstart, velocity, tstart, theta_deg)

    # initialize total radiated and deposited energy to zero
    E_rad_total = 0.0
    E_deposited = 0.0
    print(f"setting total radiated energy to {E_rad_total:.4e} J\n")
    bolide_outputs.write(f"Setting total radiated energy to {E_rad_total:.4e} J\n\n")

    # initialize arrays
    T_array           = np.zeros(1)
    p_array           = np.zeros(1)
    rho_array         = np.zeros(1)
    a_array           = np.zeros(1)
    lum_array         = np.zeros(1)
    alt_array         = np.zeros(1)
    time_array        = np.zeros(1)
    mass_array        = np.zeros(1)
    dia_array         = np.zeros(1)
    area_array        = np.zeros(1)
    KE_array          = np.zeros(1)
    v_array           = np.zeros(1)
    accel_array       = np.zeros(1)
    heatflux_array    = np.zeros(1)
    Fdrag_array       = np.zeros(1)
    E_rad_ttl_array   = np.zeros(1)
    E_deposited_array = np.zeros(1)
    q_array           = np.zeros(1)
    Mach_array        = np.zeros(1)
    T_stag_array      = np.zeros(1)
    p_stag_array      = np.zeros(1)
    T_surf_array      = np.zeros(1)

    # main loop
    while timestep < mxcycl:

        # calculate quantities
        T       = temp_altitude_calc(alt)
        p       = atmos_pressure_calc(alt, T)
        rho     = atmos_density_calc(alt, p, T)
        a       = sound_speed_calc(T, p, rho)
        M       = mach_num_calc(v, a)
        q       = dynamic_pressure_calc(rho, v)
        Fdrag   = Fdrag_calc(rho, v, area, Cd)
        KE      = KE_calc(mass, v)
        T_stag  = T_stag_calc(M, T)
        p_stag  = p_stag_calc(q, p)
        q_h     = heatflux_calc(Fdrag, v, area)
        T_surf  = T_surf_calc(epsilon, q_h)

        # Multi-stage fragmentation check
        (did_frag, area, mass, frag_stage, frag_times, frag_alts,
        n_fragments, flare_end_time) = frag_attempt(bolide_outputs, alt, t, frag_stage, frag_times, frag_alts,
                                    current_strength, n_fragments, density, mass, area, q, flare_end_time, flare_duration)

        if did_frag:
            fragmented = True
            # After first fragmentation → switch to weaker strength
            if frag_stage == 1:
                current_strength = strength_secondary
                print(f"Material strength updated to secondary value: {current_strength:.2e} Pa")
            if frag_stage == 2:
                current_strength = strength_tertiary
                print(f"Material strength updated to tertiary value: {current_strength:.2e} Pa")
            if frag_stage == 3:
                current_strength = strength_quaternary
                print(f"Material strength updated to quaternary value: {current_strength:.2e} Pa")
            if frag_stage > 3:
                current_strength = 2.e7
                print(f'Material strength updated to terminal value of {current_strength} Pa')

        # Luminosity (single term)
        tau        = tau_calc(v, mass, rho, rho_tau_scale, n_fragments if fragmented else 1)
        dmdt       = dmdt_calc(q_h, area, L_ablation)
        luminosity = lum_calc(tau, Fdrag, v, dmdt, n_frag_init, fragmented, t, flare_end_time)

        # Per-fragment reporting
        lum_per_fragment = luminosity / n_fragments if fragmented and n_fragments > 1 else None

        # trajectory update (your original)
        (x, z, R, alt, v, theta, vx, vz, acc, alpha, delta, gamma, phi, s,
        surfx, surfz) = trajectory_update(x, z, R, v, theta, vx, vz, alpha, delta, gamma, phi, s, Fdrag, dt, mass)

        # mass / size update
        mass     = mass_update_calc(mass, dmdt, dt)
        diameter = diameter_update_calc(mass, density, n_fragments if fragmented else 1, fragmented)
        area     = area_update_calc(diameter, n_fragments if fragmented else 1, fragmented)

        # energy accumulation
        E_rad_total = E_rad_total_calc(E_rad_total, luminosity, dt)
        E_deposited = E_deposited_calc(E_deposited, Fdrag, v, dmdt, dt)

        # store to arrays (your original indexing)
        T_array[timestep-1]           = T
        p_array[timestep-1]           = p
        rho_array[timestep-1]         = rho
        a_array[timestep-1]           = a
        lum_array[timestep-1]         = luminosity
        alt_array[timestep-1]         = alt
        time_array[timestep-1]        = t
        mass_array[timestep-1]        = mass
        dia_array[timestep-1]         = diameter
        area_array[timestep-1]        = area
        KE_array[timestep-1]          = KE
        v_array[timestep-1]           = v
        accel_array[timestep-1]       = acc
        heatflux_array[timestep-1]    = q_h
        Fdrag_array[timestep-1]       = Fdrag
        E_rad_ttl_array[timestep-1]   = E_rad_total
        E_deposited_array[timestep-1] = E_deposited
        q_array[timestep-1]           = q
        Mach_array[timestep-1]        = M
        T_stag_array[timestep-1]      = T_stag
        p_stag_array[timestep-1]      = p_stag
        T_surf_array[timestep-1]      = T_surf

        # advance time
        t = t + dt

        # write outputs
        wr_out(bolide_outputs, luminosity, z, t, mass, diameter, area, KE, v, acc, q_h, Fdrag,
               E_rad_total, E_deposited, q, M, T_stag, p_stag, T_surf, T, p, rho, a, theta_deg, tau,
               frag_stage, lum_per_fragment)
        
        wr_power_vs_time(bolide_pvt, luminosity, t)

        wr_power_vs_alt(bolide_pva, luminosity, alt)

        # end conditions
        if t > tstop:
            print('STOP simulation: time exceeded end time\n\n')
            break
        if mass <= 0.:
            print('STOP simulation: mass <= 0 kg\n\n')
            break
        if alt <= 0.:
            print('STOP simulation: altitude <= 0 km\n\n')
            break
        if v <= 0.:
            print('STOP simulation: v <= 0 m/s\n\n')
            break

        # extend arrays
        T_array           = np.vstack([T_array, new_row])
        p_array           = np.vstack([p_array, new_row])
        rho_array         = np.vstack([rho_array, new_row])
        a_array           = np.vstack([a_array, new_row])
        lum_array         = np.vstack([lum_array, new_row])
        alt_array         = np.vstack([alt_array, new_row])
        time_array        = np.vstack([time_array, new_row]) 
        mass_array        = np.vstack([mass_array, new_row])
        dia_array         = np.vstack([dia_array, new_row])
        area_array        = np.vstack([area_array, new_row])
        KE_array          = np.vstack([KE_array, new_row])
        v_array           = np.vstack([v_array, new_row])
        accel_array       = np.vstack([accel_array, new_row])
        heatflux_array    = np.vstack([heatflux_array, new_row])
        Fdrag_array       = np.vstack([Fdrag_array, new_row])
        E_rad_ttl_array   = np.vstack([E_rad_ttl_array, new_row])
        E_deposited_array = np.vstack([E_deposited_array, new_row])
        q_array           = np.vstack([q_array, new_row])
        Mach_array        = np.vstack([Mach_array, new_row])
        T_stag_array      = np.vstack([T_stag_array, new_row])
        p_stag_array      = np.vstack([p_stag_array, new_row])
        T_surf_array      = np.vstack([T_surf_array, new_row])

        # update timestep
        timestep = timestep + 1
        print(f'timestep updated: new time will be {t:.3e} sec\n\n')

    # final stats
    peak_power, peak_rad_i = peak_find(bolide_outputs, lum_array)
    E_event                = E_event_calc(bolide_outputs, E_rad_total, E_deposited) # your original

    # print out if bolide fragmented or not
    if fragmented == True:
        print(f'bolide fragmented at {frag_alts[0]/1000.:.1f} km and {frag_times[0]:.2f} sec')
    else:
        print(f'bolide did not fragment')

    bolide_outputs.close()
    bolide_pvt.close()
    bolide_pva.close()

    return (header, strength_initial, T_array, p_array, rho_array, a_array, lum_array,
            alt_array, time_array, mass_array, dia_array, area_array,
            frag_times, frag_alts, KE_array, v_array, accel_array, heatflux_array,
            Fdrag_array, E_rad_ttl_array, q_array, Mach_array, T_stag_array,
            p_stag_array, T_surf_array, peak_power, peak_rad_i, E_rad_total, E_event)

##############################################################

# Main Program
event_name = input('Enter event name: ')
print(f'Event name input: {event_name}')
inputpath  = f"{basename}{event_name}/bolum_in.dat"
outputpath = f"{basename}{event_name}"

# read inputs
(diameter, velocity, theta_deg,
init_strength, strength_secondary, strength_tertiary, strength_quaternary,
density, porosity, Cd, L_ablation,
n_fragments, n_frag_init, flare_duration, rho_tau_scale,
zstart, xstart, tstart, dt, tstop) = rdinput(inputpath)

# write out initial parameters
write_init_params(outputpath, diameter, velocity, theta_deg,
init_strength, strength_secondary, strength_tertiary, strength_quaternary,
density, porosity, Cd, L_ablation,
n_fragments, n_frag_init, flare_duration, rho_tau_scale,
zstart, xstart, tstart, dt, tstop)

# run simulation
(header, strength, T_array, p_array, rho_array, a_array,
lum_array, alt_array, time_array, mass_array, dia_array,
area_array, frag_times, frag_alts, KE_array, v_array,
accel_array, heatflux_array, Fdrag_array, E_rad_ttl_array,
q_array, Mach_array, T_stag_array, p_stag_array, T_surf_array,
peak_power, peak_rad_i, E_rad_tot, E_event) = bolide_luminosity_model(outputpath, diameter, velocity, theta_deg,
                                                                      init_strength, strength_secondary, strength_tertiary, strength_quaternary,
                                                                      density, porosity, Cd, L_ablation,
                                                                      n_fragments, n_frag_init, flare_duration, rho_tau_scale,
                                                                      zstart, xstart, tstart, dt, tstop)
                                                             
# end simulation and plot output
# error checking
print(f"total cycles:     {time_array.size}\n")

print(f"STOP all done bolide simulation complete\n")
plot_outputs(outputpath, header, time_array.flatten(), lum_array.flatten(), alt_array.flatten(),
             frag_times, frag_alts, v_array.flatten(), accel_array.flatten(), q_array.flatten(),
             T_stag_array.flatten(), p_stag_array.flatten(), T_surf_array.flatten(),
             E_rad_ttl_array.flatten(), peak_power, peak_rad_i, E_rad_tot, E_event)