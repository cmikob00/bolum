import numpy as np
import matplotlib.pyplot as plt

# Constants
g = 9.81  # gravity (m/s^2)
R = 287.05         # J/kg·K, specific gas constant for dry air
sigma_sb = 5.670374419e-8  # Stefan–Boltzmann constant
gamma = 1.4        # adiabatic index for air
rho_air0 = 1.225  # sea-level atmospheric density (kg/m^3)
H = 7160.0  # scale height of atmosphere (m)
lum_efficiency = 0.03  # luminous efficiency
strength_baseline = 1.e6  # baseline strength in Pascals
epsilon = 0.9      # surface emissivity

def atmospheric_density(h):
    # calculate atmospheric density
    rho = rho_air0 * np.exp(-h / H)

    return rho

def atmospheric_pressure(h):
    # simple exponential model for pressure
    p = 101325 * np.exp(-h / H)

    return p

def temperature_at_altitude(h):
    # Linear lapse rate model: T = T0 - L * h
    T0 = 288.15  # K
    L = 0.0065   # K/m

    return max(T0 - L * h, 200)  # floor to avoid unrealistic temps

def dynamic_pressure(rho, v):
    # calculate dynamic pressure
    q = 0.5 * rho * v**2

    return q

def compute_strength(density, porosity):
    # estimate compressive strength based on density and porosity.
    strength = strength_baseline * (density / 3000.0) * (1.0 - porosity)

    print(f"strength of bolide is {strength} Pa")

    return strength

def bolide_luminosity_model(diameter, velocity, angle_deg, density, porosity, dt, n_fragments, flare_duration):
    # printing out initial conditions
    angle = np.radians(angle_deg)
    print(f"initial angle of entry in radians = {angle} radians")
    area = np.pi * (diameter / 2)**2
    print(f"area of bolide = {area} m^2")
    mass = (4/3) * np.pi * (diameter / 2)**3 * density
    print(f"mass of bolide = {mass} kg")

    h = 100.e3  # initial altitude [m]
    print(f"initial altitude = {h} meters")
    v = velocity
    print(f"initial velocity = {v} m/s")
    t = 0.0
    print(f"initial time = {t} seconds")

    # initializing arrays
    lum_curve = []
    altitudes = []
    KE_array = []
    Mach_array = []
    T_stag_array = []
    p_stag_array = []
    T_surf_array = []
    accel_array = []
    time_array = []

    fragmented = False
    frag_altitude = None
    frag_time = None
    flare_end_time = None

    strength = compute_strength(density, porosity)

    while h > 0 and v > 0:
        rho = atmospheric_density(h)
        p = atmospheric_pressure(h)
        T = temperature_at_altitude(h)
        a = np.sqrt(gamma * R * T)
        q = dynamic_pressure(rho, v)

        if not fragmented and q > strength:
            frag_altitude = h
            fragmented = True
            frag_time = t
            flare_end_time = t + flare_duration
            area = area * n_fragments
            mass = mass / n_fragments

        drag = 0.5 * rho * v**2. * area
        acc = drag / mass

        dv = acc * dt
        dh = -v * np.sin(angle) * dt

        dE = drag * v * dt
        luminosity = lum_efficiency * dE / dt

        if fragmented and t <= flare_end_time:
            luminosity = luminosity * 10

        # Store all arrays
        altitudes.append(h)
        lum_curve.append(luminosity)
        KE_array.append(0.5 * mass * v**2)
        Mach_array.append(v / a)
        T_stag_array.append(T * (1 + 0.5 * (gamma - 1) * (v / a)**2))
        p_stag_array.append(p + q)
        T_surf_array.append((drag / (area * epsilon * sigma_sb))**0.25)
        accel_array.append(acc)
        time_array.append(t)

        v = v - dv
        print(f"updated velocity = {v} m/s")
        h = h + dh
        print(f"updated altitude = {h} meters")
        t = t + dt
        print(f"updated time = {t} seconds")

        if t > 60:
            break

    return (
        np.array(altitudes),
        np.array(lum_curve),
        frag_altitude,
        frag_time,
        np.array(KE_array),
        np.array(Mach_array),
        np.array(T_stag_array),
        np.array(p_stag_array),
        np.array(T_surf_array),
        np.array(accel_array),
        np.array(time_array)
    )

# Inputs
diameter = 230.0        # meters
velocity = 11200.0    # m/s
angle = 90.0          # degrees
density = 3000.0      # kg/m^3
porosity = 0.2        # 0 to 1
dt = 0.01             # seconds
n_fragments = 10      # number of fragments produced
flare_duration = 0.2  # seconds

alts, lums, frag_h, frag_t, KEs, Machs, T_stags, p_stags, T_surfs, accels, times = bolide_luminosity_model(
    diameter, velocity, angle, density, porosity, dt, n_fragments, flare_duration)

# Plot
plt.figure()
plt.plot(lums, alts / 1000.)
plt.xscale('log')
plt.xlabel('Luminosity (Watts)')
plt.ylabel('Altitude (km)')
plt.axhline(frag_h / 1000., linestyle='--', color='green')
plt.title(f'Bolide Light Curve\nFragmentation at {frag_h/1000:.1f} km' if frag_h else 'No Fragmentation')
plt.grid(True)
plt.show()

plt.figure()
plt.plot(times, lums)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time (seconds)')
plt.ylabel('Luminosity (Watts)')
plt.axvline(frag_t, linestyle='--', color='green')
plt.title(f'Bolide Light Curve\nFragmentation at {frag_h/1000:.1f} km' if frag_h else 'No Fragmentation')
plt.grid(True)
plt.show()
