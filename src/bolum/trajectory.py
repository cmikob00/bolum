'''
Bolum trajectory model file
'''

from bolum.constants import *

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
    R_E = (surfx**2. + surfz**2.)**0.5
    bolide_outputs.write(f"Initial Earth surface x-comp = {surfz:.4e} meters\n")
    bolide_outputs.write(f"Initial Earth surface z-comp = {surfz:.4e} meters\n")
    bolide_outputs.write(f"Earth Radius                 = {R_E:.4e} meters\n\n")

    return x, z, R, alt, v, t, theta, vx, vz, alpha, delta, gamma, phi, s, surfx, surfz, R_E

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
    R_E   = (surfx**2. + surfz**2.)**0.5
    # print(f"surfx = {surfx} m")
    # print(f"surfz = {surfz} m")
    # print(f"Radius of Earth = {R_E} m")

    return x, z, R, alt, v, theta, vx, vz, acc, alpha, delta, gamma, phi, s, surfx, surfz, R_E