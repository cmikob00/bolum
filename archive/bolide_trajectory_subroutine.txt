# bolide trajectory subroutine

import math

# trajectory calculations

# initial values
alt   = 50000.
alpha = 3.
dt    = 0.1
v     = 20000.
theta = 15.

# constants
R_Earth = 6378.14e3
C_Earth = 2. * math.pi * R_Earth
rad2deg = 180 / math.pi
deg2rad = math.pi / 180
g       = 9.81 * (R_Earth / (R_Earth + alt))**2.
drag    = 5.

# more initial values
alpha = alpha * deg2rad
x = (R_Earth + alt) * math.sin(alpha)
z = (R_Earth + alt) * math.cos(alpha)
print('alpha', alpha)
print('x', x)
print('z', z)

s = C_Earth * (alpha / (2. * math.pi))
print('s', s)

R = (x**2. + z**2.)**0.5
print('R', R)

theta = theta * deg2rad
delta = (math.pi / 2.) - theta
gamma = delta + alpha
phi = (math.pi / 2.) - (alpha + delta)

# calc velocity components
vx = v * math.cos(theta)
vz = -v * math.sin(theta)

# calc accels
gx = -g * math.sin(alpha)
gz = -g * math.cos(alpha)

dragx = -drag * math.cos(theta)
dragz = drag * math.sin(theta)

ax = gx + dragx
az = gz + dragz

a = (ax**2. + az**2.)**0.5

print('ax', ax)
print('az', az)
print('a', a)

# vel calcs
v0x = vx
v0z = vz

vx = v0x + ax * dt
vz = v0z + az * dt

v = (vx**2. + vz**2.)**0.5

print('vx', vx)
print('vz', vz)
print('v', v)

# pos calcs
x0 = x
z0 = z

x = x0 + v0x * dt + 0.5 * ax * dt**2.
z = z0 + v0z * dt + 0.5 * az * dt**2.

dx = x - x0
dz = z - z0

d = (dx**2. + dz**2.)**0.5

print('x', x)
print('z', z)
print('d', d)

# find new alt using law of cosines
alt0 = alt
R0 = R
delta0 = delta
gamma0 = gamma

R = (R0**2. + d**2. - 2. * R0 * d * math.cos(gamma))**0.5

alt = R - R_Earth
print('alt', alt)

# find delta alpha with law of sines
alpha0 = alpha
sindalpha = math.sin(gamma) * (d / R)
dalpha = math.asin(sindalpha)
alpha = alpha0 + dalpha
print('dalpha', dalpha)
print('alpha', alpha)

# calc final values
# new theta - v vector angle wrt horizontal
theta0 = theta
tantheta = abs(vz) / vx
theta = math.atan(tantheta)
print(theta)
print('theta', theta * rad2deg)

dtheta = theta - theta0
print('dtheta', dtheta)

# new s - downrange distance
s0 = s
s = C_Earth * (alpha / (2 * math.pi))
ds = s - s0
print('s', s)
print('ds', ds)

# new phi - flight path angle
phi0 = phi
phi = (math.pi / 2.) - (alpha + delta)
print('phi', phi * rad2deg)

# Earth surface calcs
surfx = R_Earth * math.cos((math.pi / 2.) - alpha)
surfz = R_Earth * math.sin((math.pi / 2.) - alpha)
print('surfx', surfx)
print('surfz', surfz)