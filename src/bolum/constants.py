'''
bolum constants file
'''

import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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

# other simulation parameters
mxcycl  = 100000           # maximum number of cycles
new_row = np.zeros(1)