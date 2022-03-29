#-------------------------------
# Constant definition and initialization
#-------------------------------

import os
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib.lines import Line2D
import matplotlib as mpl

# Settings for plots
fontsiz = 20
mpl.rcParams.update({'font.size': fontsiz})
mpl.rcParams['xtick.labelsize'] = fontsiz
mpl.rcParams['ytick.labelsize'] = fontsiz
plt.rcParams['figure.figsize'] = [12, 8]
mpl.rcParams['xtick.major.size'] = fontsiz/2
mpl.rcParams['ytick.major.size'] = fontsiz/2
mpl.rcParams['xtick.minor.size'] = fontsiz/4
mpl.rcParams['ytick.minor.size'] = fontsiz/4


# Create some folders if needed
for path in ["figures", "fluxes"]:
    if not os.path.exists(path):
        os.system("mkdir "+path)

#-------
# Parameters
#-------

# Cosmological parameters
Om_rad = 9.2e-5 # Omega radiation
Om_m = 0.315    # Omega matter
Om_dm = 0.265   # Omega dark matter
Om_L = 0.685    # Omega cosmological constant
year_sec = 3.154e7 # 1yr in seconds
hlittle = 0.68      # Reduced Hubble constant
H0 = hlittle*3.2407e-18 # Hubble constant, s^-1
G = 6.67259e-8  # Gravitational constant, cm^3 g^-1 s^-2
rho_c = 3.*H0**2./(8.*np.pi*G) # Critical density, g cm^-3
zeq = 2.4e4*Om_m*hlittle**2. -1. # Matter-radiation equilibrium redshift

# Unit conversions
Msun = 1.989e33 # Solar mass, g
c = 29979245800.0 # Speed of light, cm/s
MpcToCm = 3.086e24 # 1 Mpc = MpcToCm cm
GpToCm = 1.e3*MpcToCm
gr_to_GeV= 5.60958884e+23  # g to GeV
time_conversion = 1.519267407e+24 # s to GeV^(-1)
leng_conversion = 5.06773058e+13 # cm to GeV^(-1)
G_N = (6.67408e-11*pow(leng_conversion*100.,3.)/(gr_to_GeV*1000.)/pow(time_conversion,2.)) # G_N for Black-body approximation

# Colors for plots
cols = ["r", "m", "b", "g","m", "orange"]

# Inital PBH mass (in g) evaporating at the age of the universe
Mevap = 7.8e14  # estimation from BlackHawk results (other estimates give 5.e14)

# Folder to place BlackHawk spectra
folder_blackhawk = "BlackHawkData/"

#-------
# Miscelanea functions
#-------

# Write in Latex scientific notation
def scinot(x):
    exp = int(math.floor(math.log10(abs(x))))
    prefactor = x / 10**exp
    if prefactor == 1.:
        if exp==0.:
            return r"$1$"
        else:
            return r"$10^{"+str(exp)+r"}$"
    else:
        return "{:.1f}".format(prefactor)+r"$ \times 10^{"+str(exp)+r"}$"

# Given a value, find its nearest element from an array
def find_nearest(array, value, axis=0):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin(axis)
    return idx
