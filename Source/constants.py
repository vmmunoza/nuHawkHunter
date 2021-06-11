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
mpl.rcParams.update({'font.size': 12})

# Create some folders if needed
for path in ["figures", "folder_fluxes"]:
    if not os.path.exists(path):
        os.system("mkdir "+path)

#-------
# Parameters
#-------

# Cosmological parameters
Om_rad = 9.2e-5
Om_m = 0.315
Om_dm = 0.265
Om_L = 0.685
year_sec = 3.154e7 # 1yr in seconds
hlittle = 0.68
H0 = hlittle*3.2407e-18 # s^-1
G = 6.67259e-8  # cm^3 g^-1 s^-2
rho_c = 3.*H0**2./(8.*np.pi*G) # g cm^-3
zeq = 2.4e4*Om_m*hlittle**2. -1. # Matter-radiation equilibrium redshift

# Unit conversions
Msun = 1.989e33 # g
c = 29979245800.0 #cm/s
MpcToCm = 3.086e24 # 1 Mpc = MpcToCm cm
GpToCm = 1.e3*MpcToCm
gr_to_GeV= 5.62e23
mass_conversion = 5.60958884e+23 # g to GeV
time_conversion = 1.519267407e+24 # s to GeV^(-1)
leng_conversion = 5.06773058e+13 # cm to GeV^(-1)
G_N = (6.67408e-11*pow(leng_conversion*100.,3.)/(mass_conversion*1000.)/pow(time_conversion,2.)) # G_N for Black-body approximation
sectoyear = 3.17098e-8

# Colors for plots
cols = ["r", "m", "purple", "b", "g", "orange"]

# Inital PBH mass (in g) evaporating at the age of the universe
Mevap = 4.e14

#-------
# Miscelanea functions
#-------

def scinot(x):      # write in Latex scientific notation
    exp = int(math.floor(math.log10(abs(x))))
    prefactor = x / 10**exp
    if prefactor == 1.:
        if exp==0.:
            return r"$1$"
        else:
            return r"$10^{"+str(exp)+r"}$"
    else:
        return "{:.1f}".format(prefactor)+r"$ \times 10^{"+str(exp)+r"}$"

def find_nearest(array, value, axis=0):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin(axis)
    return idx
