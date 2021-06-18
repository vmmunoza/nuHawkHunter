#-------------------------------
# Statistical functions for computing the bounds
#-------------------------------

import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from Source.event_rate import *
from scipy.optimize import fsolve, brentq

years = 10.
#chi2_th = 1.    # 68.27 CL (1 sigma)
chi2_th = 2.71  # 90 % CL
#chi2_th = 3.84  # 95 % CL
#chi2_th = 4.    # 95.45 % CL (2 sigma)

# Poisson chi2
# s: signal, b: background, d: data
def Chi_sq_bin(s,b,d):
    return 2.*( s + b -d + d*np.log(d/(s+b)) )
    #return  (s + b - d)**2./b   # Gaussian approximation

"""
# Compute the chi2 for a grid of PBH parameters for a given experiment
def compute_chi2_2D(Mpbhs, fpbhs, exp):

    data_final = []

    Eback, eventback = back_rate(exp)
    eventback = eventback*years

    #Ebin1, Ebin2 = energybins(exp)
    #Ebin = np.linspace(Ebin1, Ebin2)
    #eventback = backint(Ebin)

    #if exp=="SK":
    #    Edat, eventdat = EdatSK, datSK
    #else:
    #    Edat, eventdat = Eback, eventback

    eventdat = eventback    # for forecasts, take data as background

    for Mpbh in Mpbhs:

        folder = "fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh, exp)
        Evec, events = np.loadtxt(folder, unpack=True)

        for fpbh in fpbhs:

            signal = fpbh*events*years

            chi2_tot = 0

            for idx, element in enumerate(Ebin):
                chi2_tot = chi2_tot + Chi_sq_bin(signal[idx],eventback[idx],eventdat[idx])

            data_final.append([Mpbh,fpbh,chi2_tot])

    np.savetxt("data/chi2_PHB_"+exp+".txt",data_final)
    return data_final
"""

# Compute the chi2 for a grid of PBH parameters for a given experiment
def compute_chi2_2D_mod(Mpbhs, fpbhs, exp, is_DM=True):

    data_final = []

    Eback, eventback = back_rate(exp)
    bin = 1.
    if (exp=="DARWIN") or (exp=="ARGO"):
        bin = 5.e-3
    Eback, eventback = binned_events(Eback, eventback, bin)
    eventback = eventback*years

    eventdat = eventback    # for forecasts, take data as background

    fpbh_bounds = []

    for Mpbh in Mpbhs:

        folder = "fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh, exp)
        Evec, events = np.loadtxt(folder, unpack=True)

        chi2_fpbh = []

        for fpbh in fpbhs:

            signal = fpbh*events*years

            chi2_tot = 0.

            for i, element in enumerate(Eback):
                chi2_tot += Chi_sq_bin(signal[i], eventback[i], eventdat[i])

            chi2_fpbh.append(chi2_tot)

            data_final.append([Mpbh,fpbh,chi2_tot])


        #chi2int = interpolate.interp1d(np.log10(fpbhs), chi2_fpbh, fill_value="extrapolate")
        chi2int = interpolate.interp1d(fpbhs, chi2_fpbh, fill_value="extrapolate")
        #fpbhvec = np.logspace(np.log10(fpbhs[0]), np.log10(fpbhs[-1]))
        minchi2 = 0. #np.amin(chi2int(np.log10(fpbhvec)))

        if is_DM:
            fpbh_bounds.append( fsolve( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), 1.e-2  ) )
            #fpbh_bounds.append( 10.**fsolve( lambda logfpbh: chi2int(logfpbh) - (minchi2 + chi2_th), -2  ) )
            #fpbh_bounds.append( brentq( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), fpbhs[0], 1.e5  ) )
            #fpbh_bounds.append( 10.**brentq( lambda logfpbh: chi2int(logfpbh) - (minchi2 + chi2_th), np.log10(fpbhs[0]), 2.  ) )
        else:
            fpbh_bounds.append( fsolve( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), 1.e-16  ) )

    np.savetxt("outputs/chi2_PHB_"+exp+".txt",data_final)
    return data_final, fpbh_bounds


def grid_val(ptx,pty,ptz):

    my_df = pd.DataFrame.from_dict({'y':pty, 'x':ptx, 'fxy': ptz})

    def bivariate_interpolation(df, FXY= 'fxy', Y='y', X='x'):
        df_copy = df.copy().sort_values(by=[Y, X], ascending=True)
        x = np.array(df_copy[X].unique(), dtype='float64')
        y = np.array(df_copy[Y].unique(), dtype='float64')
        z = np.array(df_copy[FXY].values, dtype='float64')
        Z = z.reshape(len(y), len(x))

        interp_spline = interpolate.RectBivariateSpline(y, x, Z, kx=1, ky=1)
        #x, y = np.meshgrid(x, y)
        #interp_spline = interpolate.Rbf(y, x, Z, kind="linear")
        return interp_spline


    interpolador= bivariate_interpolation(my_df)

    def chi2_func(y,x):
        return interpolador(y, x, grid=False)

    FunVec=np.vectorize(chi2_func)

    y_value=np.unique(pty)
    x_value=np.unique(ptx)
    X_grid, Y_grid = np.meshgrid(x_value, y_value)
    Z_grid= FunVec(Y_grid,X_grid)

    return X_grid, Y_grid, Z_grid

# Bound from Poisson chi2 for one single bin
def bound_anal_chi2_Poi(back, sig):
    #sig is with fpbh=1
    chi2_th = 2.71
    snr = sig/back
    return brentq( lambda fpbh: 2.71/(2.*back) - (fpbh*snr - np.log(1. + fpbh*snr)), 1.e-7, 1.e5  )

# Bound from approximated chi2 for one single bin for low snr
def bound_anal_chi2_apr(back, sig):
    #sig is with fpbh=1
    chi2_th = 2.71
    snr = sig/back
    return np.sqrt(chi2_th(back))/snr
    #return brentq( lambda fpbh: chi2_th/(back) - (fpbh*snr)**2., 1.e-7, 1.e5  )
