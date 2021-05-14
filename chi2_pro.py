import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from event_rate import *
from scipy.optimize import fsolve, brentq

years = 10.
chi2_th = 2.71

# Poisson chi2
# s: signal, b: background, d: data
def Chi_sq_bin(s,b,d):
    return 2.*( s + b -d + d*np.log(d/(s+b)) )

# Compute the chi2 for a grid of PBH parameters for a given experiment
def compute_chi2_2D(Mpbhs, fpbhs, exp):

    data_final = []

    Eback, eventback = back_rate(exp)
    eventback = eventback*years
    backint = interp1d(Eback, eventback)
    Ebin = Eback    # Binear?

    """Ebin1, Ebin2 = energybins(exp)
    Ebin = np.linspace(Ebin1, Ebin2)
    eventback = backint(Ebin)"""

    """if exp=="SK":
        Edat, eventdat = EdatSK, datSK
    else:
        Edat, eventdat = Eback, eventback"""

    eventdat = eventback    # for forecasts, take data as background

    for Mpbh in Mpbhs:

        folder = "folder_fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh, exp)
        Evec, events = np.loadtxt(folder, unpack=True)
        evint = interp1d(Evec, events)

        for fpbh in fpbhs:

            signal = fpbh*events*years

            chi2_tot = 0

            for idx, element in enumerate(Ebin):
                chi2_tot = chi2_tot + Chi_sq_bin(signal[idx],eventback[idx],eventdat[idx])

            data_final.append([Mpbh,fpbh,chi2_tot])

    np.savetxt("data/chi2_PHB_"+exp+".txt",data_final)
    return data_final


# Compute the chi2 for a grid of PBH parameters for a given experiment
def compute_chi2_2D_mod(Mpbhs, fpbhs, exp):

    data_final = []

    Eback, eventback = back_rate(exp)
    eventback = eventback*years
    backint = interp1d(Eback, eventback)
    Ebin = Eback    # Binear?

    """Ebin1, Ebin2 = energybins(exp)
    Ebin = np.linspace(Ebin1, Ebin2)
    eventback = backint(Ebin)"""

    """if exp=="SK":
        Edat, eventdat = EdatSK, datSK
    else:
        Edat, eventdat = Eback, eventback"""

    eventdat = eventback    # for forecasts, take data as background

    fpbh_bounds = []

    for Mpbh in Mpbhs:

        folder = "folder_fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh, exp)
        Evec, events = np.loadtxt(folder, unpack=True)
        evint = interp1d(Evec, events)

        chi2_fpbh = []

        for fpbh in fpbhs:

            signal = fpbh*events*years

            chi2_tot = 0

            for idx, element in enumerate(Ebin):
                chi2_tot += Chi_sq_bin(signal[idx],eventback[idx],eventdat[idx])

            chi2_fpbh.append(chi2_tot)

            data_final.append([Mpbh,fpbh,chi2_tot])


        chi2int = interpolate.interp1d(fpbhs, chi2_fpbh, fill_value="extrapolate")
        fpbhvec = np.logspace(np.log10(fpbhs[0]), np.log10(fpbhs[-1]))
        minchi2 = np.amin(chi2int(fpbhvec))
        fpbh_bounds.append( fsolve( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), 1.e-2  ) )
        #fpbh_bounds.append( brentq( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), fpbhs[0], fpbhs[-1]  ) )

    np.savetxt("data/chi2_PHB_"+exp+".txt",data_final)
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
