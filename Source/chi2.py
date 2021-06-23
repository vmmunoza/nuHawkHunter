#-------------------------------
# Statistical functions for computing the bounds
#-------------------------------

from Source.event_rate import *
from scipy.optimize import fsolve, brentq

# Years for the exposure
years = 10.

# Chi2 treshold
#chi2_th = 1.    # 68.27 CL (1 sigma)
chi2_th = 2.71  # 90 % CL
#chi2_th = 3.84  # 95 % CL
#chi2_th = 4.    # 95.45 % CL (2 sigma)

# Poisson chi2
# s: signal, b: background, d: data
def chi2_Poisson(s,b,d):
    return 2.*( s + b -d + d*np.log(d/(s+b)) )
    #return  (s + b - d)**2./b   # Gaussian approximation

# From a list of PBH masses, compute the chi2 and the bound on the abundance for a given experiment
def pbh_bounds_chi2(Mpbhs, fpbhs, exp, is_DM=True):

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
                chi2_tot += chi2_Poisson(signal[i], eventback[i], eventdat[i])

            chi2_fpbh.append(chi2_tot)


        #chi2int = interp1d(np.log10(fpbhs), chi2_fpbh, fill_value="extrapolate")
        chi2int = interp1d(fpbhs, chi2_fpbh, fill_value="extrapolate")
        #fpbhvec = np.logspace(np.log10(fpbhs[0]), np.log10(fpbhs[-1]))
        minchi2 = 0. #np.amin(chi2int(np.log10(fpbhvec)))

        # If is DM, obtain f_pbh, else obtain \beta'
        if is_DM:
            fpbh_bounds.append( fsolve( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), 1.e-2  ) )
            #fpbh_bounds.append( 10.**fsolve( lambda logfpbh: chi2int(logfpbh) - (minchi2 + chi2_th), -2  ) )
            #fpbh_bounds.append( brentq( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), fpbhs[0], 1.e5  ) )
            #fpbh_bounds.append( 10.**brentq( lambda logfpbh: chi2int(logfpbh) - (minchi2 + chi2_th), np.log10(fpbhs[0]), 2.  ) )
        else:
            fpbh_bounds.append( fsolve( lambda fpbh: chi2int(fpbh) - (minchi2 + chi2_th), 1.e-16  ) )

    return fpbh_bounds

# Bound from Poisson chi2 for one single bin
def bound_anal_chi2_Poi(back, sig):
    # sig is with fpbh=1
    snr = sig/back
    return brentq( lambda fpbh: chi2_th/(2.*back) - (fpbh*snr - np.log(1. + fpbh*snr)), 1.e-7, 1.e5  )

# Bound from approximated chi2 for one single bin, valid for low snr
def bound_anal_chi2_apr(back, sig):
    #sig is with fpbh=1
    snr = sig/back
    return np.sqrt(chi2_th/back)/snr
    #return brentq( lambda fpbh: chi2_th/(back) - (fpbh*snr)**2., 1.e-7, 1.e5  )
