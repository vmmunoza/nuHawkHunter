#-------------------------------
# Cross section definition
#-------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


# IBD cross section constants, from Strumia and Vissani 2003
m_e = 0.511  # MeV
m_p = 938.272 # MeV
m_n = 939.565   # MeV
np_dif =  m_n - m_p # 1.293 MeV
Mnp = (m_p + m_n)/2. # 938.9 MeV
G_F = 1.166e-5  # GeV^-2
costhetaC = 0.9746  # cosine of Cabibbo angle
M_V = np.sqrt(0.71) # MeV
M_A = 1.    # MeV
m_pi = 139. # MeV
xi = 3.706
g1_0 = -1.270
Eth = ((m_n + m_e)**2. - m_p**2.)/(2.*m_p)  # nu threshold in cross section, 1.8057 MeV
delta = (m_n**2. - m_p**2. - m_e**2.)/(2.*m_p)

E_nu_th = np_dif + m_e   # nu threshold 1.804 MeV


# From Strumia and Vissani 2003
def dsigmadE_IBD(E_nu, E_e):

    # Mandelstam variables
    s = 2.*m_p*E_nu + m_p**2.
    u = s + m_e**2. - 2.*m_p*(E_nu + E_e)
    t = m_n**2. - m_p**2. - 2.*m_p*(E_nu - E_e)

    # f and g dimensionless form factors
    f1 = ( 1.- (1.+xi)*t/(4.*Mnp**2.) )/(1. - t/(4.*Mnp**2.) )/(1. - t/(M_V**2.) )**2.
    f2 = xi/(1. - t/(4.*Mnp**2.) )/(1. - t/(M_V**2.) )**2.
    g1 = g1_0/(1. - t/(M_A**2.) )**2.
    g2 = 2.*Mnp**2.*g1/(m_pi**2. - t)

    # Complete expression
    #A = (t-m_e)**2./16.*( 4.*f1**2.*(4.*Mnp**2. + t + m_e**2.) + 4.*g1**2.*(-4.*Mnp**2. + t + m_e**2.) + f2**2.*(t**2./Mnp**2. + 4.*t + 4.*m_e**2.) + 4.*m_e**2.*t*g2**2./Mnp**2. + 8. )
    #B =
    #C =

    # NLO approx, Strumia and Vissani 2003 eqs. 10
    A = Mnp**2.*( f1**2. - g1**2. )*(t - m_e**2.) - Mnp**2.*np_dif**2.*( f1**2. + g1**2. ) -2.*m_e**2.*Mnp*np_dif*g1*( f1 + f2 )
    B = t*g1*(f1+f2)
    C = ( f1**2. + g1**2. )/4.

    """
    EnuCM = (s-m_p**2.)/(2.*np.sqrt(s))
    EeCM = (s - m_n**2. + m_e**2.)/(2.*np.sqrt(s))
    peCM = np.sqrt( (s - ( m_n - m_e )**2.)*( s - ( m_n + m_e )**2. ) )/(2.*np.sqrt(s))
    #peCM[peCM == np.nan] = 1.e-1

    E_1 = E_nu - delta - EnuCM*(EeCM + peCM)/m_p
    E_2 = E_nu - delta - EnuCM*(EeCM - peCM)/m_p
    #print(E_1, E_2)
    """

    Msquare = A - (s-u)*B + (s-u)**2.*C
    dsigmadt = G_F**2.*costhetaC**2./(2.*np.pi*(s-m_p**2.)**2.)*Msquare
    # Allowed range of energies for E_nu and E_e
    all_range = 1.#np.heaviside( E_nu - Eth, 0. )#*np.heaviside( E_e - E_1, 0. )*np.heaviside( E_2 - E_e, 0. )
    dsigmadEe = 2.*m_p*dsigmadt*all_range

    return dsigmadEe

"""
E_ovec = np.linspace(2.,50)
plt.loglog(E_ovec, dsigmadE_IBD(E_ovec, 1.))
plt.loglog(E_ovec, dsigmadE_IBD(E_ovec, 10.))
plt.loglog(E_ovec, dsigmadE_IBD(E_ovec, 20.))

plt.show()
exit()
"""

# From 1712.06985, in cm^2, correction from Beacom DSNB review
def sigmaIBD(E_e):
    return 9.52e-44*( E_e*np.sqrt(E_e**2. - m_e**2.) )*(1. - 7.*(E_e + np_dif)/m_p )

# nu_e Argon cross section for DUNE
# From Denton and Suliga mail
EEAr, sigAr = np.loadtxt("data/crosssections/XS_nue_40Ar.txt", unpack=True, delimiter=";")
sigmaAr = interp1d(EEAr, sigAr, fill_value="extrapolate")

# nu_ebar Carbon cross section for JUNO background
EEC, sigC = np.loadtxt("data/crosssections/XS_nue_12C_NC.txt", unpack=True, delimiter=";")
sigmaC = interp1d(EEC, sigC, fill_value="extrapolate")

def gauss_prof(res, Ee, E_o, offset=0.):
    #deltaE = res*np.sqrt(E_o)  # in MeV
    #deltaE = np.sqrt(res**2.*E_o + offset**2.*E_o**2.)  # in MeV
    #deltaE = res*np.sqrt(E_o) + offset*E_o  # in MeV
    deltaE = -0.123 + 0.376*np.sqrt(E_o) + 0.0349*E_o
    return 1./np.sqrt( 2.*np.pi*deltaE**2. )*np.exp( -1./2.*( Ee-E_o )**2./deltaE**2. )

gauss_prof = np.vectorize(gauss_prof)

def Enu_from_Ee(Ee):
    return 1./2.*(np_dif + m_p - np.sqrt(np_dif**2. - 2.*np_dif*m_p - 4.*Ee*m_p + m_p**2.))

EnuIBDtab, sigIBDtab, EeIBD = np.loadtxt("data/crosssections/XS_IBD.txt", unpack=True)
sigIBDtab *= 1.e-41 # units in cm^2

sigIBD = interp1d(EeIBD, sigIBDtab, fill_value="extrapolate")
EnuIBD = interp1d(EeIBD, EnuIBDtab, fill_value="extrapolate")

if __name__=="__main__":

    """plt.loglog(EeIBD, sigIBDtab, "r-")
    plt.loglog(EeIBD, sigmaIBD(EeIBD), "b-")
    plt.show()
    exit()"""

    """plt.plot(EeIBD, EnuIBDtab, "r-")
    plt.plot(EeIBD, 1./2.*(np_dif + m_p - np.sqrt(np_dif**2. - 2.*np_dif*m_p - 4.*EeIBD*m_p + m_p**2.)), "b--")
    plt.plot(EeIBD, EeIBD + np_dif, "g:")
    plt.show()
    exit()"""

    def fitfunc(E, A, alp):
        return A*E**alp

    def fitfunc2(E, A, alp):
        #return A*np.log10(E) + alp*np.log10(E)**2.
        return A*np.log10(E)**alp

    Eii = np.linspace(EEC[0], EEC[-1])

    #fitpars = np.polyfit(EE, np.log10(sig), 2)
    #fit = np.poly1d(fitpars)
    fitparsAr, firtcov = curve_fit(fitfunc, EEAr, np.log10(sigAr))#, p0=[-4.e1, -0.03])
    fitparsC, firtcov = curve_fit(fitfunc, EEC, np.log10(sigC))#, p0=[-4.e1, -0.03])

    #plt.plot(Eii, 10.**fit(Eii))
    plt.plot(Eii, sigmaAr(Eii), "r-")
    plt.plot(Eii, sigmaC(Eii), "r--")
    plt.plot(Eii, 10.**fitfunc(Eii, *fitparsAr), "b:")
    plt.plot(Eii, 10.**fitfunc(Eii, *fitparsC), "b--")
    #plt.plot(Eii, 10.**(-44.8*Eii**(-0.035)),"g:")
    plt.yscale("log")
    plt.xscale("log")
    print(fitparsAr, fitparsC)
    print(np.abs(10.**fitfunc(Eii, *fitparsAr) - sigmaAr(Eii))/sigmaAr(Eii))
    print(np.abs(10.**fitfunc2(Eii, *fitparsC) - sigmaC(Eii))/sigmaC(Eii))

    #EE, sig = np.loadtxt("data/crosssections/nuAr_XS_Kaes.txt", unpack=True, delimiter=";")
    #sigmaAr = interp1d(EE, sig)
    #plt.plot(EE, sig, "m-", lw=2, alpha=0.7)
    #print(np.abs(sig-sigmaAr(EE))/sig)

    #tab = np.loadtxt("data/crosssections/xs_nue_Ar40.txt", unpack=True)   # from code SNOwGLoBES https://github.com/SNOwGLoBES/snowglobes/tree/master/xscns
    #EE, sig = 10.**tab[0]*1.e3, tab[1]*1e-38*1.e-3
    #sigmaAr = interp1d(EE, sig)
    #plt.plot(EE, sig, "b-", lw=2, alpha=0.7)

    plt.show()