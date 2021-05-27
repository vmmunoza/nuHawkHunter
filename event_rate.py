
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d, Rbf
from flux_stuff import *
from cross_sections import *

# TAKE ONLY electron, now it is taking all species / 6

# Backgrounds

# SK
EbackSK, backSK = np.loadtxt("data/Background_SKII.csv", unpack=True, delimiter=",")
EdatSK, datSK = np.loadtxt("data/DataPoints_SKII.csv", unpack=True, delimiter=",")
backSK, datSK = backSK/7.82, datSK/7.82     # Exposure 176 kt-years, so 7.82 yrs

# HK (with Gd)
# From 1805.04163, figure 188 right
EbackHK, backHK = np.loadtxt("data/Back_HK.txt", unpack=True, delimiter=";")
backHK = backHK/2./10.     # Change units

# JUNO
# From 1507.05613, figure 5-2
EbackJUNO, backJUNO = np.loadtxt("data/Back_JUNO.txt", unpack=True, delimiter=";")
backJUNO = backJUNO/10.     # Change units  (CHECK THIS NORMALIZATION)

# DUNE
# From 1804.03157, figure 6
#EbackDUNE, backDUNE = np.loadtxt("data/Back_DUNE.txt", unpack=True, delimiter=";")
#backDUNE = backDUNE/2./20.     # Change units from 1804.03157

# Atmospheric flux (for DUNE and JUNO), in MeV^-1 s^-1 (check)
Eatmos, fluxatmos_e, fluxatmos_antie = np.loadtxt("data/flux_atmos_e.dat", unpack=True)
fluxatmos_e, fluxatmos_antie = fluxatmos_e*year_sec, fluxatmos_antie*year_sec   # to yr^-1


# From 1804.03157
def event_rate(E_o, E_nu, flux, exp):

    fluxint = interp1d(E_nu, flux, fill_value="extrapolate")

    if (exp=="SK") or (exp=="HK"):

        if exp=="HK":
            ntot = 2.5e34   # 187 kton
            eps = 0.67      # 1804.03157
        if exp=="SK":
            ntot = 2.5e34*22.5/187  # Same than HK, 22.5 kton
            eps = 0.67              # Check

        """
        E_e = E_o
        #return ntot*eps*integrate.simps( flux*dsigmadE_IBD(E_nu, E_e)*E_nu, np.log(E_nu) )
        #return ntot*eps*integrate.simps( flux*dsigmadE_IBD(E_nu, E_e), E_nu )
        #Enuvec = np.linspace(Eth, E_nu[-1], 2000)
        #return ntot*eps*integrate.simps( fluxint(Enuvec)*dsigmadE_IBD(Enuvec, E_e), Enuvec )
        #return ntot*eps*integrate.quad( lambda logEnu: fluxint(np.exp(logEnu))*dsigmadE_IBD(np.exp(logEnu), E_e)*np.exp(logEnu), np.log(Eth), np.log(E_nu[-1]) )[0]
        return ntot*eps*integrate.quad( lambda Enu: fluxint(Enu)*dsigmadE_IBD(Enu, E_e), Eth, E_nu[-1] )[0]
        """
        return ntot*eps*fluxint(E_o + np_dif)*sigmaIBD(E_o)

    if exp=="JUNO":
        ntot = 1.2e33   # 17 kton
        eps = 0.5       # 1507.05613, for signal, reactor and atm. CC
        # E_visible -> Eo = E_e + m_e
        return ntot*eps*fluxint(E_o + np_dif - m_e)*sigmaIBD(E_o - m_e)

    if exp=="DUNE":

        ntot = 6.02e32  # 40 kton
        eps = 0.86
        return ntot*eps*fluxint(E_o)*sigmaAr(E_o)

def back_rate(exp):

    if exp=="SK":
        return EbackSK, backSK

    if exp=="HK":
        return EbackHK, backHK
        #return EbackSK, backSK*187/22.5

    if exp=="JUNO":
        maxind = 23
        EbackJUNO, fluxe, fluxebar = Eatmos[2:maxind], fluxatmos_e[2:maxind], fluxatmos_antie[2:maxind]
        backCC = event_rate(EbackJUNO, EbackJUNO, fluxebar, exp)
        ntot_C, eps_NC = 4.505e33, 0.011
        backNC = 0.#ntot_C*eps_NC*fluxe*sigmaC(EbackJUNO)
        #print(backNC/backCC)
        return EbackJUNO, backNC + backCC

    if exp=="DUNE":
        # Consider energies above 19 MeV to avoid solar backgrounds. Take up to  ~70 MeV, check
        EbackDUNE, fluxatmoslow = Eatmos[5:17], fluxatmos_e[5:17]
        backDUNE = np.array([event_rate(Eo, EbackDUNE, fluxatmoslow, exp) for Eo in EbackDUNE])
        return EbackDUNE, backDUNE

def binned_events(Eback, events, bin=1.):
    # Binned events
    #bin = 1.
    Ebin = np.arange(Eback[0], Eback[-1], step=bin)[:-1]
    eventsint = interp1d(Eback, events, fill_value="extrapolate")
    eventsbin = []
    for Eb in Ebin:
        eventsbin.append( integrate.quad( eventsint, Eb, Eb+bin )[0] )
    return Ebin, np.array(eventsbin)

def compute_events(Mpbhs, fpbhs, exp, plotevents=0):

    Eback, eventback0 = back_rate(exp)
    Ebackbin, eventback = binned_events(Eback, eventback0)

    for mm, Mpbh in enumerate(Mpbhs):

        folder = "folder_fluxes/{:.1e}/flux.txt".format(Mpbh)
        E_nu, flux = np.loadtxt(folder, unpack=True)
        # From GeV to MeV, seconds to years
        E_nu, flux = 1.e3*E_nu, 1.e-3*flux/6.*year_sec    # El 1/6 es para coger solo electron, hacerlo bien!

        """plt.loglog(Eatmos, fluxatmos_e,"b")
        plt.loglog(E_nu, flux,"r")
        plt.xlim(1,100)
        plt.ylim(1e-4,1e10)
        plt.show()
        exit()"""

        if Mpbh<Mevap:
            fpbhlabel = r" g, $\beta'=$"
        else:
            fpbhlabel = r" g, $f_{\rm PBH}=$"

        events = []

        #Evec = np.linspace(Eback[0], Eback[-1], 200)
        #for E_o in Evec:
        for E_o in Eback:
            events.append( event_rate(E_o, E_nu, flux, exp) )
        events = np.array(events)

        Ebackbin, eventsbin = binned_events(Eback, events)
        #Ebackbin, eventsbin = binned_events(Evec, events)

        np.savetxt("folder_fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh,exp), np.transpose([Ebackbin, eventsbin]) )

        """print("\nMpbh",np.log10(Mpbh))
        print("E",Eback[0], Eback[-1])
        print("back",binned_events(Eback, eventback0, Eback[-2]-Eback[0])[1]*10)
        print("events",binned_events(Eback, events, Eback[-2]-Eback[0])[1]*10*5.5e-4)"""

        if plotevents:
            plt.plot(Ebackbin, fpbhs[mm]*eventsbin, color=cols[mm], label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]))

        #return Eback, events

    if plotevents:
        plt.xlim(Ebackbin[0], Ebackbin[-1])
        plt.fill_between(Ebackbin, np.zeros_like(Ebackbin), eventback, "grey", alpha=0.5, label="Background")
        if exp=="SK":
            plt.scatter(EdatSK, datSK, color="b", marker="o", label = r"SK-II Data")


if __name__=="__main__":

    plotevents = 1

    #Mpbhs =  [1.e12, 1.e13, 1.e14]
    Mpbhs =  [1e15]#, 2e15, 4.e15]
    Mpbhs =  [1e15, 2e15, 4.e15]
    #Mpbhs = np.linspace(1.e15, 1.e16,10)

    #fpbhs = np.ones_like(Mpbhs)    # this is beta prime or fpbh
    fpbhs = 1.e-2*np.ones_like(Mpbhs)
    cols = ["r", "purple", "b", "g", "orange"]


    #exp = "SK"
    exp = "HK"
    #exp = "JUNO"
    #exp = "DUNE"

    compute_events(Mpbhs, fpbhs, exp, plotevents)

    #plt.ylim(1.e-4, 3.e1)
    plt.yscale("log")
    plt.legend()
    plt.xlabel('$E{\\rm \,\, [MeV]}$')
    plt.ylabel('${\\rm d}N/d E \,\, [{\\rm MeV}^{-1}{\\rm yr}^{-1}]$')
    plt.title(exp)
    plt.savefig("figures/events_"+exp+".png", bbox_inches='tight')
    plt.show()
    plt.gcf().clear()
