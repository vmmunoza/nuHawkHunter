#-------------------------------
# Compute the event rate at different experiments
#-------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d, Rbf
from Source.flux_stuff import *
from Source.cross_sections import *

#-- Load backgrounds --#

# 1 for computing the number of events with the gaussian profile resolution, 0 assumes a delta function
use_res = 1

# SK
EbackSK, backSK = np.loadtxt("data/backevents/Background_SKII.csv", unpack=True, delimiter=",")
EdatSK, datSK = np.loadtxt("data/backevents/DataPoints_SKII.csv", unpack=True, delimiter=",")
backSK, datSK = backSK/7.82, datSK/7.82     # Exposure 176 kt-years, so 7.82 yrs

# HK (with Gd)
# From 1805.04163, figure 188 right
EbackHK, backHK = np.loadtxt("data/backevents/Back_HK.txt", unpack=True, delimiter=";")
backHK = backHK/2./10.     # Change units
EbackHK-=1. # moving energy bin
bckHK = interp1d(EbackHK, backHK)

EbackSKinvmu, bacSKinvmu = np.loadtxt("data/backevents/Back_SK_invmuon.txt", unpack=True, delimiter=";")
bacSKinvmu = bacSKinvmu*bckHK(EbackSKinvmu[23])/bacSKinvmu[23]
"""plt.plot(EbackSKinvmu, bacSKinvmu,"r-")
plt.plot(EbackHK, backHK,"b:")
plt.show()
exit()"""

# JUNO
# From 1507.05613, figure 5-2
#EbackJUNO, backJUNO = np.loadtxt("data/backevents/Back_JUNO.txt", unpack=True, delimiter=";")
#backJUNO = backJUNO/10.     # Change units
EbackJUNO_NC, backJUNO_NC = np.loadtxt("data/backevents/Back_JUNO_NC.txt", unpack=True, delimiter=";") # Only NC back
backJUNO_NC = backJUNO_NC/10.     # Change units

# DUNE
# From 1804.03157, figure 6
#EbackDUNE, backDUNE = np.loadtxt("data/backevents/Back_DUNE.txt", unpack=True, delimiter=";")
#backDUNE = backDUNE/2./20.     # Change units from 1804.03157

# Atmospheric flux (for DUNE and JUNO), in MeV^-1 s^-1
Eatmos, fluxatmos_e, fluxatmos_antie = np.loadtxt("data/backfluxes/flux_atmos_e.dat", unpack=True)
fluxatmos_e, fluxatmos_antie = fluxatmos_e*year_sec, fluxatmos_antie*year_sec   # to yr^-1


# Function to compute the event rate for an experiment
# E_o: observed energy (of positron, gamma, etc)
# E_nu: neutrino energy array
# flux: neutrino flux evaluated at E_nu
# exp: experiment
def event_rate(E_o, E_nu, flux, exp):

    fluxint = interp1d(E_nu, flux, fill_value="extrapolate")
    #fluxint1 = interp1d(np.log(E_nu), np.log(flux), fill_value="extrapolate")
    #fluxint = lambda E: np.exp(fluxint1(np.log(E)))

    # IBD channel
    if (exp=="SK") or (exp=="HK") or (exp=="JUNO"):

        if exp=="HK":
            ntot = 2.5e34   # 187 kton
            eps = 0.67      # 1804.03157
            # E_visible -> Eo = E_e
            Ee = E_o
            res = 0.5

        elif exp=="SK":
            ntot = 2.5e34*22.5/187  # Same than HK, but with 22.5 kton
            eps = 0.74              # From 1804.03157, without neutron capture efficiency. Check
            # E_visible -> Eo = E_e
            Ee = E_o
            res = 0.5

        elif exp=="JUNO":
            ntot = 1.2e33   # 17 kton
            eps = 0.5       # 1507.05613, for signal, reactor and atm. CC
            # E_visible -> Eo = E_e + m_e
            Ee = E_o - m_e
            res = 0.03


        """
        E_e = E_o
        #return ntot*eps*integrate.simps( flux*dsigmadE_IBD(E_nu, E_e)*E_nu, np.log(E_nu) )
        #return ntot*eps*integrate.simps( flux*dsigmadE_IBD(E_nu, E_e), E_nu )
        #Enuvec = np.linspace(Eth, E_nu[-1], 2000)
        #return ntot*eps*integrate.simps( fluxint(Enuvec)*dsigmadE_IBD(Enuvec, E_e), Enuvec )
        #return ntot*eps*integrate.quad( lambda logEnu: fluxint(np.exp(logEnu))*dsigmadE_IBD(np.exp(logEnu), E_e)*np.exp(logEnu), np.log(Eth), np.log(E_nu[-1]) )[0]
        return ntot*eps*integrate.quad( lambda Enu: fluxint(Enu)*dsigmadE_IBD(Enu, E_e), Eth, E_nu[-1] )[0]
        """


        """if type(E_o) is np.ndarray:
            evntrt = []


            for Eo in E_o:
                #Ee = np.linspace(m_e, Eo*10.)
                print(Eo, fluxint(Ee + np_dif)*sigmaIBD(Ee)*gauss_prof(Ee, Eo))
                evntrt.append( ntot*eps*integrate.quad( lambda Ee: fluxint(Ee + np_dif)*sigmaIBD(Ee)*gauss_prof(Ee, Eo), m_e, Eo*10. )[0] )
                #evntrt.append( ntot*eps*integrate.simps( fluxint(Ee + np_dif)*sigmaIBD(Ee)*gauss_prof(Ee, Eo), Ee ) )
            return np.array(evntrt)

        else:
            #Ee = np.linspace(m_e, E_o*10.)
            #return ntot*eps*integrate.simps( fluxint(Enu_from_Ee(Ee))*sigmaIBD(Ee)*gauss_prof(Ee, E_o), Ee )
            return ntot*eps*integrate.quad( lambda Ee: fluxint(Enu_from_Ee(Ee))*sigmaIBD(Ee)*gauss_prof(Ee, E_o), m_e, E_o*10. )[0]"""

        if use_res:

            Ee = np.linspace(m_e, E_o*10.,1500)
            Enu = EnuIBD(Ee)
            #Enu = Ee + np_dif
            return ntot*eps*integrate.simps( fluxint(Enu)*sigIBD(Ee)*gauss_prof(res, Ee, E_o), Ee )
            #return ntot*eps*integrate.quad( lambda Ee: fluxint(Ee + np_dif)*sigmaIBD(Ee)*gauss_prof(res, Ee, E_o), m_e, E_o*10. )[0]

        else:

            # Get Enu from Ee = (Enu - Delta) (1 - Enu/mp), correction at same order than cross section, see Beacom DSNB review
            #Enu = 1./2.*(np_dif + m_p - np.sqrt(np_dif**2. - 2.*np_dif*m_p - 4.*Ee*m_p + m_p**2.))
            # Enu = Ee + np_dif
            Enu = EnuIBD(Ee)
            return ntot*eps*fluxint(Enu)*sigIBD(Ee)

    # nu-Ar channel
    if exp=="DUNE":
        ntot = 6.02e32  # 40 kton
        eps = 0.86
        res = 0.11  # 2106.05013
        offset = 0.02

        if use_res:

            Enu = np.linspace(5., E_o*10.,500)
            return ntot*eps*integrate.simps( fluxint(Enu)*sigmaAr(Enu)*gauss_prof(res, Enu, E_o, offset), Enu )

        else:

            return ntot*eps*fluxint(E_o)*sigmaAr(E_o)


    # CEnuNS
    if (exp=="DARWIN") or (exp=="ARGO"):
        if exp=="ARGO":
            # Use Argon
            Z = 18.
            A = 40.
            mfrac = 1.  # fraction, used?
            m_uma = 39.948
            m_fid = 360.e6  # g
        if exp=="DARWIN":
            # Use Xenon
            Z = 54.
            A = 132
            mfrac = 1. # take the most abundant
            m_uma = 131.293
            m_fid = 20.e6#40.e6   # g
            #Atope = np.array([128., 129., 130., 131., 132., 134., 136.])
            #mfrac = np.array([0.019, 0.264, 0.041, 0.212, 0.269, 0.104, 0.089])

        eps = 1.
        ntot = m_fid/(m_uma*m_p*1.e-3/gr_to_GeV)*A

        mT = m_p*Z + m_n*(A-Z)
        E_r = E_o
        Enuvec = np.logspace( np.log10(E_nu_min_CE(E_r, mT)), np.log10(E_nu[-1]), 500 )
        return ntot*eps*integrate.simps( sigma_diff_CEnuNS(Enuvec, E_r, A, Z, mT)*fluxint(Enuvec)*np.heaviside(E_r_max(Enuvec, mT)-E_r, 0.), Enuvec)


# Compute the event rates for the backgrounds
# Returns an array of observed energies (somewhat arbitrary in principle) and the background events
def back_rate(exp):

    if exp=="SK":
        return EbackSK - np_dif, backSK

    if exp=="HK":
        #EbackHKnew, backHKnew = EbackHK[2:], backHK[2:] # avoid spallation background above 15 MeV
        #return EbackHKnew, backHKnew
        #return EbackSK, backSK*187/22.5
        #bckHKinvmu = interp1d(EbackSKinvmu, bacSKinvmu, fill_value="extrapolate")
        maxind = 21
        Enuatm, fluxe, fluxebar = Eatmos[2:maxind], fluxatmos_e[2:maxind], fluxatmos_antie[2:maxind]
        backHKCC = np.array([event_rate(Enu - np_dif, Enuatm, fluxebar, exp) for Enu in Enuatm])
        bckHKCC = interp1d(Enuatm, backHKCC, fill_value="extrapolate")
        return EbackSKinvmu - np_dif, bacSKinvmu + bckHKCC(EbackSKinvmu)

    if exp=="JUNO":
        maxind = 21#23
        Enuatm, fluxe, fluxebar = Eatmos[2:maxind], fluxatmos_e[2:maxind], fluxatmos_antie[2:maxind]
        #backCC = event_rate(Enuatm - np_dif + m_e, Enuatm, fluxebar, exp)

        backCC = []
        for Enu in Enuatm:
            backCC.append( event_rate(Enu - np_dif + m_e, Enuatm, fluxebar, exp) )
        backCC = np.array(backCC)

        #ntot_C, eps_NC = 4.505e33, 0.011
        #backNC = 0.#ntot_C*eps_NC*fluxe*sigmaC(Enuatm)
        #print(backNC/backCC)
        backNCfit = interp1d(EbackJUNO_NC, backJUNO_NC, fill_value=(0.,0.), bounds_error=False)
        backNC = backNCfit(Enuatm)
        return Enuatm - np_dif + m_e, backNC + backCC

    if exp=="DUNE":
        # Consider energies above 19 MeV to avoid solar backgrounds. Take up to  ~70 MeV, check
        maxind = 21
        EbackDUNE, fluxatmoslow = Eatmos[5:maxind], fluxatmos_e[5:maxind]
        backDUNE = np.array([event_rate(Eo, EbackDUNE, fluxatmoslow, exp) for Eo in EbackDUNE])
        return EbackDUNE, backDUNE

    if (exp=="DARWIN") or (exp=="ARGO"):
        # For CEnuNS, observed energies of the order of keV
        backfolder = "data/backfluxes/"
        Eatm, atm_nue = np.loadtxt(backfolder+"atmnue_noosc_fluka_flux.dat",unpack=True)
        Eatm, atm_nuebar = np.loadtxt(backfolder+"atmnuebar_noosc_fluka_flux.dat",unpack=True)
        Eatm, atm_numu = np.loadtxt(backfolder+"atmnumu_noosc_fluka_flux.dat",unpack=True)
        Eatm, atm_numubar = np.loadtxt(backfolder+"atmnumubar_noosc_fluka_flux.dat",unpack=True)
        Eatm*=1.e3 # to MeV
        atmflux = (atm_nue + atm_nuebar + atm_numu + atm_numubar)/1.e7 # GeV^-1 m^-2 s^-1 -> MeV^-1 cm^-2 s^-1
        atmint = interp1d(Eatm, atmflux, fill_value=0., bounds_error=False)
        Ehep, sol_hep = np.loadtxt(backfolder+"HEPNeutrinoFlux.dat",unpack=True)
        EB8, sol_B8_1, sol_B8_2, sol_B8_3 = np.loadtxt(backfolder+"B8NeutrinoFlux.dat",unpack=True)
        sol_B8 = sol_B8_1 + sol_B8_2 + sol_B8_3
        hepint = interp1d(Ehep, sol_hep, fill_value=0., bounds_error=False)
        B8int = interp1d(EB8, sol_B8, fill_value=0., bounds_error=False)
        EbackCE = np.logspace(np.log10(EB8[0]), np.log10(Eatm[-1]), 100)
        # Sum backgrounds and correct normalization, see table III of 1812.05550 or Table 2 of 1208.5723
        fluxbacks = atmint(EbackCE) + B8int(EbackCE)*4.59e6 + hepint(EbackCE)*8.31e3
        EobsCE = np.logspace(-4., -1.,200)
        #EobsCE = np.logspace(np.log10(5.e-3), -1.,100)
        backCE = np.array([event_rate(Eo, EbackCE, fluxbacks, exp) for Eo in EobsCE])

        return EobsCE, backCE*year_sec

# Bin energies and events
def binned_events(Eback, events, bin=1.):

    Ebin = np.arange(Eback[0], Eback[-1], step=bin)[:-1]
    #eventsint = interp1d(Eback, events, fill_value="extrapolate")
    eventsint1 = interp1d(np.log(Eback), np.log(events))#, fill_value=0.)
    eventsint = lambda E: np.exp(eventsint1(np.log(E)))
    eventsbin = []
    for Eb in Ebin:
        eventsbin.append( integrate.quad( eventsint, Eb, Eb+bin )[0] )
    return Ebin, np.array(eventsbin)/bin

# Compute the event rate for the signal and backgrounds for a range of PBH masses
def compute_events(Mpbhs, fpbhs, exp, as_DM, plotevents=0, binevents=1):

    Eobs, eventback0 = back_rate(exp)
    #eventback0=eventback0/1.e3*10. # factor for plotting ARGO, CHANGE!
    #eventback0=eventback0/1.e6/20. # factor for plotting ARGO, CHANGE!
    if binevents:
        bin = 1.
        if (exp=="DARWIN") or (exp=="ARGO"):
            bin = 5.e-3
        Ebackbin, eventback = binned_events(Eobs, eventback0, bin)

    for mm, Mpbh in enumerate(Mpbhs):

        fileflux = "fluxes/{:.1e}/flux_isDM_{}.txt".format(Mpbh, as_DM)
        E_nu, flux = np.loadtxt(fileflux, unpack=True)
        # From GeV to MeV, seconds to years
        E_nu, flux = 1.e3*E_nu, 1.e-3*flux*year_sec

        if as_DM:
            fpbhlabel = r" g, $f_{\rm PBH}=$"
        else:
            fpbhlabel = r" g, $\beta'=$"

        events = []

        for E_o in Eobs:
            events.append( event_rate(E_o, E_nu, flux, exp) )
        events = np.array(events)
        #events = 6*events/1.e3*10.  # factor for ARGO, CHANGE!
        #events = 6*events/1.e6/20.  # factor for ARGO, CHANGE!

        if binevents:
            Ebackbin, eventsbin = binned_events(Eobs, events, bin)
            np.savetxt("fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh,exp), np.transpose([Ebackbin, eventsbin]) )
        else:
            np.savetxt("fluxes/{:.1e}/event_rate_{}.txt".format(Mpbh,exp), np.transpose([Eobs, events]) )

        """
        print("\nMpbh",np.log10(Mpbh))
        print("E",Eback[0], Eback[-1],"step",Eback[-2]-Eback[0])
        factor = 10*20./17.  # 10 years, mass correction
        bin = Eback[-2]-Eback[0]
        bin = 40-12
        facfpbh = 1.3e-1  #5.5e-4
        print("back",binned_events(Eback, eventback0, bin)[1]*factor)
        print("events",binned_events(Eback, events, bin)[1]*factor*6)   # 6 dof, fpbh juno paper
        """

        if plotevents:
            if binevents:
                plt.step(Ebackbin, fpbhs[mm]*eventsbin, color=cols[mm], label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]))
            else:
                plt.plot(Eobs, fpbhs[mm]*events, color=cols[mm], label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]))

    if plotevents:
        if binevents:
            plt.xlim(Ebackbin[0], Ebackbin[-1])
            plt.fill_between(Ebackbin, np.zeros_like(Ebackbin), eventback, color="b", alpha=0.3, label="Background",step="pre")
        else:
            plt.xlim(Eobs[0], Eobs[-1])
            plt.fill_between(Eobs, np.zeros_like(Eobs), eventback0, color="b", alpha=0.3, label="Background")
        if exp=="SK":
            plt.scatter(EdatSK, datSK, color="b", marker="o", label = r"SK-II Data")
