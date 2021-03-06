#-------------------------------
# Compute the event rate at different experiments
#-------------------------------


import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d, Rbf
from Source.fluxes import *
from Source.cross_sections import *

#---------------------
# Load backgrounds from external files
#---------------------

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

# JUNO
# From 1507.05613, figure 5-2
EbackJUNO_NC, backJUNO_NC = np.loadtxt("data/backevents/Back_JUNO_NC.txt", unpack=True, delimiter=";") # Only NC back
backJUNO_NC = backJUNO_NC/10.     # Change units

# DUNE
# From 1804.03157, figure 6
#EbackDUNE, backDUNE = np.loadtxt("data/backevents/Back_DUNE.txt", unpack=True, delimiter=";")
#backDUNE = backDUNE/2./20.     # Change units from 1804.03157

# Atmospheric flux (for DUNE and JUNO), in MeV^-1 s^-1
Eatmos, fluxatmos_e, fluxatmos_antie = np.loadtxt("data/backfluxes/flux_atmos_e.dat", unpack=True)
fluxatmos_e, fluxatmos_antie = fluxatmos_e*year_sec, fluxatmos_antie*year_sec   # to yr^-1

#---------------------
# Experiment class
#---------------------

# name: name of the experiment
# typeexp: experiment type, possible options are: WaterCherenkov, WaterCherenkovGd, LiXeCEnuNS, LiquidScintillator, LiquidArgon, LiArCEnuNS
# ntot: total number of targets, or fiducial mass if use_mfid==1
# eps: detection efficiency of the experiment
# res: energy resolution parameters (see energy_resolution function definition)
# lat: latitude of the experiment (relevant for atmospheric background)
# bin: energy bin size
# use_res: 1 for computing the number of events with the gaussian profile resolution, 0 assumes a delta function for the resolution (and the res is not employed)
# use_mfid: if use_mfid==1, compute ntot from the fiducial mass
class experiment():
    def __init__(self, name, typeexp, ntot, eps, res, bin, lat, use_res=1, use_mfid=0):

        self.name = name
        self.type = typeexp
        self.ntot = ntot
        self.eps = eps
        self.res = res
        self.bin = bin
        self.lat = lat
        self.use_res = use_res
        self.use_mfid = use_mfid

        if self.type=="LiArCEnuNS":
            # Liquid Argon
            self.Z = 18.
            self.A = 40.
            self.m_uma = 39.948
            self.mT = m_p*self.Z + m_n*(self.A-self.Z)
        if self.type=="LiXeCEnuNS":
            # Liquid Xenon
            self.Z = 54.
            self.A = 132
            self.m_uma = 131.293
            self.mT = m_p*self.Z + m_n*(self.A-self.Z)

        # If use_mfid==1, compute ntot from the fiducial mass
        if use_mfid==1:
            self.ntot = self.ntot*(gr_to_GeV*1.e3)/(self.m_uma*m_p)

        # Define the channel given the experiment type, among IBD, nuAr and CEnuNS
        if self.type=="WaterCherenkov" or self.type=="WaterCherenkovGd" or self.type=="LiquidScintillator":
            self.channel = "IBD"
        elif self.type=="LiquidArgon":
            self.channel = "nuAr"
        elif self.type=="LiXeCEnuNS" or self.type=="LiArCEnuNS":
            self.channel = "CEnuNS"
        else:
            print("Experiment type not valid.\n")

        # Compute background event rate
        self.Eback, self.backrate = self.back_rate()
        # Compute binned background event rate
        self.Eback_bin, self.backrate_bin = self.binned_events(self.Eback, self.backrate)

        # Relevant energy range
        self.Emin, self.Emax = self.Eback[0], self.Eback[-1]

    # Latitude correction for atmospheric background, interpolated from the values shown in arXiv:astro-ph/0701305, Table III
    def lat_factor(self):
        lats = [1.5, 36.5, 40., 45., 63.7]
        s_lats = [0.8, 1., 1.25, 1.5, 2.]
        s_lat = interp1d(lats, s_lats)
        return s_lat(self.lat)

    # Energy resolution delta_E (in MeV, E_o in MeV)
    def energy_resolution(self, E):
        Eterm, sqrtterm, offset = self.res
        deltaE = Eterm*E + sqrtterm*np.sqrt(E) + offset
        #deltaE = np.sqrt(self.res**2.*E + self.offset**2.*E**2.)  # in MeV
        return deltaE

    # Gaussian resolution profile
    def gauss_prof(self, Etrue, E_o):
        deltaE = self.energy_resolution(Etrue)
        return 1./np.sqrt( 2.*np.pi*deltaE**2. )*np.exp( -1./2.*( Etrue-E_o )**2./deltaE**2. )

    # Background rate
    def back_rate(self):

        if self.type=="WaterCherenkov":
            return EbackSK - np_dif, backSK

        elif self.type=="WaterCherenkovGd":
            lat = self.lat_factor() # correction for latitude
            minind, maxind = 2, 21
            Enuatm, fluxe, fluxebar = Eatmos[minind:maxind], fluxatmos_e[minind:maxind], fluxatmos_antie[minind:maxind]
            backHKCC = np.array([self.event_rate(Enu - np_dif, Enuatm, fluxebar) for Enu in Enuatm])*lat
            bckHKCC = interp1d(Enuatm, backHKCC, fill_value="extrapolate")
            return EbackSKinvmu - np_dif, bacSKinvmu + bckHKCC(EbackSKinvmu)

        elif self.type=="LiquidScintillator":
            lat = self.lat_factor() # correction for latitude
            minind, maxind = 2, 21
            Enuatm, fluxe, fluxebar = Eatmos[minind:maxind], fluxatmos_e[minind:maxind]*lat, fluxatmos_antie[minind:maxind]*lat
            backCC = np.array( [self.event_rate(Enu - np_dif + m_e, Enuatm, fluxebar) for Enu in Enuatm ] )

            # NC contribution, negligible
            #ntot_C, eps_NC = 4.505e33, 0.011
            #backNC = 0.#ntot_C*eps_NC*fluxe*sigmaC(Enuatm)

            backNCfit = interp1d(EbackJUNO_NC, backJUNO_NC, fill_value=(0.,0.), bounds_error=False)
            backNC = backNCfit(Enuatm)
            return Enuatm - np_dif + m_e, backNC + backCC

        elif self.type=="LiquidArgon":
            lat = self.lat_factor() # correction for latitude
            # Consider energies above 19 MeV to avoid solar backgrounds. Take up to  ~70 MeV, check
            # Include lower energies for computing the solar background events, since these are relevant due to energy resolution, later take only energy bins above 19 MeV
            minind, maxind = 5, 21
            EbackDUNE, fluxatmoslow = Eatmos[:maxind], fluxatmos_e[:maxind]*lat
            backDUNE = np.array([self.event_rate(Eo, EbackDUNE, fluxatmoslow) for Eo in EbackDUNE])
            # Include solar background
            backfolder = "data/backfluxes/"
            Ehep, sol_hep = np.loadtxt(backfolder+"HEPNeutrinoFlux.dat",unpack=True)
            EB8, sol_B8_1, sol_B8_2, sol_B8_3 = np.loadtxt(backfolder+"B8NeutrinoFlux.dat",unpack=True)
            sol_B8 = sol_B8_1 #+ sol_B8_2 + sol_B8_3
            hepint = interp1d(Ehep, sol_hep, fill_value=0., bounds_error=False)
            B8int = interp1d(EB8, sol_B8, fill_value=0., bounds_error=False)
            # Sum backgrounds and correct normalization, see table III of 1812.05550 or Table 2 of 1208.5723
            backsolarDUNE = ( B8int(EbackDUNE)*4.59e6 + hepint(EbackDUNE)*8.31e3 )*year_sec
            backDUNE += np.array([self.event_rate(Eo, EbackDUNE, backsolarDUNE) for Eo in EbackDUNE])
            return EbackDUNE[minind:], backDUNE[minind:]

        elif self.channel=="CEnuNS":
            lat = self.lat_factor() # correction for latitude
            # For CEnuNS, observed energies of the order of keV
            backfolder = "data/backfluxes/"
            Eatm, atm_nue = np.loadtxt(backfolder+"atmnue_noosc_fluka_flux.dat",unpack=True)
            Eatm, atm_nuebar = np.loadtxt(backfolder+"atmnuebar_noosc_fluka_flux.dat",unpack=True)
            Eatm, atm_numu = np.loadtxt(backfolder+"atmnumu_noosc_fluka_flux.dat",unpack=True)
            Eatm, atm_numubar = np.loadtxt(backfolder+"atmnumubar_noosc_fluka_flux.dat",unpack=True)
            Eatm*=1.e3 # to MeV
            atmflux = (atm_nue + atm_nuebar + atm_numu + atm_numubar)*lat/1.e7 # latitude correction, GeV^-1 m^-2 s^-1 -> MeV^-1 cm^-2 s^-1
            atmint = interp1d(Eatm, atmflux, fill_value=0., bounds_error=False)
            Ehep, sol_hep = np.loadtxt(backfolder+"HEPNeutrinoFlux.dat",unpack=True)
            EB8, sol_B8_1, sol_B8_2, sol_B8_3 = np.loadtxt(backfolder+"B8NeutrinoFlux.dat",unpack=True)
            sol_B8 = sol_B8_1 #+ sol_B8_2 + sol_B8_3
            hepint = interp1d(Ehep, sol_hep, fill_value=0., bounds_error=False)
            B8int = interp1d(EB8, sol_B8, fill_value=0., bounds_error=False)
            EbackCE = np.logspace(np.log10(EB8[0]), np.log10(Eatm[-1]), 100)
            # Sum backgrounds and correct normalization, see table III of 1812.05550 or Table 2 of 1208.5723
            fluxbacks = atmint(EbackCE) + B8int(EbackCE)*4.59e6 + hepint(EbackCE)*8.31e3
            #EobsCE = np.logspace(-4., -1.,200)
            EobsCE = np.logspace(np.log10(5.e-3), -1.,100)
            backCE = np.array([self.event_rate(Eo, EbackCE, fluxbacks) for Eo in EobsCE])
            backCE = backCE/6.  # I multiply by 6 in self.event_rate for the 6 dof, so here we must divide it (this is very ad hoc)
            return EobsCE, backCE*year_sec

        else:
            print("Experiment not valid: background not found.\n")

    # "Raw" event rate, without resolution integral
    # E_o: observed energy (of positron, photon, etc)
    # E_nu: neutrino energy array
    # flux: neutrino flux evaluated at E_nu
    def nonint_event_rate(self, E_o, E_nu, flux):

        infidxs = np.argwhere(~np.isinf(flux)).reshape(-1)
        E_nu, flux = E_nu[infidxs], flux[infidxs]

        fluxint = interp1d(E_nu, flux, fill_value="extrapolate")

        # Inverse Beta Decay channel (IBD)
        if self.channel == "IBD":

            # Get Enu from Ee = (Enu - Delta) (1 - Enu/mp), correction at same order than cross section, see Beacom DSNB review
            #Enu = 1./2.*(np_dif + m_p - np.sqrt(np_dif**2. - 2.*np_dif*m_p - 4.*Ee*m_p + m_p**2.))
            # Enu = Ee + np_dif

            # In Liquid Scintillator, E_visible -> Eo = E_e + m_e
            if self.type == "LiquidScintillator":
                E_e = E_o - m_e
            else:
                E_e = E_o
            Enu = EnuIBD(E_e)

            events = self.ntot*self.eps*fluxint(Enu)*sigIBD(E_e)

        # nuArgon channel
        elif self.channel == "nuAr":

            events = self.ntot*self.eps*fluxint(E_o)*sigmaAr(E_o)

        # CEnuNS channel
        elif self.channel == "CEnuNS":

            E_r = E_o
            Enuvec = np.logspace( np.log10(E_nu_min_CE(E_r, self.mT)), np.log10(E_nu[-1]), 500 )
            events = self.ntot*self.eps*integrate.simps( sigma_diff_CEnuNS(Enuvec, E_r, self.A, self.Z, self.mT)*fluxint(Enuvec)*np.heaviside(E_r_max(Enuvec, self.mT)-E_r, 0.), Enuvec)
            # CEnuNS is sensitive to 6 degrees of freedom, since only primary spectrum is relevant, just multiply by 6
            events = 6*events

        return events


    # Method to compute the event rate for an experiment
    # If use_res==1, compute resolution integral, otherwise call nonint_event_rate
    # E_o: observed energy (of positron, photon, etc)
    # E_nu: neutrino energy array
    # flux: neutrino flux evaluated at E_nu
    def event_rate(self, E_o, E_nu, flux):

        if self.use_res:
            # IBD channel
            if self.channel == "IBD":
                Ee = np.linspace(m_e, E_o*10.,1500)
                Enu = EnuIBD(Ee)
                events = integrate.simps( self.nonint_event_rate(Ee, E_nu, flux)*self.gauss_prof(Ee, E_o), Ee )

            # nuAr channel
            elif self.channel == "nuAr":
                Enu = np.linspace(5., E_o*10.,500)
                events = integrate.simps( self.nonint_event_rate(Enu, E_nu, flux)*self.gauss_prof(Enu, E_o), Enu )

            # CEnuNS channel
            elif self.channel == "CEnuNS":
                E_r = np.linspace(0.001,E_o*10,500)
                #events = integrate.simps( self.nonint_event_rate(E_r, E_nu, flux)*self.gauss_prof(E_r, E_o), E_r )
                events = integrate.quad(lambda E_r: self.nonint_event_rate(E_r, E_nu, flux)*self.gauss_prof(E_r, E_o), 0.001, E_o*10, epsrel=1e-4, limit=300 )[0]

            return events
        else:
            return self.nonint_event_rate(E_o, E_nu, flux)


    # Binned events
    def binned_events(self, Evec, events):

        Ebin = np.arange(Evec[0], Evec[-1], step=self.bin)[:-1]
        #eventsint = interp1d(Eback, events, fill_value="extrapolate")
        eventsint1 = interp1d(np.log(Evec), np.log(events))#, fill_value=0.)
        eventsint = lambda E: np.exp(eventsint1(np.log(E)))
        eventsbin = []
        for Eb in Ebin:
            eventsbin.append( integrate.quad( eventsint, Eb, Eb+self.bin, limit=400, epsrel= 1e-5 )[0] )
        return Ebin, np.array(eventsbin)#/self.bin

# Compute the event rate for the signal and backgrounds for a range of PBH masses
def compute_events(Mpbhs, fpbhs, exp, as_DM, mass_spec = 0, sig = 0, plotevents=0, binevents=1):

    # Sufix for outputs depending on the mass function
    sufx = sufix(mass_spec, sig)

    for mm, Mpbh in enumerate(Mpbhs):

        fileflux = "fluxes/{:.1e}/flux_isDM_{}".format(Mpbh, as_DM)+sufx
        E_nu, flux = np.loadtxt(fileflux, unpack=True)
        # From seconds to years
        flux = flux*year_sec

        if as_DM:
            fpbhlabel = r" g, $f_{\rm PBH}=$"
        else:
            fpbhlabel = r" g, $\beta'=$"

        # Define an array of observed energies, taken for convenience as the background energy vector
        Eobs = exp.Eback

        # Compute event rate for PBH signal
        events = np.array( [exp.event_rate(E_o, E_nu, flux) for E_o in Eobs] )

        if binevents:
            Ebackbin, eventsbin = exp.binned_events(Eobs, events)
            np.savetxt("fluxes/{:.1e}/event_rate_{}".format(Mpbh,exp.name)+sufx, np.transpose([Ebackbin, eventsbin]) )
        else:
            np.savetxt("fluxes/{:.1e}/event_rate_{}".format(Mpbh,exp.name)+sufx, np.transpose([Eobs, events]) )

        if plotevents:
            if binevents:
                plt.step(Ebackbin, fpbhs[mm]*eventsbin, color=cols[mm], label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]))
            else:
                plt.plot(Eobs, fpbhs[mm]*events, color=cols[mm], label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]))

    if plotevents:
        if binevents:
            plt.xlim(Ebackbin[0], Ebackbin[-1])
            plt.fill_between(exp.Eback_bin, np.zeros_like(exp.Eback_bin), exp.backrate_bin, color="grey", alpha=0.3, label="Background",step="pre")
        else:
            plt.xlim(Eobs[0], Eobs[-1])
            plt.fill_between(exp.Eback, np.zeros_like(exp.Eback), exp.backrate, color="grey", alpha=0.3, label="Background")
        if exp.name=="SK":
            plt.scatter(EdatSK, datSK, color="b", marker="o", label = r"SK-II Data")

#---
# Experiments
#---

# Defintion of the different future experiments
# See the paper for further information regarding the employed values of the parameters
# The experiment class allow to define customized experiments using some of the available channels but with different settings

"""
# Super Kamiokande
SK = experiment(name = "SK",
                typeexp = "WaterCherenkov",
                ntot = 1.5e33,   # 22.5 kton
                eps = 0.74,      # From 1804.03157, without neutron capture efficiency. Check
                res = [0.0349, 0.376, -0.123], #[0., 0.5, 0.],
                bin = 1.,
                lat = 36.5)
"""

# Hyper Kamiokande
HK = experiment(name = "HK",
                typeexp = "WaterCherenkovGd",
                ntot = 2.5e34,   # 187 kton
                eps = 0.67,      # 1804.03157
                res = [0.0349, 0.376, -0.123],  # 1010.0118
                bin = 2.,
                lat = 36.5)

# JUNO
JUNO = experiment(name = "JUNO",
                  typeexp = "LiquidScintillator",
                  ntot = 1.2e33,   # 17 kton
                  eps = 0.5,       # 1507.05613, for signal, reactor and atm. CC
                  res = [0., 0.03, 0.], # 1507.05613
                  bin = 4.,
                  lat = 22.22)

# DUNE
DUNE = experiment(name = "DUNE",
                  typeexp = "LiquidArgon",
                  ntot = 6.02e32,   # 40 kton
                  eps = 0.86,       # 1608.07853, 1804.03157
                  res = [0.2, 0., 0.], # 2002.03005
                  bin = 4.,
                  lat = 44.25)

# DARWIN
DARWIN = experiment(name = "DARWIN",
                    typeexp = "LiXeCEnuNS",
                    ntot = 40.e6,   # fiducial mass in g, 0.04 kton (fiducial mass rather than ntot)
                    eps = 1.,
                    res = [0.06, 0.32/10.**(3./2.), 0.], # 10.**(3./2.) factor for keV -> MeV conversion
                    lat = 42.25,    # Gran Sasso latitude
                    bin = 5.e-3,    # 5keV resolution
                    use_mfid = 1)   # use_mfid = 1: ntot value is fiducial mass, not number of targets

# ARGO
ARGO = experiment(name = "ARGO",
                  typeexp = "LiArCEnuNS",
                  ntot = 370.e6,   # fiducial mass in g, 0.04 kton (rather than ntot)
                  eps = 1.,
                  res = [0.06, 0.32/10.**(3./2.), 0.], # 10.**(3./2.) factor for keV -> MeV conversion
                  lat = 42.25,    # Gran Sasso latitude
                  bin = 5.e-3,    # 5keV resolution
                  use_mfid = 1)   # use_mfid = 1: ntot value is fiducial mass, not number of targets
