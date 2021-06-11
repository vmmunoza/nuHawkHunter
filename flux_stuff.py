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

for path in ["figures", "folder_fluxes"]:
    if not os.path.exists(path):
        os.system("mkdir "+path)

#---
# Miscelanea
#---

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

# cosmo
Om_rad = 9.2e-5
Om_m = 0.315
Om_dm = 0.265
Om_L = 0.685
year_sec = 3.154e7 # 1yr in seconds
hlittle = 0.68
H0 = hlittle*3.2407e-18 # s^-1
G = 6.67259e-8  # cm^3 g^-1 s^-2
rho_c = 3.*H0**2./(8.*np.pi*G) # g cm^-3

Msun = 1.989e33 # g
c = 29979245800.0 #cm/s
MpcToCm = 3.086e24 # 1 Mpc = MpcToCm cm
GpToCm = 1.e3*MpcToCm
gr_to_GeV= 5.62e23

# Black-body approximation
mass_conversion = 5.60958884e+23 # g to GeV
time_conversion = 1.519267407e+24 # s to GeV^(-1)
leng_conversion = 5.06773058e+13 # cm to GeV^(-1)
G_N = (6.67408e-11*pow(leng_conversion*100.,3.)/(mass_conversion*1000.)/pow(time_conversion,2.))

cols = ["r", "m", "purple", "b", "g", "orange"]

# Hawking temperature, in GeV, Mpbh in g
def Tpbh(Mpbh):
    return 1.06e13/Mpbh

def gamma(E, Mpbh):
    if E>4.*Tpbh(Mpbh):
        return E**2.*27.*np.pi*G_N**2.*Mpbh**2.
    else:
        return E**2.*2.*G_N**2.*Mpbh**2.

def blackbody(E, Mpbh):
    if E>Tpbh(Mpbh)*100.:
        fdistr =  np.exp( -E/Tpbh(Mpbh) )   # Avoid numerical problems
    else:
        fdistr = ( np.exp( E/Tpbh(Mpbh) ) +1. )**(-1.)
    dNdtdE = gamma(E, Mpbh)*fdistr/(2.*np.pi)
    return 6.*dNdtdE*mass_conversion**2.*time_conversion    # 6 degrees of freedom, 3 flavours and particle-antiparticle (as in Blackhawk data, I think)

def dNdEdt_extension(Elim,interp,E,Mpbh):
    if E>Elim:
        return blackbody(E, Mpbh)
    else:
        return interp(E)

gamma = np.vectorize(gamma)
blackbody = np.vectorize(blackbody)
dNdEdt_extension = np.vectorize(dNdEdt_extension)

Mevap = 4.e14
# n_pbh at z = 0
# beta and f_pbh defined outside of the function for the moment (i.e., f_PBH=1 or beta'=1 for this formula, scale it later)
def n_pbh(Mpbh):
    if Mpbh<Mevap:
        # See Carr 2010
        return (GpToCm)**(-3.)/(7.98e-29)*(Mpbh/Msun)**(-3./2.)
    else:
        return Om_dm*rho_c/Mpbh




# Approx lifetime, in seconds (Carr et al 2010)
def tau_pbh(M):
    return 407.*(M/1e10)**3.

zeq = 2.4e4*Om_m*hlittle**2. -1.

# Approx lifetime, in seconds (Carr et al 2010)
def tau_pbh(M):
    return 407.*(M/1.e10)**3.

zeq = 2.4e4*Om_m*hlittle**2. -1.

# COSMOLOGY

def hubble(z):      # s^{-1}
    return H0*np.sqrt( Om_rad*(1.+z)**(4.) + Om_m*(1.+z)**(3.) + Om_L )

def dtdz(z):
    return (hubble(z)*(1.+z))**(-1.)

# z as a function of t
def redshift(t):
    # Lookback time (Mo et al book)
    #t = 2./3./H0*Om_m**(-1./2.)*(1.+z)**(-3./2.)
    #z = (3.*H0*t/2.*Om_m**(1./2.))**(-2./3.) - 1.
    z = ( t/( 2./3./H0*Om_m**(-1./2.) ) )**(-2./3.) - 1.
    if z > zeq:
        z = 0.5*1.e10*t**(-1./2.) - 1.
    return z

def time_from_z(z):
    #t = 2./3./H0*Om_m**(-1./2.)*(1.+z)**(-3./2.)
    t = 2./3./H0*Om_L**(-1./2.)*np.log( (np.sqrt( Om_L/(1.+z)**3. ) + np.sqrt( Om_m + Om_L/(1.+z)**3. ))/np.sqrt( Om_m) )
    if z > zeq:
        t = ((1.+z)/1.e10)**(-2.)
    return t

def time_int_z(zmin, zmax):
    return integrate.quad(lambda a: 1./(a*hubble(a**(-1) - 1.)), 1./(1.+zmax), 1./(1.+zmin))[0] # s
    #return integrate.quad(lambda zz: 1./((1.+zz)*hubble(zz)), zmin, zmax)[0] # s
    #return integrate.quad(lambda zz: dtdz(zz), zmin, zmax)[0] # s

def time_int_a(amin, amax):
    #return integrate.quad(lambda a: dtdz(a**(-1.) - 1.), amin, amax)[0] # s
    return integrate.quad(lambda a: 1./(a*hubble(a**(-1) - 1.)), amin, amax)[0] # s

def z_from_t(t):
    return 1./(fsolve(lambda aa: time_int_a(0., aa)-t, 1e-3))[0] - 1.

time_int_a = np.vectorize(time_int_a)
z_from_t = np.vectorize(z_from_t)
redshift = np.vectorize(redshift)

sectoyear = 3.17098e-8

# age of the universe
ageuniverse = time_int_a(0., 1.)


tvec = np.logspace(-30, np.log10(time_int_a(0., 1.)), 500)
z_from_t_int = interp1d(tvec, z_from_t(tvec))#, fill_value="extrapolate")

"""plt.loglog(tvec, z_from_t_int(tvec),"r-")
plt.loglog(tvec, z_from_t(tvec),"b:")
plt.show()
exit()"""

"""
timess = np.logspace(-30, 17)
plt.loglog( timess, z_from_t(timess), "r-" )
plt.loglog( timess, redshift(timess), "b:" )
plt.show()
exit()
"""
#print(z_from_t(1e9/sectoyear), z_from_t(time_int_a(0., 1.)), time_int_a(0., 1.)*sectoyear/1e9)
#exit()

"""
for z in [1e10, 1000, 100]:
    #print(time_int_z(0.,z), time_from_z(z))
    tt = time_int_a(0., 1./(1.+z))
    print( z, z_from_t(tt) )





#zvec = np.logspace(10,0)
#plt.loglog(zvec, time_int(0.,zvec)*sectoyear)
#plt.show()

tvec = np.logspace(-30, 15)
plt.loglog(tvec, z_from_t(tvec))
plt.ylabel("z")
plt.xlabel("t [s]")
plt.show()

exit()
"""


redshift = np.vectorize(redshift)

# Evaporation redshift, valid for the matter dominated epoch at z>>1 or for EdS
def zevap(M):
    return redshift(tau_pbh(M))

# Compute the diffuse flux
def flux(zmin, zmax, Mpbh, E_vec, spec_int, aprox=0):
    flux_vec = []
    logonepluszz = np.linspace( np.log(1.+zmin), np.log(1.+zmax) , 100)
    onepluszz = np.exp(logonepluszz)
    for E in E_vec:
        if aprox:
            rate = blackbody(E*onepluszz, Mpbh)
        else:
            rate = dNdEdt_extension(E_vec[-1],spec_int,E*onepluszz,Mpbh)
        integral = integrate.simps( dtdz(onepluszz-1.)*rate*onepluszz*onepluszz, logonepluszz )
        flux_vec.append( n_pbh(Mpbh)*integral*c )
    return np.array(flux_vec)

# Compute the approximated flux (does not work well, it has to be revised)
def flux_approx(zmin, zmax, Mpbh, E):  #
    if E>4.*Tpbh(Mpbh):
        integral = blackbody(E, Mpbh)*dtdz(0.)
    else:
        logonepluszz = np.linspace( np.log(1.+zmin), np.log(1.+zmax) , 100)
        onepluszz = np.exp(logonepluszz)
        #integral = integrate.simps( dtdz(onepluszz-1.)*rate*onepluszz*onepluszz, logonepluszz )
        integral = E**2.*27.*np.pi*G_N**2.*Mpbh**2.*integrate.simps( dtdz(onepluszz-1.)*onepluszz**2.*onepluszz*onepluszz, logonepluszz )
    return n_pbh(Mpbh)*integral*c

flux_approx = np.vectorize( flux_approx )

#---
# Galactic component (from 2010.16053)
#---

# Some NFW parameters
r_dot=8.5 # kpc
r_h=200. # kpc
r_s = 20.   # kpc
rhocentral = 0.4 # GeV/cm^3

# Radial distance as a function of the cosinus and line-of-sight distance
def galactocentric_d(l, cos):
    return np.sqrt(r_dot**2.-2.*r_dot*l*cos+l**2.)

# NFW profile
def NFW_profile(l,cos):

    r = galactocentric_d(l, cos)
    return rhocentral*(r_dot/r)*((1.+r_dot/r_s)/(1.+r/r_s))**2.

# Maximum line-of-sight distance
def l_max(cos):
    return np.sqrt(r_h**2. - r_dot**2.*(1.-cos**2.)) + r_dot*cos

# D-factor, \int dl \rho_NFW(r(l))
def D_factor():

    #---------------------------------------- 1st integrate in l:
    cos_vec=np.linspace(-1.,1.,100)
    Integral_l=[]
    for cos in cos_vec:
        #I=integrate.quad(NFW_profile,0.,l_max(cos),args=(cos), epsabs= 1e-10, epsrel= 1e-10,limit=500)[0]
        lvec = np.linspace(0., l_max(cos), 100)
        I = integrate.simps(NFW_profile(lvec, cos), lvec)
        #loglvec = np.linspace(np.log(1.e-5), np.log(l_max(cos)), 100)
        #I = integrate.simps(NFW_profile(np.exp(loglvec), cos)*np.exp(loglvec), loglvec)
        Integral_l.append(I)

    #---------------------------------------- now integrate in cosine:

    #Int_cos= interp1d(cos_vec,Integral_l)
    #Int_total=integrate.quad(lambda cos: Int_cos(cos),cos_vec[0],cos_vec[-1], epsabs= 1e-10, epsrel= 1e-10,limit=500)[0]
    Galactic_contrib = integrate.simps( Integral_l, cos_vec)/2.  # (factor of 2 comes from 2pi of integral over 4pi factor of Eq 4.)

    #--- we must add a kpc to cm factor to match units of flux:
    return MpcToCm*1.e-3*Galactic_contrib/gr_to_GeV     # this has units of g/cm^2

galacticfactor = D_factor()

# Galactic flux (f_PBH=1, scale it afterwards)
def galactic_flux(Mpbh, energyrate):

    return energyrate*galacticfactor/Mpbh


# Routine to compute fluxes for a range of PBH masses
def compute_flux(Mpbhs, fpbhs, plot_fluxes = 0):

    onefile = []

    fluxes_max = []
    for mm, Mpbh in enumerate(Mpbhs):

        folder = "BlackHawkData/{:.1e}/".format(Mpbh)

        if not os.path.exists("folder_fluxes/{:.1e}".format(Mpbh)):
            os.system("mkdir "+"folder_fluxes/{:.1e}".format(Mpbh))

        #-----------
        # INST FILES
        #-----------

        data_primary = np.genfromtxt(folder+"instantaneous_primary_spectra.txt",skip_header = 2)
        data_secondary = np.genfromtxt(folder+"instantaneous_secondary_spectra.txt",skip_header = 2)
        E_prim = data_primary[:,0]
        E_sec = data_secondary[:,0]

        tot_sec = 0.
        for i in [2,3,4]:   # three neutrino species
            tot_sec += data_secondary[:,i+1]

        """plt.loglog(E_sec, data_secondary[:, 3],label="e")
        plt.loglog(E_sec, data_secondary[:, 4],label=r"$\mu$")
        plt.loglog(E_sec, data_secondary[:, 5],label=r"$\tau$")
        plt.ylim(1e18,1.e25)
        plt.xlim(1e-3,1)
        plt.legend()
        plt.show()
        exit()"""

        #flux_max = max(tot_sec)
        #plt.loglog(E_sec,tot_sec, linestyle="--", color=cols[mm])
        #plt.loglog(E_prim,data_primary[:,6],label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+" g", color=cols[mm])

        spec_prim = interp1d(E_prim,data_primary[:,6],fill_value="extrapolate")
        spec_sec = interp1d(E_sec,tot_sec,fill_value="extrapolate")
        #spec_prim = np.vectorize(spec_prim)
        #spec_sec = np.vectorize(spec_sec)

        zmin = 0.
        if Mpbh<Mevap:
            zmin = zevap(Mpbh)
        zmax = (1.+zmin)*1.e5 - 1.
        print("Mass: {:.1e}, z min: {:.1e}".format( Mpbh, zmin ) )

        flux_prim = flux(zmin, zmax, Mpbh, E_prim, spec_prim)
        flux_sec = flux(zmin, zmax, Mpbh, E_sec, spec_sec)

        if Mpbh>Mevap:
            flux_galac = galactic_flux(Mpbh, spec_sec(E_sec))
            flux_sec += flux_galac

        fluxes_max.append( max(flux_sec*fpbhs[mm]) )

        if Mpbh<Mevap:
            fpbhlabel = r" g, $\beta'=$"
        else:
            fpbhlabel = r" g, $f_{\rm PBH}=$"

        if plot_fluxes:

            #plt.loglog( E_prim, fpbhs[mm]*flux_prim, color = cols[mm], label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]) )
            plt.loglog( E_sec, fpbhs[mm]*flux_sec, color = cols[mm], linestyle="-", label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+fpbhlabel+scinot(fpbhs[mm]) )
            #plt.loglog( Evec, fpbhs[mm]*flux_anal(Evec, Mpbh), color = cols[mm], linestyle=":" )#, label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+r" g, $f_{\rm PBH}=$"+scinot(fpbhs[mm]) )
            #plt.loglog( Evec, fpbhs[mm]*flux_approx(Evec, Mpbh), linestyle="--", color = cols[mm] )

        #np.savetxt("BlackHawkData/{:.1e}/flux.txt".format(Mpbh), np.transpose([E_sec, flux_sec]) )

        np.savetxt("folder_fluxes/{:.1e}/flux.txt".format(Mpbh), np.transpose([E_sec, flux_sec]) )

        # For masses above ~2.e15, instantaneous flux is equal to the total one
        #if Mpbh<Mevap:
        #if Mpbh<3.e15:
        if Mpbh<3.e33:

            #-----------
            # TOT FILES
            #-----------

            spec_tot_prim = np.genfromtxt(folder + "neutrinos_primary_spectrum.txt",skip_header = 1)
            spec_tot_e = np.genfromtxt(folder + "nu_e_secondary_spectrum.txt",skip_header = 1)
            spec_tot_mu = np.genfromtxt(folder + "nu_mu_secondary_spectrum.txt",skip_header = 1)
            spec_tot_tau = np.genfromtxt(folder + "nu_tau_secondary_spectrum.txt",skip_header = 1)
            Evec = spec_tot_e[0,1:]
            timevec = spec_tot_e[1:,0]
            Evec_prim = spec_tot_prim[0,1:]
            timevec_prim = spec_tot_prim[1:,0]
            spec_tot = spec_tot_e + spec_tot_mu + spec_tot_tau

            for it, t in enumerate(timevec):
                if timevec[it-1]==t:
                    finindex = it
                    break

            """plt.loglog(Evec, spec_tot[150,1:],"m:",lw=4,alpha=0.5,label=r"$t=${:.1e}".format(timevec[2]))
            plt.loglog(Evec, spec_tot[150,1:],"r:",lw=4,alpha=0.5,label=r"$t=${:.1e}".format(timevec[150]))
            plt.loglog(Evec, spec_tot[finindex,1:],"g:",lw=4,alpha=0.5,label=r"$t=${:.1e}".format(timevec[finindex]))
            plt.loglog(Evec, spec_tot[-1,1:],"b:",lw=4,alpha=0.5,label=r"$t=${:.1e}".format(timevec[-1]))
            plt.loglog(Evec_prim, spec_tot_prim[2,1:],"m")
            plt.loglog(Evec_prim, spec_tot_prim[150,1:],"r")
            plt.loglog(Evec_prim, spec_tot_prim[finindex,1:],"g")
            plt.loglog(Evec_prim, spec_tot_prim[-1,1:],"b")
            #plt.ylim(1.e-29, 1e30)
            plt.ylim(1.e-10, 1e25)
            plt.xlim(1., 1e15)
            plt.legend()
            plt.xlabel('$E{\\rm \,\, (GeV)}$')
            plt.ylabel('${\\rm d}N/d E dt \,\, ({\\rm GeV}^{-1}{\\rm s}^{-1})$')
            plt.savefig("figures/dNdEdt_test.pdf", bbox_inches='tight')
            plt.show()
            exit()"""

            #xx, yy = np.meshgrid(Evec, timevec)
            #d2NdEdt_ts = Rbf(np.log10(xx), np.log10(yy), np.log10(np.transpose(spec_tot[1:,1:])),kind="linear")
            #print(d2NdEdt_ts)
            #exit()

            timevec = timevec[timevec<=ageuniverse]

            reds = z_from_t_int(timevec)
            #reds = z_from_t(timevec)
            #reds = redshift(timevec)

            d2NdEdt_ts = []
            wi = []
            for it, time in enumerate(timevec):
                d2NdEdt = spec_tot[1+it,1:]
                d2NdEdt_time = interp1d(Evec,d2NdEdt,fill_value="extrapolate")
                d2NdEdt_time_prim = interp1d(Evec_prim,spec_tot_prim[1+it,1:],fill_value="extrapolate")

                #wi.append( blackbody(Evec*(1.+redshift(time)), Mpbh) )
                #we = dNdEdt_extension(Evec[-1]*(1.+reds[it]),d2NdEdt_time,Evec*(1.+reds[it]),Mpbh)
                we = d2NdEdt_time(Evec*(1.+reds[it]))

                #we = np.zeros_like(Evec)

                for iE, EE in enumerate(Evec):
                    if EE*(1.+reds[it])>Evec[-1]:
                        #we[iE] = blackbody(EE*(1.+reds[it]), Mpbh)
                        we[iE] = d2NdEdt_time_prim(EE*(1.+reds[it]))

                    #if EE*(1.+redshift(time))>4.*Tpbh(Mpbh):
                    #if EE*(1.+reds[it])>4.*Tpbh(Mpbh):
                    #    we[iE] = min(we[iE], blackbody(EE*(1.+reds[it]), Mpbh))

                    """if EE*(1.+reds[it])>Evec[-1]:
                        we[iE] = blackbody(EE*(1.+reds[it]), Mpbh)
                    else:
                        we[iE] = d2NdEdt_time(EE*(1.+reds[it]))"""
                    #if EE*(1.+redshift(time))>Evec[-1]:
                    #    we[iE] = blackbody(EE*(1.+redshift(time)), Mpbh)

                d2NdEdt_ts.append(  we )
            d2NdEdt_ts = np.array(d2NdEdt_ts)


            finite_differences = 0

            """plt.loglog(Evec, d2NdEdt_ts[150,:],":",lw=4,alpha=0.5)
            plt.loglog(Evec, d2NdEdt_ts[250,:],":",lw=4,alpha=0.5)
            plt.loglog(Evec, wi[150])
            plt.loglog(Evec, wi[250])
            plt.show()
            exit()"""

            flux_tot = []
            for j, EE in enumerate(Evec):
                #integrand = d2NdEdt_ts[:,j]*(1.+redshift(timevec))
                integrand = d2NdEdt_ts[:,j]*(1.+reds)
                if finite_differences:
                    integral = 0.
                    for it, t in enumerate(timevec[:-2]):
                        integral += (timevec[it+1]-timevec[it])*integrand[it]*np.heaviside(ageuniverse-t,0.)
                else:
                    integral = integrate.simps( integrand[:finindex]*timevec[:finindex]*np.heaviside(ageuniverse-timevec[:finindex],0.), np.log(timevec[:finindex]) )

                #print(np.mean(d2NdEdt_ts[:,j]*(1.+redshift(timevec))*timevec))
                flux_tot.append( n_pbh(Mpbh)*integral*c )
                #print(integral)
            flux_tot = np.array(flux_tot)

            if Mpbh>Mevap:
                #spec_tot_today = spec_sec(Evec)
                ind = find_nearest(spec_tot[1:,0], ageuniverse, axis=0)
                spec_tot_today = spec_tot[1+ind,1:]
                flux_galac = galactic_flux(Mpbh, spec_tot_today)
                flux_tot += flux_galac

            if plot_fluxes:
                if Mpbh>Mevap:
                    flux_sec = flux_tot-flux_galac
                    # Change units to MeV
                    plt.loglog( Evec*1.e3, fpbhs[mm]*flux_sec/1.e3, linestyle="--", color = cols[mm])
                    plt.loglog( Evec*1.e3, fpbhs[mm]*flux_galac/1.e3, linestyle=":", color = cols[mm])
                    plt.loglog( Evec*1.e3, fpbhs[mm]*flux_tot/1.e3, linestyle="-", color = cols[mm])
                else:
                    plt.loglog( Evec*1.e3, fpbhs[mm]*flux_tot/1.e3, color = cols[mm], linestyle="-" )   # Change units to MeV

            np.savetxt("folder_fluxes/{:.1e}/flux.txt".format(Mpbh), np.transpose([Evec, flux_tot]) )

            onefile.extend(list(flux_tot))

    masses = []
    for Mpbh in Mpbhs:
        masses.extend(list(np.tile(Mpbh, len(Evec))))

    np.savetxt("data/totalflux.txt", np.transpose([np.array(masses), np.tile(Evec, len(Mpbhs)), np.array(onefile)]) )

    return fluxes_max
