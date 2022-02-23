#-------------------------------
# Compute the neutrino fluxes
#-------------------------------

from Source.evaporation import *

#-------
# Generic stuff
#-------

# Number density of PBHs n_pbh at z = 0 (cm^-3)
# beta and f_pbh defined outside of the function (i.e., f_PBH=1 or beta'=1 for this formula, scale it later)
def n_pbh(Mpbh, as_DM, mass_spec, sig=0):
    # For monochromatic
    if mass_spec==0:
        if as_DM:
            return Om_dm*rho_c/Mpbh
        else:
            # See Carr 2010
            return (GpToCm)**(-3.)/(7.98e-29)*(Mpbh/Msun)**(-3./2.)
    # For lognormal, employ normalization such that BlackHawk parameter amplitude_lognormal = 1
    elif mass_spec==1:
        Mmean = Mpbh*np.exp(sig**2./2.)
        if as_DM:
            return Om_dm*rho_c/Mmean
        else:
            return (GpToCm)**(-3.)/(7.98e-29)*(Mpbh/Msun)**(-3./2.)*np.exp(-9./8.*sig**2.)


# Modify the flux for each flavour due to neutrino oscillations (oscillation angles from 2006.11237)
def flux_oscillations(F_nue, F_numu, F_nutau):
    Fprime_nue = 0.542396*F_nue + 0.176253*F_numu + 0.281351*F_nutau
    Fprime_numu = 0.176253*F_nue + 0.451522*F_numu + 0.372224*F_nutau
    Fprime_nutau = 0.281351*F_nue + 0.372224*F_numu + 0.346425*F_nutau
    return Fprime_nue, Fprime_numu, Fprime_nutau

# Sufix for text files, depending on if mass function is monochromatic (mc) or lognormal (ln), specifying sigma in the later case
def sufix(mass_spec, sig = 0):
    if mass_spec==0:
        sufix = "_mc.txt"
    elif mass_spec==1:
        sufix = "_ln_sig_{:.1e}.txt".format(sig)
    return sufix

#-------
# Galactic component (from 2010.16053)
#-------

# Some NFW parameters
r_dot = 8.5 # kpc
r_h = 200. # kpc
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

    #-------------------------------------------- 1st integrate in l:
    cos_vec=np.linspace(-1.,1.,100)
    Integral_l=[]
    for cos in cos_vec:
        #I=integrate.quad(NFW_profile,0.,l_max(cos),args=(cos), epsabs= 1e-10, epsrel= 1e-10,limit=500)[0]
        lvec = np.linspace(0., l_max(cos), 100)
        I = integrate.simps(NFW_profile(lvec, cos), lvec)
        #loglvec = np.linspace(np.log(1.e-5), np.log(l_max(cos)), 100)
        #I = integrate.simps(NFW_profile(np.exp(loglvec), cos)*np.exp(loglvec), loglvec)
        Integral_l.append(I)

    #-------------------------------------------- now integrate in cosine:

    #Int_cos= interp1d(cos_vec,Integral_l)
    #Int_total=integrate.quad(lambda cos: Int_cos(cos),cos_vec[0],cos_vec[-1], epsabs= 1e-10, epsrel= 1e-10,limit=500)[0]
    Galactic_contrib = integrate.simps( Integral_l, cos_vec)/2.  # (factor of 2 comes from 2pi of integral over 4pi factor of Eq 4.)

    #------- we must add a kpc to cm factor to match units of flux:
    return MpcToCm*1.e-3*Galactic_contrib/gr_to_GeV     # this has units of g/cm^2

galacticfactor = D_factor()

# Galactic flux (f_PBH=1, scale it afterwards)
def galactic_flux(Mpbh, energyrate, mass_spec, sig=0):

    if mass_spec==0:
        galflux = energyrate*galacticfactor/Mpbh
    elif mass_spec==1:
        Mmean = Mpbh*np.exp(sig**2./2.)
        galflux = energyrate*galacticfactor/Mmean
    return galflux

#-------
# Extragalactic flux (only for instantaneous)
#-------

# Compute the diffuse flux
def flux(zmin, zmax, Mpbh, E_vec, spec_int, as_DM, mass_spec, sig=0, aprox=0):
    flux_vec = []
    logonepluszz = np.linspace( np.log(1.+zmin), np.log(1.+zmax) , 100)
    onepluszz = np.exp(logonepluszz)
    for E in E_vec:
        if aprox:
            rate = blackbody(E*onepluszz, Mpbh)
        else:
            rate = dNdEdt_extension(E_vec[-1],spec_int,E*onepluszz,Mpbh)
        integral = integrate.simps( dtdz(onepluszz-1.)*rate*onepluszz*onepluszz, logonepluszz )
        flux_vec.append( n_pbh(Mpbh, as_DM, mass_spec, sig)*integral*c )
    return np.array(flux_vec)

# Compute the approximated flux (not employed, does not work well, it has to be revised)
def flux_approx(zmin, zmax, Mpbh, E, as_DM, mass_spec):  #
    if E>4.*Tpbh(Mpbh):
        integral = blackbody(E, Mpbh)*dtdz(0.)
    else:
        logonepluszz = np.linspace( np.log(1.+zmin), np.log(1.+zmax) , 100)
        onepluszz = np.exp(logonepluszz)
        #integral = integrate.simps( dtdz(onepluszz-1.)*rate*onepluszz*onepluszz, logonepluszz )
        integral = E**2.*27.*np.pi*G_N**2.*Mpbh**2.*integrate.simps( dtdz(onepluszz-1.)*onepluszz**2.*onepluszz*onepluszz, logonepluszz )
    return n_pbh(Mpbh, as_DM, mass_spec, sig)*integral*c

flux_approx = np.vectorize( flux_approx )


#-------
# Routine to compute fluxes for a range of PBH masses
#-------

# as_DM: 1 for PBHs as DM and use f_PBH, otherwise, it uses beta'
# mass_spec: mass spectrum, 0 for monochromatic, 1 for lognormal
# sig: variance for lognormal mass spectrum (only used if mass_spec=1)
# use_inst: if 1, use instantaneous Blackhawk tables (not recommended for masses <2e15), otherwise it employs total tables
# flux normalization assumes fpbh=1 (for PBHs as DM, as_DM=1) or beta'=1 (for as_DM=0)
def compute_flux(Mpbhs, as_DM, mass_spec = 0, sig = 0, use_inst = 0):

    sufx = sufix(mass_spec, sig)

    # This will be filled by PBH mass, energy and flux, not used in this code, for exporting only
    onefile = []

    for mm, Mpbh in enumerate(Mpbhs):

        print("Mass: {:.1e} g".format( Mpbh ) )

        folder = folder_blackhawk+"{:.1e}/".format(Mpbh)

        if not os.path.exists("fluxes/{:.1e}".format(Mpbh)):
            os.system("mkdir "+"fluxes/{:.1e}".format(Mpbh))

        if use_inst:

            #---------------
            # Instantaneous spectra
            #---------------

            data_primary = np.genfromtxt(folder+"instantaneous_primary_spectra"+sufx, skip_header = 2)
            data_secondary = np.genfromtxt(folder+"instantaneous_secondary_spectra"+sufx, skip_header = 2)
            E_prim = data_primary[:,0]
            Evec = data_secondary[:,0]

            #tot_sec = 0.
            #for i in [2,3,4]:   # three neutrino species
            #    tot_sec += data_secondary[:,i+1]

            spec_tot_e, spec_tot_mu, spec_tot_tau = flux_oscillations(data_secondary[:,3], data_secondary[:,4], data_secondary[:,5])
            tot_sec = spec_tot_e/2. # 1/2 for taking only neutrino (or antineutrino)

            spec_prim = interp1d(E_prim,data_primary[:,6],fill_value="extrapolate")
            spec_sec = interp1d(Evec,tot_sec,fill_value="extrapolate")

            zmin = max([0., zevap(Mpbh)])
            # Take an arbitrary large maximum z
            zmax = (1.+zmin)*1.e5 - 1.

            flux_prim = flux(zmin, zmax, Mpbh, E_prim, spec_prim, as_DM, mass_spec, sig)/1.e3     # Change units to MeV (factor 1.e3)
            flux_sec = flux(zmin, zmax, Mpbh, Evec, spec_sec, as_DM, mass_spec, sig)/1.e3     # Change units to MeV (factor 1.e3)

            if Mpbh>=Mevap:
                flux_galac = galactic_flux(Mpbh, spec_sec(Evec), mass_spec, sig)/1.e3     # Change units to MeV (factor 1.e3)
                flux_sec += flux_galac

            # Change units to MeV (factor 1.e3)
            Evec = Evec*1.e3

            np.savetxt("fluxes/{:.1e}/flux_isDM_{}".format(Mpbh, as_DM)+sufx, np.transpose([Evec, flux_sec]) )
            if Mpbh>=Mevap:
                np.savetxt("fluxes/{:.1e}/flux_galac_isDM_{}".format(Mpbh, as_DM)+sufx, np.transpose([Evec, flux_galac]) )
                np.savetxt("fluxes/{:.1e}/flux_extragalac_isDM_{}".format(Mpbh, as_DM)+sufx, np.transpose([Evec, flux_sec - flux_galac]) )

            onefile.extend(list(flux_sec))

            # For masses above ~2.e15, instantaneous flux is equal to the total one
        else:

            #---------------
            # Total spectra
            #---------------

            # Get the neutrino spectra
            # Secondary spectra in BlackHawk already include the primary ones
            spec_tot_e = np.genfromtxt(folder + "nu_e_secondary_spectrum"+sufx,skip_header = 1)
            spec_tot_mu = np.genfromtxt(folder + "nu_mu_secondary_spectrum"+sufx,skip_header = 1)
            spec_tot_tau = np.genfromtxt(folder + "nu_tau_secondary_spectrum"+sufx,skip_header = 1)
            Evec = spec_tot_e[0,1:]
            timevec = spec_tot_e[1:,0]

            # Take into account oscillations
            spec_tot_e[1:,1:], spec_tot_mu[1:,1:], spec_tot_tau[1:,1:] = flux_oscillations(spec_tot_e[1:,1:], spec_tot_mu[1:,1:], spec_tot_tau[1:,1:])
            # Consider only electronic neutrinos
            spec_tot = spec_tot_e
            spec_tot[1:,1:] = spec_tot[1:,1:]/2.# 1/2 for taking only neutrino (or antineutrino)

            # BlackHawk files often provide a repeated time at some rows, because the precision at writing the output is not enough
            # In such cases, take only the first row among those with the same time value
            indexes = []
            for it, t in enumerate(timevec[:-2]):
                if timevec[it+1]<t: print("Array not sorted")
                if timevec[it+1]>ageuniverse:
                    break
                elif not timevec[it+1]==t:
                    indexes.append( it+1 )

            # Compute the redshifted spectrum
            d2NdEdt_ts = []
            for it in indexes:
                time_it = timevec[it]

                d2NdEdt = spec_tot[1+it,1:]

                # To avoid numerical problems in the interpolation
                d2NdEdt[d2NdEdt==0.] = 1.e-300

                logd2NdEdt_time = interp1d(np.log(Evec),np.log(d2NdEdt), fill_value=-300, bounds_error=False)
                rateredshift = np.exp(logd2NdEdt_time(np.log(Evec*(1.+z_from_t_int(timevec[it])))))

                d2NdEdt_ts.append( rateredshift )

            d2NdEdt_ts = np.array(d2NdEdt_ts)

            timevec = timevec[indexes]

            reds = z_from_t_int(timevec)

            # Integrate the flux until z=0 or until PBHs are completely evaporated
            flux_tot = []
            for j, EE in enumerate(Evec):
                integrand = d2NdEdt_ts[:timevec.shape[0],j]*(1.+reds)
                # Introduce a step function to finish the integral at the current age of the universe
                integral = integrate.simps( integrand*timevec*np.heaviside(ageuniverse-timevec,0.), np.log(timevec) )
                flux_tot.append( n_pbh(Mpbh, as_DM, mass_spec, sig)*integral*c )

            # Change units to MeV (factor 1.e3)
            Evec, flux_tot = Evec*1.e3, np.array(flux_tot)/1.e3

            # If PBHS are DM, include galactic contribution
            if Mpbh>=Mevap:
                # Find the spectrum evaluated at the current age of the universe
                ind = find_nearest(spec_tot[1:,0], ageuniverse, axis=0)
                spec_tot_today = spec_tot[1+ind,1:]/1.e3    # Change units to MeV (factor 1.e3)
                flux_galac = galactic_flux(Mpbh, spec_tot_today, mass_spec, sig)
                flux_tot += flux_galac

            np.savetxt("fluxes/{:.1e}/flux_isDM_{}".format(Mpbh, as_DM)+sufx, np.transpose([Evec, flux_tot]) )
            if Mpbh>=Mevap:
                np.savetxt("fluxes/{:.1e}/flux_galac_isDM_{}".format(Mpbh, as_DM)+sufx, np.transpose([Evec, flux_galac]) )
                np.savetxt("fluxes/{:.1e}/flux_extragalac_isDM_{}".format(Mpbh, as_DM)+sufx, np.transpose([Evec, flux_tot - flux_galac]) )

            onefile.extend(list(flux_tot))

    masses = []
    for Mpbh in Mpbhs:
        masses.extend(list(np.tile(Mpbh, len(Evec))))

    # File with PBH mass, energy and flux, not used in this code, for exporting only
    np.savetxt("fluxes/totalflux_Mpbh_from_{:.1e}_to_{:.1e}".format(Mpbhs[0], Mpbhs[-1])+sufx, np.transpose([np.array(masses), np.tile(Evec, len(Mpbhs)), np.array(onefile)]) )
