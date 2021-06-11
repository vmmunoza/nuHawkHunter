
import matplotlib.pyplot as plt
import numpy as np
import math
#!python
import pylab
from pylab import arange,pi,sin,cos,sqrt
from scipy import integrate
from scipy.interpolate import interp1d, Rbf
from flux_stuff import *
import os

for path in ["figures", "folder_fluxes"]:
    if not os.path.exists(path):
        os.system("mkdir "+path)

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

    np.savetxt("totalflux.txt", np.transpose([np.array(masses), np.tile(Evec, len(Mpbhs)), np.array(onefile)]) )

    return fluxes_max

if __name__=="__main__":

    plot_fluxes = 1
    plot_DM = 0

    if plot_DM:
        #Mpbhs =  [1e15, 5e15]#, 2e15, 4e15]
        Mpbhs =  [1e15, 2e15, 4e15]
        fpbh = 1.e-4
        sufix = "DM"
    else:
        Mpbhs = [1.e12, 1.e13, 1.e14]
        fpbh = 1.e-20
        sufix = "evaporated"

    #Mpbhs = np.linspace(1.e15, 8.e15, 15)

    fpbhs = fpbh*np.ones_like(Mpbhs)    # this is beta prime for evaporated PBHs
    cols = ["r", "m", "purple", "b", "g", "orange"]

    fluxes_max = compute_flux(Mpbhs, fpbhs, plot_fluxes)

    if plot_fluxes:

        backfolder = "data/backfluxes/"

        # Plot backgrounds
        Eatm, atm_nue = np.loadtxt(backfolder+"atmnue_noosc_fluka_flux_norm.dat",unpack=True)
        Eatm, atm_nuebar = np.loadtxt(backfolder+"atmnuebar_noosc_fluka_flux_norm.dat",unpack=True)
        Eatm, atm_numu = np.loadtxt(backfolder+"atmnumu_noosc_fluka_flux_norm.dat",unpack=True)
        Eatm, atm_numubar = np.loadtxt(backfolder+"atmnumubar_noosc_fluka_flux_norm.dat",unpack=True)
        atmflux = atm_nue + atm_nuebar + atm_numu + atm_numubar
        EB8, sol_B8_1, sol_B8_2, sol_B8_3 = np.loadtxt(backfolder+"B8NeutrinoFlux.dat",unpack=True)
        sol_B8 = sol_B8_1 + sol_B8_2 + sol_B8_3
        Ehep, sol_hep = np.loadtxt(backfolder+"HEPNeutrinoFlux.dat",unpack=True)
        EO15, sol_O15 = np.loadtxt(backfolder+"O15NeutrinoFlux.dat",unpack=True)
        EN13, sol_N13 = np.loadtxt(backfolder+"N13NeutrinoFlux.dat",unpack=True)
        Epp, sol_pp = np.loadtxt(backfolder+"PPNeutrinoFlux.dat",unpack=True)

        atmint = interp1d(Eatm, atmflux, fill_value=0., bounds_error=False)
        B8int = interp1d(EB8, sol_B8, fill_value=0., bounds_error=False)
        hepint = interp1d(Ehep, sol_hep, fill_value=0., bounds_error=False)
        O15int = interp1d(EO15, sol_O15, fill_value=0., bounds_error=False)
        N13int = interp1d(EN13, sol_N13, fill_value=0., bounds_error=False)
        ppint = interp1d(Epp, sol_pp, fill_value=0., bounds_error=False)

        Ebacks = np.logspace(np.log10(Epp[0]), np.log10(Eatm[-1]), 500)
        # Sum backgrounds and correct normalization, see table III of 1812.05550 or Table 2 of 1208.5723
        backmax = atmint(Ebacks) + B8int(Ebacks)*4.59e6 + hepint(Ebacks)*8.31e3 + O15int(Ebacks)*1.56e8 + N13int(Ebacks)*2.17e8 + ppint(Ebacks)*6.03e10

        minback = []
        for back in [atmflux , sol_B8 , sol_hep , sol_O15 , sol_N13]:
            minback.append( np.min(back) )
        minback = np.min(np.array(minback))

        plt.fill_between( Ebacks, np.zeros_like(Ebacks), backmax, color = "b", alpha=0.3)

        plt.text(120., 2.e-3, "Atm.")
        plt.text(20., 1.e1, r"hep")
        plt.text(6., 5.e6, r"$^8$B")
        plt.text(1.7, 1.e7, r"$^{15}$O")

        customlegend = []
        for n, Mpbh in enumerate(Mpbhs):
            if Mpbh<Mevap:
                fpbhlabel = r", $\beta'=$"
            else:
                fpbhlabel = r", $f_{\rm PBH}=$"
            customlegend.append( Line2D([0], [0], color=cols[n], lw=4, label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+" g"))#+fpbhlabel+scinot(fpbhs[mm])))

        if plot_DM:

            customlegend.append( Line2D([0], [0], color="black", linestyle=":", label="Galactic"))
            customlegend.append( Line2D([0], [0], color="black", linestyle="--", label="Extragalactic"))
            customlegend.append( Line2D([0], [0], color="black", linestyle="-", label="Total"))

        customlegend.append( Line2D([0], [0], color="b", lw=6, linestyle="-", alpha=0.3, label="Backgrounds"))

        plt.xlim(1., 2.e2)
        plt.ylim( 1.e-5, 1.e8 )
        #plt.xlim(0.2, 1.e2)
        #plt.ylim( 1.e-2, 1.e4 )
        plt.legend(handles=customlegend, fontsize=10)#, loc="lower left")
        plt.xlabel('$E{\\rm \,\, [MeV]}$')
        if plot_DM:
            plt.title(r"$f_{\rm PBH}=$"+scinot(fpbh))
        else:
            plt.title(r"$\beta'=$"+scinot(fpbh))
        plt.ylabel('$\Phi \,\, [{\\rm MeV}^{-1}{\\rm s}^{-1}{\\rm cm}^{-2}]$')
        plt.savefig("figures/fluxes_"+sufix+".pdf", bbox_inches='tight')
        plt.show()
        plt.gcf().clear()
