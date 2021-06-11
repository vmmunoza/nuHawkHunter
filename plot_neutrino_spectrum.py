## 1 - Necessary importations

import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d


import matplotlib as mpl
mpl.rcParams.update({'font.size': 12})

def scinot(x):      # write in Latex scientific notation
    exp = int(math.floor(math.log10(abs(x))))
    prefactor = x / 10**exp
    if prefactor == 1.:
        return r"$10^{"+str(exp)+r"}$"
    else:
        return "{:.0f}".format(prefactor)+r"$ \times 10^{"+str(exp)+r"}$"

# Black-body approximation

mass_conversion = 5.60958884e+23 # g to GeV
time_conversion = 1.519267407e+24 # s to GeV^(-1)
leng_conversion = 5.06773058e+13 # cm to GeV^(-1)
G_N = (6.67408e-11*pow(leng_conversion*100.,3.)/(mass_conversion*1000.)/pow(time_conversion,2.))

def gamma(E, Mpbh):
    Tpbh = 1.06e13/Mpbh   # GeV, Mpbh in g
    if E>4.*Tpbh:
        return E**2.*27.*np.pi*G_N**2.*Mpbh**2.
    else:
        return E**2.*2.*G_N**2.*Mpbh**2.

gamma = np.vectorize(gamma)

def blackbody(E, Mpbh):
    Tpbh = 1.e13/Mpbh   # GeV, Mpbh in g
    dNdtdE = gamma(E, Mpbh)*( np.exp( E/Tpbh ) +1. )**(-1.)/(2.*np.pi)
    return 6.*dNdtdE*mass_conversion**2.*time_conversion



## 2 - Folder definition

# Here put the BlackHawk path
BH_path = ""
redshift = 0.
#Mpbhs = [1.e8, 1.e10, 1.e11, 1.e12, 1.e15]
#Mpbhs = [1.e10, 1.e11, 1.e12]
Mpbhs = [1.e13, 1.e14, 1.e15]
#Mpbhs = [1.e15, 2.e15]
cols = ["r", "b", "g", "m", "y"]

result_folder = BH_path + "results/"
# Here put a path to some figures folder
fig_folder = BH_path+"figures/"

for mm, Mpbh in enumerate(Mpbhs):


    # Here put the name of the 'destination_folder'
    data_folder = "{:.1e}".format(Mpbh)
    folder = "BlackHawkData/" + data_folder + "/"

    # 3 - Recovering data

    #data_spectrum = np.genfromtxt(folder+"BH_spectrum.txt",skip_header = 2)
    data_primary = np.genfromtxt(folder+"instantaneous_primary_spectra.txt",skip_header = 2)
    data_secondary = np.genfromtxt(folder+"instantaneous_secondary_spectra.txt",skip_header = 2)

    # 5 - Plotting options (primary data)

    labels_primary=np.array(np.zeros(16),str)
    labels_primary[5]="$\\nu,\overline{\\nu}$"

    # 7 - Plotting options (secondary data)
    epoch = 1 # 0: BBN epoch, 1: present epoch

    if epoch == 0:
        nb_fin_part = 11

        # Put 1 to plot the particle spectrum

        nu_e_secondary=1
        nu_mu_secondary=1
        nu_tau_secondary=1


        part_show_secondary=np.zeros(nb_fin_part)
        part_show_secondary=np.array(part_show_secondary,int)

        part_show_secondary[3] = nu_e_secondary
        part_show_secondary[4] = nu_mu_secondary
        part_show_secondary[5] = nu_tau_secondary


        labels_secondary=np.array(np.zeros(nb_fin_part),str)

        labels_secondary[3]="$\\nu_e,\overline{\\nu}_e$"
        labels_secondary[4]="$\\nu_\mu,\overline{\\nu}_\mu$"
        labels_secondary[5]="$\\nu_\\tau,\overline{\\nu}_\\tau$"


    elif epoch == 1:
        nb_fin_part = 6

        # Put 1 to plot the particle spectrum

        nu_e_secondary=1
        nu_mu_secondary=1
        nu_tau_secondary=1

        part_show_secondary=np.zeros(nb_fin_part)
        part_show_secondary=np.array(part_show_secondary,int)
        part_show_secondary[2] = nu_e_secondary
        part_show_secondary[3] = nu_mu_secondary
        part_show_secondary[4] = nu_tau_secondary

        labels_secondary=np.array(np.zeros(nb_fin_part),str)

        labels_secondary[2]="$\\nu_e,\overline{\\nu}_e$"
        labels_secondary[3]="$\\nu_\mu,\overline{\\nu}_\mu$"
        labels_secondary[4]="$\\nu_\\tau,\overline{\\nu}_\\tau$"

    # 8 - Plotting secondary spectra

    #plt.figure(3)
    #plt.clf()
    #plt.axes([0.125,0.2,0.90 - 0.125,0.90-0.2])
    flux_max = 0.
    for i in range(nb_fin_part):
        if part_show_secondary[i]:
            flux_max = max(flux_max,max(data_secondary[:,i+1]))
    plt.ylim(flux_max/1e+10,flux_max*10.)
    #plt.xlim(data_primary[0,0],data_primary[-1,0])
    """
    for i in range(nb_fin_part):
        if part_show_secondary[i]:
            plt.loglog(data_secondary[:,0]/(1.+redshift),data_secondary[:,i+1],label = labels_secondary[i],linewidth = 2,linestyle="--")
    """
    tot_sec=0.
    for i in [2,3,4]:   # three neutrino species
        tot_sec += data_secondary[:,i+1]

    Esec = data_secondary[:,0]*1.e3 # in MeV
    intsec = interp1d(Esec, tot_sec/1.e3, fill_value="extrapolate")     # in MeV
    intprim = interp1d(data_primary[:,0]*1.e3,data_primary[:,6]/1.e3,fill_value="extrapolate")     # in MeV

    plt.loglog(Esec,intprim(Esec),linestyle=":",linewidth = 2, color=cols[mm])
    plt.loglog(Esec,intsec(Esec)-intprim(Esec),"--",linewidth = 2, color=cols[mm])
    plt.loglog(Esec,intsec(Esec),"-",linewidth = 2, color=cols[mm])

    #plt.loglog(data_secondary[:,0],tot_sec,"--",linewidth = 2, color=cols[mm])
    #plt.loglog(data_primary[:,0],data_primary[:,6],label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+" g",linewidth = 2, color=cols[mm])
    #plt.loglog(data_primary[:,0], blackbody(data_primary[:,0], Mpbh), linestyle=":",linewidth = 2, color=cols[mm])

customlegend = []
for n, Mpbh in enumerate(Mpbhs):
    customlegend.append( Line2D([0], [0], color=cols[n], lw=4, label = r"$M_{\rm PBH}=$"+scinot(Mpbh)+" g"))

customlegend.append( Line2D([0], [0], color="black", linestyle=":", label="Primary"))
customlegend.append( Line2D([0], [0], color="black", linestyle="--", label="Secondary"))
customlegend.append( Line2D([0], [0], color="black", linestyle="-", label="Total"))



#plt.ylim(1.e2,1.e25)
#plt.xlim(data_primary[0,0],data_primary[-1,0])
plt.ylim(1.e16,1.e23)
plt.xlim(1.e0,1.e5)
plt.xlabel('$E{\\rm \,\, [MeV]}$')
#plt.ylabel('${\\rm d}N/d E dt \,\, [{\\rm MeV}^{-1}{\\rm s}^{-1}]$')
plt.ylabel(r'$\frac{dN}{dE dt} \,\, [{\rm MeV}^{-1}{\rm s}^{-1}]$')
#plt.ylabel('${\\rm d}^2n/{\\rm d}t{\\rm d}E\,\, ({\\rm GeV}^{-1}\cdot{\\rm s}^{-1}\cdot{\\rm cm}^{-3})$')
#title = '${\\rm Secondary\,\, spectra\,\,}-\,a^*_i=%.1f,\,M=%.3e\,{\\rm g}$'%(float(data_folder[:7]),float(data_folder[8:]))
#plt.title(title,fontsize = 30,y = 1.02)
plt.grid()
plt.legend(handles=customlegend)

fig_name = fig_folder +  "neutrino_spectra_blackhawk.pdf"

plt.savefig(fig_name)
plt.show()
## 9 - Close all figures
plt.close("all")
