import os, fileinput, sys, time
import numpy as np
from scipy.interpolate import interp1d

time_ini = time.time()
path = os.getcwd()

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)


#Mpbhs = [1.e12, 1.e13, 1.e14]
#Mpbhs = [1.e15, 2.e15, 4.e15]
Mpbhs = np.logspace(np.log10(1.e12), np.log10(7.e14), 15)
#Mpbhs = np.logspace(np.log10(8.e14), np.log10(1.e16), 15)

# Choose mass spectrum: 0 = monochromatic, 1 = lognormal
spectrum_choice = 0

# For lognormal:
sigma = 1.
BHnumber = 50  # number of masses

param_file = path + '/parameters.txt'

#--- MAIN ---#

if not os.path.exists("BlackHawkData"):
    os.system("mkdir BlackHawkData")

for M in Mpbhs:

    # Sufix for text files, depending on if mass function is monochromatic (mc) or lognormal (ln), specifying sigma in the later case
    if spectrum_choice==0:
        sufix = "_mc.txt"
    elif spectrum_choice==1:
        sufix = "_ln_sig_{:.1e}.txt".format(sigma)

    Tpbh = 1.06e13/M   # GeV, Mpbh in g
    Emin = 1.e-2*Tpbh
    Emax = 1.e5*Tpbh
    Enumber = 150

    replaceAll(param_file, "destination_folder =", "destination_folder = {:.1e} \t\t\t\t\t\t\t\t\t # name of the output folder in results/ \n".format(M))
    replaceAll(param_file, "spectrum_choice = ", "spectrum_choice = {:} \t\t\t\t\t\t\t\t\t # form of the BH distribution: 0=Dirac, 1=log-normal, 2=power-law, 3=critical collapse, 4=peak theory, 5=uniform -1=user-defined \n".format(spectrum_choice))

    if spectrum_choice==0:
        replaceAll(param_file, "Mmin = ", "Mmin = {:.1e} \t\t\t\t\t\t\t\t\t # lowest BH mass in g (larger than the Planck mass) \n".format(M))
        replaceAll(param_file, "BHnumber = ", "BHnumber = {:} \t\t\t\t\t\t\t\t\t # number of BH masses (should be the number of tabulated masses if spectrum_choice=5) \n".format(1))


    if spectrum_choice==1:
        replaceAll(param_file, "BHnumber = ", "BHnumber = {:} \t\t\t\t\t\t\t\t\t # number of BH masses (should be the number of tabulated masses if spectrum_choice=5) \n".format(BHnumber))
        replaceAll(param_file, "Mmin = ", "Mmin = {:.1e} \t\t\t\t\t\t\t\t\t # lowest BH mass in g (larger than the Planck mass) \n".format(M/10.))
        replaceAll(param_file, "Mmax = ", "Mmax = {:.1e} \t\t\t\t\t\t\t\t\t # highest BH mass in g (larger than the Planck mass) \n".format(M*10.))
        replaceAll(param_file, "amplitude_lognormal = ", "amplitude_lognormal = {:.1e} \t\t\t\t\t\t\t\t\t # amplitude of the log-normal (mass density) distribution in g.cm^-3 \n".format( 1. ) ) #M*np.exp(-sigma**2./2.) ))
        replaceAll(param_file, "stand_dev_lognormal = ", "stand_dev_lognormal = {:.1e} \t\t\t\t\t\t\t\t\t # dimensionless variance of the log-normal distribution \n".format(sigma))
        replaceAll(param_file, "crit_mass_lognormal = ", "crit_mass_lognormal = {:.1e} \t\t\t\t\t\t\t\t\t # # characteristic mass of the log-normal distribution in g \n".format(M))


    replaceAll(param_file, "Emin = ", "Emin = {:.1e} \t\t\t\t\t\t\t\t\t # minimal energy in GeV of the primary particles \n".format(Emin))
    replaceAll(param_file, "Emax = ", "Emax = {:.1e} \t\t\t\t\t\t\t\t\t # maximal energy in GeV of the primary particles \n".format(Emax))
    replaceAll(param_file, "Enumber = ", "Enumber = {:} \t\t\t\t\t\t\t\t\t # number of primary particles energies to be simulated \n".format(Enumber))

    # Run
    #os.system("./BlackHawk_inst.x parameters.txt ")
    os.system("./BlackHawk_tot.x parameters.txt ")

    # Save only relevant data
    if not os.path.exists("BlackHawkData/{:.1e}/".format(M)):
        os.system("mkdir BlackHawkData/{:.1e}/".format(M))
    #os.system("cp results/{:.1e}/instantaneous_primary_spectra.txt BlackHawkData/{:.1e}/instantaneous_primary_spectra".format(M,M)+sufix)
    #os.system("cp results/{:.1e}/instantaneous_secondary_spectra.txt BlackHawkData/{:.1e}/instantaneous_secondary_spectra".format(M,M)+sufix)
    os.system("cp results/{:.1e}/neutrinos_primary_spectrum.txt BlackHawkData/{:.1e}/neutrinos_primary_spectrum".format(M,M)+sufix)
    os.system("cp results/{:.1e}/nu_e_secondary_spectrum.txt BlackHawkData/{:.1e}/nu_e_secondary_spectrum".format(M,M)+sufix)
    os.system("cp results/{:.1e}/nu_mu_secondary_spectrum.txt BlackHawkData/{:.1e}/nu_mu_secondary_spectrum".format(M,M)+sufix)
    os.system("cp results/{:.1e}/nu_tau_secondary_spectrum.txt BlackHawkData/{:.1e}/nu_tau_secondary_spectrum".format(M,M)+sufix)
    if spectrum_choice==1:
        os.system("cp results/{:.1e}/BH_spectrum.txt BlackHawkData/{:.1e}/BH_spectrum".format(M,M)+sufix)



print("Minutes elapsed:",(time.time()-time_ini)/60.)
