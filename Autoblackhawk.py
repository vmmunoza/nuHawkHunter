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


Mpbhs = [1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14]
Mpbhs = [1.e15, 2.e15, 4.e15, 8.e15]
Mpbhs = np.linspace(1.e15, 8.e15, 15)
#Mpbhs = [1.e15, 2.e15]

param_file = path + '/parameters.txt'

#--- MAIN ---#

if not os.path.exists("BlackHawkData"):
    os.system("mkdir BlackHawkData")

for M in Mpbhs:

    Tpbh = 1.06e13/M   # GeV, Mpbh in g
    Emin = 1.e-2*Tpbh
    Emax = 1.e5*Tpbh
    Enumber = 150

    replaceAll(param_file, "destination_folder =", "destination_folder = {:.1e} \t\t\t\t\t\t\t\t\t # name of the output folder in results/ \n".format(M))
    replaceAll(param_file, "Mmin = ", "Mmin = {:.1e} \t\t\t\t\t\t\t\t\t # lowest BH mass in g (larger than the Planck mass) \n".format(M))
    replaceAll(param_file, "Emin = ", "Emin = {:.1e} \t\t\t\t\t\t\t\t\t # minimal energy in GeV of the primary particles \n".format(Emin))
    replaceAll(param_file, "Emax = ", "Emax = {:.1e} \t\t\t\t\t\t\t\t\t # maximal energy in GeV of the primary particles \n".format(Emax))
    replaceAll(param_file, "Enumber = ", "Enumber = {:} \t\t\t\t\t\t\t\t\t # number of primary particles energies to be simulated \n".format(Enumber))

    # Run
    os.system("./BlackHawk_inst.x parameters.txt ")
    os.system("./BlackHawk_tot.x parameters.txt ")

    # Save only relevant data
    if not os.path.exists("BlackHawkData/{:.1e}/".format(M)):
        os.system("mkdir BlackHawkData/{:.1e}/".format(M))
    os.system("cp results/{:.1e}/instantaneous_primary_spectra.txt BlackHawkData/{:.1e}/instantaneous_primary_spectra.txt".format(M,M))
    os.system("cp results/{:.1e}/instantaneous_secondary_spectra.txt BlackHawkData/{:.1e}/instantaneous_secondary_spectra.txt".format(M,M))
    os.system("cp results/{:.1e}/neutrinos_primary_spectrum.txt BlackHawkData/{:.1e}/neutrinos_primary_spectrum.txt".format(M,M))
    os.system("cp results/{:.1e}/nu_e_secondary_spectrum.txt BlackHawkData/{:.1e}/nu_e_secondary_spectrum.txt".format(M,M))
    os.system("cp results/{:.1e}/nu_mu_secondary_spectrum.txt BlackHawkData/{:.1e}/nu_mu_secondary_spectrum.txt".format(M,M))
    os.system("cp results/{:.1e}/nu_tau_secondary_spectrum.txt BlackHawkData/{:.1e}/nu_tau_secondary_spectrum.txt".format(M,M))



print("Minutes elapsed:",(time.time()-time_ini)/60.)
