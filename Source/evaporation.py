#-------------------------------
# Hawking evaporation functions
#-------------------------------

from Source.cosmo import *


# Hawking temperature, in GeV, Mpbh in g
def Tpbh(Mpbh):
    return 1.06e13/Mpbh

# Grayfactor approximation
def gamma(E, Mpbh):
    if E>4.*Tpbh(Mpbh):
        return E**2.*27.*np.pi*G_N**2.*Mpbh**2.
    else:
        return E**2.*2.*G_N**2.*Mpbh**2.

# Blackbody approximation for spectrum rate
def blackbody(E, Mpbh):
    if E>Tpbh(Mpbh)*100.:
        fdistr =  np.exp( -E/Tpbh(Mpbh) )   # Avoid numerical problems
    else:
        fdistr = ( np.exp( E/Tpbh(Mpbh) ) +1. )**(-1.)
    dNdtdE = gamma(E, Mpbh)*fdistr/(2.*np.pi)
    return 3.*dNdtdE*mass_conversion**2.*time_conversion     # 3 flavours (as in Blackhawk data, I think)
    # 6 degrees of freedom, 3 flavours and particle-antiparticle (as in Blackhawk data, I think)

# Extend interpolated spectrum rate 'interp' with the blackbody approximation for high energies > Elim
def dNdEdt_extension(Elim,interp,E,Mpbh):
    if E>Elim:
        return blackbody(E, Mpbh)
    else:
        return interp(E)

# Approx lifetime in seconds given the PBH mass in g (Carr et al 2010)
def tau_pbh(M):
    return 407.*(M/1e10)**3.

# Evaporation redshift, valid for the matter dominated epoch at z>>1 or for EdS
def zevap(M):
    return redshift(tau_pbh(M))

gamma = np.vectorize(gamma)
blackbody = np.vectorize(blackbody)
dNdEdt_extension = np.vectorize(dNdEdt_extension)
