#-------------------------------
# Cosmology functions
#-------------------------------

from Source.constants import *


# Hubble rate, s^{-1}
def hubble(z):
    return H0*np.sqrt( Om_rad*(1.+z)**(4.) + Om_m*(1.+z)**(3.) + Om_L )

# dt/dz, s
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

# time as a function of z
def time_from_z(z):
    #t = 2./3./H0*Om_m**(-1./2.)*(1.+z)**(-3./2.)
    t = 2./3./H0*Om_L**(-1./2.)*np.log( (np.sqrt( Om_L/(1.+z)**3. ) + np.sqrt( Om_m + Om_L/(1.+z)**3. ))/np.sqrt( Om_m) )
    if z > zeq:
        t = ((1.+z)/1.e10)**(-2.)
    return t

# time from zmin to zmax
def time_int_z(zmin, zmax):
    return integrate.quad(lambda a: 1./(a*hubble(a**(-1) - 1.)), 1./(1.+zmax), 1./(1.+zmin))[0] # s
    #return integrate.quad(lambda zz: 1./((1.+zz)*hubble(zz)), zmin, zmax)[0] # s
    #return integrate.quad(lambda zz: dtdz(zz), zmin, zmax)[0] # s

# time from amin to amax
def time_int_a(amin, amax):
    #return integrate.quad(lambda a: dtdz(a**(-1.) - 1.), amin, amax)[0] # s
    return integrate.quad(lambda a: 1./(a*hubble(a**(-1) - 1.)), amin, amax)[0] # s

# redshift as a function of time
def z_from_t(t):
    return 1./(fsolve(lambda aa: time_int_a(0., aa)-t, 1e-3))[0] - 1.

time_int_a = np.vectorize(time_int_a)
z_from_t = np.vectorize(z_from_t)
redshift = np.vectorize(redshift)


# Age of the universe
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
