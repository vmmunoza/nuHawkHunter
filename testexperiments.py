#-------------------------------
# Compute the event rate at different experiments
#-------------------------------

from Source.experiments import *

Mpbhs =  [1e15, 2e15, 4e15]
fpbhs = 1.e-2*np.ones_like(Mpbhs)

list_exps = [HK]

compute_events(Mpbhs, fpbhs, HK, 1, plotevents=1)
