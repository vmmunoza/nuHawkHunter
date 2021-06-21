from Source.event_rate import *

Mpbhs =  [2e15]
#Mpbhs = [1e15, 2e15]
fpbhs = 4.e-3*np.ones_like(Mpbhs)
#fpbhs = np.ones_like(Mpbhs)

cols = ["r", "purple", "b", "g", "orange"]

def plot_events(Mpbhs, fpbhs, exp, as_DM):

    compute_events(Mpbhs, fpbhs, exp, as_DM, plotevents=1, binevents=0)

    #print(E_nu_min_CE(5.e-3, 40*m_p), E_nu_min_CE(1.e-2, 40*m_p), E_nu_min_CE(1.e-1, 40*m_p))
    #print("DARWIN:", E_nu_min_CE(5.e-3, 132*m_p), E_nu_min_CE(5.e-2, 132*m_p))
    #print("ARGO:", E_nu_min_CE(5.e-3, 40*m_p), E_nu_min_CE(5.e-2, 40*m_p))

    EE, evs = np.loadtxt("darwin_paper.csv",unpack=True,delimiter=",")
    plt.plot(EE*1.e-3, evs, "g:", label="PBH DARWIN paper")

    EE, evs = np.loadtxt("CEvENS_bkg.csv",unpack=True,delimiter=",")
    plt.plot(EE*1.e-3, evs, "k--", label="Bkg DARWIN paper")

    #plt.ylim(1.e-6, 1.e4)
    #plt.xlim(1.e-4, 1.e-1)
    plt.ylim(1.e-2, 1.e2)
    plt.xlim(3.e-3, 6.e-2)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(loc="lower right")
    plt.xlabel('$E{\\rm \,\, [MeV]}$')
    #plt.ylabel('${\\rm d}N/d E \,\, [{\\rm MeV}^{-1}{\\rm yr}^{-1}{}{\\rm ton}^{-1}]$')
    plt.ylabel('${\\rm d}N/d E \,\, [{\\rm keV}^{-1}]$')
    plt.title(exp)
    plt.tick_params(axis='both', which='both', top=True, right=True, direction="in")
    plt.grid(which="major",linestyle=":",linewidth=1)
    plt.savefig("figures/test_events_"+exp+".png", bbox_inches='tight', dpi=300)
    plt.show()
    plt.gcf().clear()

exp = "ARGO"
exp = "DARWIN"
#exp ="HK"
plot_events(Mpbhs, fpbhs, exp, 1)
