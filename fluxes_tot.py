from flux_stuff import *


# 1 for plotting, 0 otherwise
plot_fluxes = 1
# 1 for PBH masses which form DM, 0 otherwise
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
