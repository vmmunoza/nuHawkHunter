# nuHawk (or other name bacano culiao)

Al imperialismo ni un tantico así.

It makes use of neutrino fluxes computed with [BlackHawk](https://blackhawk.hepforge.org/).
It can be employed for Diffuse Supernova Neutrino Background studies by replacing the signal fluxes by the proper ones.

<img src="figures/fluxes_DM.png" style="width:60%">


## Autoblackhawk

We include the driver `Autoblackhawk.py` to run the BlackHawk code for several masses and save the neutrino files. It must be placed in the main directory of BlackHawk, and modifies the parameter file and runs the code automatically according to the PBH masses indicated.


## Jupyter Notebooks

There are several Jupyter notebooks to plot relevant quantities and compute bounds.

- `plot_neutrino_spectrum.ipynb`: plot spectrum rates from BlackHawk.

- `plot_fluxes.ipynb`: computes the neutrino fluxes from the BlackHawk files.

- `plot_events.ipynb`: plots the event rate for the PBH signals and for the backgrounds.

- `PBH_constraints.ipynb`: jupyter notebook to compute the bounds on the PBH abundance.


## Source

Here is a brief description of the scripts included in `Source`, where the relevant computations are defined:

- `constants.py`: definition of relevant constants and initialization of several packages.

- `cosmo.py`: includes some useful cosmological functions.

- `evaporation.py`: includes functions related with the Hawking evaporation.

- `flux_stuff.py`: just some relevant stuff for computing the fluxes.

- `cross_sections.py`: includes the relevant cross sections employed.

- `event_rate.py`: computes the event rate for the PBH signals and for the backgrounds.

- `chi2.py`: functions to compute and interpolate the chi2.


## Stuff to do

- Incluir el cálculo de coherent scattering para DARWIN y otros. Implementar la cross section correspondiente en `cross_sections.py` (sacarla de draft y de código de Sam) y el event rate y backgrounds en `event_rate.py`

- Calcular correctamente los backgrounds NC para JUNO, para lo que hace falta incluir más cross sections.

- Para el IBD, se usa la cross section total en vez de integrar la diferencial, porque esto me da problemas numéricos. Buena aproximación? Solucionar!

- Incluir los bounds correctos para SK, usando cálculo de Sergio

- Acabar con el capitalismo.


## Contact

For comments, questions etc. you can contact us at <vmmunoz2@uc.cl> or <pablo.villanueva.domingo@gmail.com>.
