# νHawkHunter
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.0000000.svg)](https://zenodo.org/record/0000000)
[![arXiv](https://img.shields.io/badge/arXiv-2109.XXXX-B31B1B.svg)](http://arxiv.org/abs/2109.XXXX)
![Alt text](https://img.shields.io/pypi/pyversions/python-binance.svg)

blablabla

This code makes use of neutrino fluxes computed with [BlackHawk](https://blackhawk.hepforge.org/) to explore the prospects of detect neutrinos coming from Primordial Black Holes in future neutrino telescopes.

It can also be employed for Diffuse Supernova Neutrino Background or similar studies by replacing the signal fluxes by the proper ones.

But first play some [Axie](https://www.youtube.com/watch?v=Lg5C2EbYueo)

<img src="figures/fluxes_DM.png" width="60%">

In the following we present a brief explanation of the scripts included.

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

- `experiments.py`: includes the definition of the different experiments and the methods to compute the event rate for the PBH signals and for the backgrounds.

- `chi2.py`: functions to compute and interpolate the chi2.


## Stuff to do

- Incluir los bounds correctos para SK, usando cálculo de Sergio

- Acabar con el capitalismo.


## Contact

For comments, questions etc. you can contact us at <vmmunoz2@uc.cl> or <pablo.villanueva.domingo@gmail.com>.
