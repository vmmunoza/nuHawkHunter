# nuHawk (or other name bacano culiao)

Al imperialismo ni un tantico así.

## Scripts

Here is a brief description of the codes included:

- `PBH_Constraints.ipynb`: jupyter notebook to compute the bounds on the PBH abundance.

- `flux_stuff.py`: just some relevant stuff for computing the fluxes.

- `fluxes_tot.py`: computes the neutrino fluxes from the BlackHawk files.

- `event_rate.py`: computes the event rate for the PBH signals and for the backgrounds.

- `cross_sections.py`: includes the relevant cross sections employed.

- `chi2_pro.py`: functions to compute and interpolate the chi2.

## Stuff to do

- Incluir el cálculo de coherent scattering para DARWIN y otros. Implementar la cross section correspondiente en `cross_sections.py` (sacarla de draft y de código de Sam) y el event rate y backgrounds en `event_rate.py`

- Para chi2, hacer bineado de 1 o 2 MeV?. Hasta ahora pillo las energías que tengo tabuladas sin más. Para bineado, habrá que interpolar backgrounds.

- Calcular correctamente los backgrounds NC para JUNO, para lo que hace falta incluir más cross sections (voy a mirar esto yo).

- Para el IBD, uso cross section total en vez de integrar la diferencial, porque esto me da problemas numéricos. Solucionar!

- Ordenar los scrips `flux_stuff.py` y `fluxes_tot.py`, se pueden juntar en uno quizá.

- Para masas >1e15, uso BlackHawk instantáneo. Chequear que sale igual al total, como debería ser al no evaporarse.

- Incluir los bounds correctos para SK, usando los datos (ahora en el código se usan los bakcgrounds como datos, como si fuese un forecast)

- Acabar con el capitalismo.

## Contact

For comments, questions etc. you can contact us at <vmmunoz2@uc.cl> or <pablo.villanueva.domingo@gmail.com>.
