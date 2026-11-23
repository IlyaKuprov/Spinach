# cce_nv_diamond

The inset of Onizhuk and Galli, [Figure 3(a)](https://doi.org/10.1002/adts.202100254), distinguishes Hahn-echo modulation by isolated nuclei from the effect of coupled nuclear pairs. The example calculates both orders, not the longer-time ensemble-averaged decay in the main panel. The field is 500 G along [111], the electron transition is 0/−1, and the nuclear bath is initially maximally mixed.

The input file contains only positions: the seed-5 realisation of the authors’ [diamond ensemble generator](https://github.com/MICCoMpy/PyCCE/blob/5ed9c197f7dd1267099c58c3bc7c4745764259e9/examples/nv_ensemble/nv_ensemble.py), restricted to 40 Å. The 493 carbon-13 nuclei provide 422 pairs at the 6 Å connectivity cutoff. The projected electron hyperfine is computed from those positions in the example. The carbon coordinates go to `create` with `sys.tols.prox_cutoff=6`, which builds the carbon–carbon dipolar tensors, and the cluster graph follows the retained couplings; the electron is given no coordinates. The plotted line is Eq. (29) evaluated from the singleton Hamiltonians, not digitised figure data.

Run `cce_nv_diamond()` to plot the real parts of the Spinach CCE1 and CCE2 echoes on alternating points of the 501-point time grid, together with the independent analytic product.
