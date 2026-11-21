# cce_nv_diamond

The inset of Onizhuk and Galli, [Figure 3(a)](https://doi.org/10.1002/adts.202100254), distinguishes Hahn-echo modulation by isolated nuclei from the effect of coupled nuclear pairs. The example calculates both orders, not the longer-time ensemble-averaged decay in the main panel. The field is 500 G along [111], the electron transition is 0/−1, and the nuclear bath is initially maximally mixed.

The input file contains only positions: the seed-5 realisation of the authors’ [diamond ensemble generator](https://github.com/MICCoMpy/PyCCE/blob/5ed9c197f7dd1267099c58c3bc7c4745764259e9/examples/nv_ensemble/nv_ensemble.py), restricted to 40 Å. The 493 carbon-13 nuclei provide 422 pairs at the 6 Å connectivity cutoff. The projected electron hyperfine and the full carbon–carbon dipolar tensors are computed from those positions. The plotted line is Eq. (29) evaluated from the singleton Hamiltonians, not digitised figure data.

`[answer,diagnostics]=cce_nv_diamond()` returns the 501-point CCE2 echo; `diagnostics.cce1` is the CCE1 result, `diagnostics.isolated` the independent analytic product, and `diagnostics.times` the total echo times in seconds.
