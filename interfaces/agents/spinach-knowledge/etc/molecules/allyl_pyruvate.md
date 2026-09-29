# etc/molecules/allyl_pyruvate.m

- MATLAB implementation: [etc/molecules/allyl_pyruvate.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/allyl_pyruvate.m)

## Call

`[sys,inter] = allyl_pyruvate(spins)`, for example:

```MATLAB
[sys,inter] = allyl_pyruvate({'1H','13C'});
```

`spins` must be a cell array of character vectors naming isotopes. The source defines six `13C` sites (`C1`–`C6`) and eight `1H` sites (`Ha`, `Hb`, `Hc`, `Hd1`, `Hd2`, `He1`, `He2`, `He3`); the argument selects which isotope rows and matching interaction subarrays to retain. It does not change the underlying parameter set. Each returned structure is a Spinach spin-system (`sys`) or interaction (`inter`) description.

## Parameter provenance and limits

Isotropic shifts and scalar J couplings are described in the source as obtained by spectral fitting; chemical-shift-anisotropy matrices and Cartesian coordinates are from DFT. The function supplies 13C–1H and 1H–1H scalar couplings, but explicitly supplies no 13C–13C couplings. The source therefore identifies this as a natural-abundance 13C simulation system, not a parameter-complete enriched-13C model. It provides no basis-set output.

The arrays contain site-specific scalar shifts, CSA matrices, coordinates, and couplings; selecting isotopes prunes these arrays consistently. The source does not annotate coordinate or interaction units, so do not infer units from the numeric literals alone. Spectral fitting and DFT are provenance labels, not details of the fitting protocol or computational settings; those are not specified here.

Source: [allyl_pyruvate.m in Spinach](https://spindynamics.org/wiki/index.php?title=allyl_pyruvate.m). The source credits `a.acharya@soton.ac.uk` and `ilya.kuprov@weizmann.ac.il`.
