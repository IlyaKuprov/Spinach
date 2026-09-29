# etc/molecules/zfs_sampling.m

- MATLAB implementation: [etc/molecules/zfs_sampling.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/zfs_sampling.m)

## Purpose

Constructs a discrete quadrature distribution for the zero-field-splitting parameters of gadolinium complexes with DOTA-type ligands in cryogenic water–methanol glasses. The distributions are the parametrisation cited for Figure 5 of the paper below.

## Use

```matlab
[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)
```

All three arguments are required; the source defines no defaults.

- `npoints_d` — finite real integer, at least 5; number of Gauss–Legendre nodes for the reduced D coordinate `D/D0`.
- `npoints_e` — finite real integer, at least 5; number of Gauss–Legendre nodes for `E/D`.
- `tol` — finite, non-negative real scalar. A product-grid point is dropped when its weight is strictly less than this threshold.

## Distribution and outputs

For `D/D0`, nodes span [-2, 2] and their Gauss–Legendre weights are reweighted by equal-width Gaussian peaks centred at -1 and +1. The common standard deviation is `1/(2*sqrt(2*log(2)))`, giving unit full width at half maximum. For `E/D`, nodes span [0, 1/3] and weights are reweighted by `-(E/D - 0.25)^2 + 0.0625`. The two one-dimensional rules are combined as a Cartesian product: returned `D` holds the reduced `D/D0` coordinate, `E` is formed as `D*(E/D)` and therefore represents reduced `E/D0`, not the sampled `E/D` coordinate; `W` holds the product weights. After thresholding, retained weights are renormalised to sum to one. Thus the output grid has at most `npoints_d*npoints_e` points; the actual count depends on `tol`.

The function also opens a figure and plots the two distribution profiles as part of the call. It is therefore not a data-only routine.

## Reference

The distribution parameters are matched to Figure 5 of [the cited study](https://doi.org/10.1007/BF03166762). See also the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=zfs_sampling.m).
