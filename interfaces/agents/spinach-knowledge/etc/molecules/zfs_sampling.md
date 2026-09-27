# etc/molecules/zfs_sampling.m

- Signature: `[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)`

## Purpose

Constructs a discrete probability distribution for the zero-field-splitting (ZFS) parameters of gadolinium complexes with DOTA-type ligands in cryogenic water–methanol glasses. The parameter distributions follow Figure 5 of the cited study.

## Method

Gauss–Legendre quadrature nodes for the reduced parameter D/D₀ span [-2, 2]. Their weights are modulated by an equal-width double-Gaussian profile centered at -1 and +1; the common standard deviation corresponds to unit full width at half maximum. Nodes for E/D span [0, 1/3] and are weighted by the parabolic profile `-(E/D - 0.25)² + 0.0625`. The two node sets are combined as a Cartesian product, with product weights. Points whose joint weight is below `tol` are discarded, and the retained weights are renormalized.

The function also plots the D/D₀ and E/D probability-density profiles used to construct the grid.

## Inputs

- `npoints_d` — number of Gauss–Legendre quadrature points for D/D₀; must be a finite real integer of at least 5.
- `npoints_e` — number of Gauss–Legendre quadrature points for E/D; must be a finite real integer of at least 5.
- `tol` — finite, non-negative real scalar threshold for discarding joint quadrature weights.

## Outputs

- `D` — D/D₀ coordinate at each retained grid point.
- `E` — E/D coordinate at each retained grid point.
- `W` — normalized weight of each retained grid point.

## Reference

Parameter distributions are described in Figure 5 of [the cited study](https://doi.org/10.1007/BF03166762). See also the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=zfs_sampling.m).
