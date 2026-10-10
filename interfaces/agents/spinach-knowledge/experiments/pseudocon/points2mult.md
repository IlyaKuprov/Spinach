# experiments/pseudocon/points2mult.m

- Source: [experiments/pseudocon/points2mult.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/points2mult.m)
- Wiki: [points2mult.m](https://spindynamics.org/wiki/index.php?title=points2mult.m)
- Signature: `Ilm=points2mult(xyz,mxyz,rho,L,method)`

## Purpose

Computes spherical multipole moments of a supplied spin-population distribution about a paramagnetic centre. The source associates the moments with Equation 32 of DOI [10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D).

## Inputs and coordinate convention

- `xyz` is an `N×3` array of sample coordinates `[x y z]`, in Angstroms; each row pairs with one entry of `rho`.
- `mxyz` is one real `[x y z]` centre coordinate, also in Angstroms. The implementation subtracts it from every row of `xyz`, so the spherical expansion is centred at the paramagnetic site.
- `rho` is a real column vector of sample spin populations/densities. It must have the same number of rows as `xyz`; its units and normalisation are not specified by the source.
- `L` is a real numeric vector of requested spherical-harmonic ranks. The output has one cell per entry of `L`.
- `method` selects how the discrete sum is interpreted: `'points'` for point populations such as Mulliken populations, or `'grid'` for values from a vectorised uniform cubic grid made with `ndgrid`.

## Calculation and output

For each requested rank `l`, the routine converts the recentred coordinates to spherical coordinates and sums `rho .* r.^l .* spher_harmon(l,m,theta,phi)` over samples for `m=0,...,l`. Each cell `Ilm{n}` stores `2*l+1` real components in the slots corresponding to `m=-l,...,l`: the zero component is the real-valued sum for `m=0`; for positive `m`, the positive slot receives the real part and the negative slot receives minus the imaginary part of the complex harmonic sum. This is a coefficient encoding, not a returned vector of complex harmonics.

With `method='grid'`, the result is multiplied by `dx*dy*dz`, with each spacing computed from the unique coordinate values on that axis as `(max-min)/(number of unique values-1)`. `'points'` applies no volume factor. The grid branch therefore relies on the supplied points forming a uniformly spaced Cartesian grid; the routine does not check that condition. The method selector is not checked by its input validator.

The validator requires real numeric `xyz` with three columns, real numeric column `rho` of matching length, real numeric row `mxyz` with three elements, and real numeric `L`. It does not check density normalisation or grid uniformity.
