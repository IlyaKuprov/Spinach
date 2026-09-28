# experiments/pseudocon/ipcs.m

- Signature: `[source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b]=ipcs(parameters,npoints,lambda)`

## Purpose

Reconstructs a source cube from measured pseudocontact shifts by solving either the Kuprov equation or Poisson's equation, with contrast and Tikhonov regularisation. The equations and algorithms are described in [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).

## Parameters / inputs

- `parameters.xyz` — N-by-3 coordinates in Å of nuclei with measured PCS values.
- `parameters.expt_pcs` — real column vector of measured PCS values in ppm, with one value per row of `xyz`.
- `parameters.xyz_all` — M-by-3 coordinates in Å of all molecular atoms; used for source confinement and molecule plots.
- `parameters.chi` — real 3-by-3 electron magnetic susceptibility tensor in Å^3; used by the Kuprov equation.
- `parameters.equation` — `'kuprov'` to recover the unpaired-electron probability density, or `'poisson'` to recover the right-hand side of Poisson's equation.
- `parameters.box_cent`, `parameters.box_size` — three-element vectors in Å defining the centre and dimensions of the allowed solution box.
- `parameters.margins` — six-element vector of lower and upper expansions of the measured-coordinate bounding box, ordered by axis, in Å.
- `parameters.confine` — two non-negative radii in Å: density is rejected within the first radius of any atom and accepted within the second.
- `parameters.sharpen` — real scalar weight of the contrast penalty.
- `parameters.plot` — cell array containing any of `'diagnostics'`, `'density'`, `'molecule'`, `'tightzoom'`, and `'box'`.
- `parameters.gpu` — logical scalar selecting GPU processing when true.
- Optional `parameters.guess` — initial source cube, resampled onto the current grid when supplied.
- `npoints` — integer greater than 10 giving the number of grid points along each cube dimension.
- `lambda` — real Tikhonov regularisation weight.

## Outputs

- `source_cube` — reconstructed source term on an [X Y Z] grid.
- `ranges` — Cartesian cube extents [xmin xmax ymin ymax zmin zmax] in Å.
- `pred_pcs` — PCS values in ppm predicted from the reconstructed source at the measured coordinates.
- `err_ls` — least-squares residual error in ppm^2.
- `reg_a` — contrast penalty term.
- `reg_b` — Tikhonov penalty term.

## Method

The routine derives cube bounds from the measured coordinates and supplied margins, applies the solution-box and atom-distance confinement masks, and constructs a tricubic sampling matrix for the measured nuclei. It forms the inverse operator in Fourier space: for the Kuprov equation this uses the rank-2 part of `chi`, while the Poisson option uses the inverse Laplacian. It minimises the PCS least-squares residual together with the contrast term weighted by `sharpen` and the Tikhonov term weighted by `lambda`, using `fmincon` with a trust-region-reflective method and supplied gradient/Hessian-vector calculations. The Kuprov solution is constrained non-negative; the Poisson solution is unbounded. GPU arrays are used when requested and available.

## References

- [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G)
- [Spin Dynamics Wiki: ipcs.m](https://spindynamics.org/wiki/index.php?title=ipcs.m)
