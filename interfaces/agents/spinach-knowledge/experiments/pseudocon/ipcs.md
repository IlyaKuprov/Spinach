# experiments/pseudocon/ipcs.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ipcs.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ipcs.m)

## Purpose

Reconstructs either the unpaired-electron probability density or a Poisson-equation source term from measured pseudocontact shifts (PCS). This is a three-dimensional inverse reconstruction utility, not a pulse-sequence simulator. The model equations and algorithms are documented in [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).

## Inputs and spatial model

Call as `[source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b] = ipcs(parameters,npoints,lambda)`. `parameters.xyz` is an N-by-3 array of PCS-bearing nuclear coordinates in Å, paired row-for-row with the real N-element column `parameters.expt_pcs` in ppm. `parameters.xyz_all` supplies all molecular atom coordinates as M-by-3 Å values. `parameters.chi` is a real symmetric 3-by-3 magnetic-susceptibility tensor in Å³.

`parameters.equation` selects `kuprov` (recover probability density) or `poisson` (recover the Poisson right-hand side). `parameters.box_cent` and `parameters.box_size` specify the centre and three side lengths in Å of the rectangular source-support box. The six `parameters.margins` extend the measured-nucleus bounding box on its lower and upper x, y, and z faces to define the returned grid `ranges`. `parameters.confine` gives two Å radii for excluding source points close to atoms and limiting support to the surrounding molecular region. These choices form a hard voxel mask; values outside the accepted region are fixed at zero.

`npoints` sets each grid dimension, so the reconstructed cube is `npoints` cubed; the source requires an integer greater than 10. `lambda` weights the Tikhonov term. `parameters.sharpen` weights the contrast penalty. An optional `parameters.guess` density cube is interpolated onto the working grid as the starting point. `parameters.plot` is a cell array drawn from `diagnostics`, `density`, `molecule`, `tightzoom`, and `box`; `parameters.gpu` is a logical scalar requesting GPU arrays when a GPU is available.

## Objective and solver

The forward PCS field is evaluated on the regular cube using Fourier-space operators, and `interpmat` samples that field at the measured nuclear coordinates. The least-squares component is the sum of squared differences between those predicted and measured shifts. For `kuprov`, the Fourier multiplier is the Kuprov operator divided by the Laplacian; for `poisson`, it is the inverse Laplacian. The zero-frequency singularity is set to zero. The optimisation adds a contrast penalty proportional to `sharpen` and a Tikhonov penalty based on the squared Laplacian of the source, weighted by `lambda`.

The implementation scales the optimised source variable by `1e3` for `kuprov` and `1e4` for `poisson`; the contrast and Tikhonov scale factors are respectively `1e3/npoints^3` and `0.640/npoints`. The PCS conversion uses `1e6`. `fmincon` uses the trust-region-reflective algorithm, supplied gradient and Hessian-vector product, and `1e-12` function, optimality, and step tolerances. The `kuprov` fit has a zero lower bound (non-negative density); the `poisson` source is unbounded. These are model and optimiser constraints, not evidence of a particular reconstruction quality.

## Outputs

`source_cube` is the reconstructed source on the `npoints`-per-axis grid, and `ranges` is `[xmin xmax ymin ymax zmin zmax]` in Å. `pred_pcs` contains calculated shifts at the input nuclei. `err_ls` is the data-only squared residual in ppm²; `reg_a` and `reg_b` report the contrast and Tikhonov contributions to the objective, respectively.

## References

- [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G)
- [Spin Dynamics Wiki: ipcs.m](https://spindynamics.org/wiki/index.php?title=ipcs.m)