# experiments/pseudocon/chi_eff.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/chi_eff.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=chi_eff.m)

## Purpose

Fits a magnetic-susceptibility tensor to measured pseudocontact shifts (PCS) for a supplied three-dimensional paramagnetic-centre probability-density cube. The spatial density and the nuclear coordinates are inputs; this routine varies the tensor, not the density or centre positions. This is a PCS/tensor fit, not a pulse-transfer simulation.

## Inputs

- `source_cube` is a numeric 3-D array with no negative entries.
- `ranges` is a numeric six-element vector of grid extents in Angstroms, ordered as `[xmin xmax ymin ymax zmin zmax]`.
- `nxyz` is an N-by-3 real array of nuclear coordinates in Angstroms.
- `expt_pcs` is an N-by-1 real column of measured PCS values in ppm; its row count must match `nxyz`.

The code normalises the density by multiplying its mean grid value by the product of the three range extents, then dividing `source_cube` by that quantity. It constructs a symmetric, traceless 3-by-3 tensor from five fitted parameters and minimises the squared Frobenius norm of the difference between measured PCS and `kpcs(...,'fft')` predictions. The initial parameter vector is `[0.10 0.10 0.01 0.01 0.01]`; `fminunc` uses central finite differences and parallel execution settings. After fitting, the code retains the tensor's spherical-rank-2 component and evaluates the predicted PCS again.

## Outputs

- `chi` is the fitted 3-by-3 susceptibility tensor in cubic Angstroms, reconstructed from its rank-2 component.
- `pred_pcs` is a column of predicted PCS values in ppm, one for each row of `nxyz`.

The source header describes the prediction as using an optimised `mxyz` and `chi`, but the implemented signature returns only `chi` and `pred_pcs`; the supplied density and ranges remain fixed inputs.
