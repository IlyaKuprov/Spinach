# experiments/pseudocon/ilpcs.m

- Signature: `[mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=ilpcs(nxyz,expt_pcs,ranks,mguess)`
- MATLAB source: [`experiments/pseudocon/ilpcs.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ilpcs.m)

## Purpose

Fits measured PCS values to a distributed paramagnetic-centre multipole model. The model is the one described in [10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D). This is a numerical parameter fit, not a pulse-sequence routine.

## Inputs and outputs

- `nxyz` — N-by-3 array of nuclear coordinates, in ångströms (Å), at the PCS observation sites.
- `expt_pcs` — real N-by-1 column of experimental PCS values, in ppm; N must match the number of rows of `nxyz`.
- `ranks` — row vector of unique non-negative integer multipole ranks, beginning with 0.
- `mguess` — 1-by-3 initial paramagnetic-centre coordinate in Å. The source notes that a good initial centre guess is essential for a successful fit.
- `mxyz` — fitted paramagnetic-centre coordinate, in Å.
- `chi` — fitted magnetic-susceptibility tensor, in cubic ångströms (Å³).
- `Ilm` — cell-array representation of the fitted multipole moments, in the ordering defined by the cited model.
- `pred_pcs` — PCS values predicted at the input coordinates, in ppm.
- `s_mxyz`, `s_chi`, and `s_Ilm` — estimated standard deviations for the fitted centre, susceptibility parameters, and multipole moments. These are calculated only when more than four outputs are requested.

The coordinates are used as Cartesian components in the same frame throughout the fit; this routine does not rotate or reorient them.

## Fit parameterisation and assumptions

The objective is the unweighted sum of squared differences between `expt_pcs` and values from `lpcs`. The minimiser is MATLAB `fminunc`, with central finite differences and parallel evaluation enabled. It starts from `mguess`, five susceptibility parameters initialised to 0.1, and zero initial values for the adjustable multipoles. The rank-0 moment is fixed at `0.5/sqrt(pi)`; the remaining requested multipole components are fitted.

The susceptibility tensor is parameterised by five values as a symmetric traceless matrix: its third diagonal element is set to minus the sum of the first two. There are no explicit bounds in the `fminunc` call. The source estimates parameter standard deviations from a numerically estimated residual Jacobian and a residual-variance factor with `N - n_mvars - 8` degrees of freedom, where `n_mvars` counts adjustable multipole components. The initial centre estimate and this model parameterisation therefore matter to interpreting the fit.

Input checks enforce the coordinate, PCS, initial-guess, and rank shapes/types described above. The function itself does not report an experimental fit quality beyond its fitted and predicted values and the optional parameter standard deviations.

## References

- [10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D)
- [Spin Dynamics Wiki: ilpcs.m](https://spindynamics.org/wiki/index.php?title=ilpcs.m)
