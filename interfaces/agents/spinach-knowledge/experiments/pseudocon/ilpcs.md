# experiments/pseudocon/ilpcs.m

- Signature: `[mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=ilpcs(nxyz,expt_pcs,ranks,mguess)`

## Purpose

Fits experimental pseudocontact shifts (PCS) to the distributed paramagnetic-centre multipole model described in [10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D).

## Parameters / inputs

- `nxyz` — N-by-3 nuclear coordinates in Å at which the PCS values are evaluated.
- `expt_pcs` — real column vector of experimental PCS values in ppm, with one value per coordinate row.
- `ranks` — row vector of unique non-negative integer multipole ranks, starting with 0. The rank-0 moment is fixed by normalisation, not fitted.
- `mguess` — real three-element row vector giving the initial paramagnetic-centre coordinates in Å. A good initial guess is important for successful fitting.

## Outputs

- `mxyz` — fitted paramagnetic-centre coordinates [x y z] in Å.
- `chi` — fitted magnetic-susceptibility tensor in Å^3.
- `Ilm` — cell array of fitted multipole moments corresponding to `ranks`; the moment components are packed in the ordering defined in the cited paper (for example, rank 0 is `N/2/sqrt(pi)`).
- `pred_pcs` — fitted PCS values in ppm at the input coordinates.
- `s_mxyz`, `s_chi`, `s_Ilm` — estimated standard deviations for the centre coordinates, susceptibility tensor elements, and multipole moments, respectively. These are calculated when one of these later outputs is requested.

## Method

The routine minimises the sum of squared residuals between `expt_pcs` and predictions from `lpcs`, varying the centre coordinates, five independent elements of the traceless susceptibility tensor, and the non-fixed multipole components. It uses `fminunc` with central finite differences and parallel evaluation enabled. If uncertainty outputs are requested, it estimates the residual Jacobian at the fitted point and uses it to calculate parameter standard deviations.

## References

- [10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D)
- [Spin Dynamics Wiki: ilpcs.m](https://spindynamics.org/wiki/index.php?title=ilpcs.m)
