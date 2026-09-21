# experiments/pseudocon/ippcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/ippcs.m`
- Signature: `[mxyz,chi,pred_pcs,s_mxyz,s_chi]=ippcs(nxyz,mguess,expt_pcs)`
- Total lines: 119

## Purpose

Fits the point electron model PCS to the experimental pseudocon- tact shift coordinates and values. Syntax: [exyz,chi,pred_pcs]=ippcs(nxyz,mguess,expt_pcs)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is measured, in Angstroms.
- mguess -initial guess for the unpaired electron coordina-
- tes as [x y z], in Angstroms.
- expt_pcs -pseudocontact shift in ppm at each nucleus.

## Outputs

- mxyz -optimized paramagnetic centre coordinates as [x y z],
- in Angstroms.
- chi -optimized magnetic susceptibility tensor in cubic
- Angstroms.
- pred_pcs -predicted pseudocontact shift at each nucleus with
- the optimized mxyz and chi, ppm.
- s_mxyz -standard deviations of paramagnetic centre
- coordinates as [x y z], in Angstroms.
- s_chi -standard deviations of magnetic susceptibility
- tensor elements in cubic Angstroms.
- Note: a good initial guess for the paramagnetic centre location is
- essential for a successful fit.

## Implementation structure

- Fits the point electron model PCS to the experimental pseudocon-
- tact shift coordinates and values. Syntax:
- [exyz,chi,pred_pcs]=ippcs(nxyz,mguess,expt_pcs)
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is measured, in Angstroms.
- mguess -initial guess for the unpaired electron coordina-
- tes as [x y z], in Angstroms.
- expt_pcs -pseudocontact shift in ppm at each nucleus.
- mxyz -optimized paramagnetic centre coordinates as [x y z],
- in Angstroms.
- chi -optimized magnetic susceptibility tensor in cubic
- Angstroms.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `onCleanup()`, `optimset()`, `ppcs()`, `vec_res_sq()`, `fminunc()`, `jacobianest()`, `sum_res_sq()`, `inv()`, `iscell()`.
