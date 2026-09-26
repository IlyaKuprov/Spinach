# experiments/pseudocon/ippcs.m

- Signature: `[mxyz,chi,pred_pcs,s_mxyz,s_chi]=ippcs(nxyz,mguess,expt_pcs)`

## Purpose

Fits the point electron model PCS to the experimental pseudocon- tact shift coordinates and values. Syntax: [exyz,chi,pred_pcs]=ippcs(nxyz,mguess,expt_pcs)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

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
