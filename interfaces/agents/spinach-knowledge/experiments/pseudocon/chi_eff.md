# experiments/pseudocon/chi_eff.m

- Signature: `[chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)`

## Purpose

Finds the optimal magnetic susceptibility tensor that a user-supplied paramagnetic centre probability density must have in order to fit the PCS data supplied. Syntax: [chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

## Parameters / inputs

- source_cube -paramagnetic centre probability density cube
- ranges -a six-element vector giving the extents of the
- probability density cube in Angstroms as
- [xmin xmax ymin ymax zmin zmax]
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is measured, in Angstroms.
- expt_pcs -pseudocontact shift in ppm at each nucleus.

## Outputs

- chi -optimised magnetic susceptibility tensor in cubic
- Angstroms.
- pred_pcs -predicted pseudocontact shift at each nucleus with
- the optimised mxyz and chi, ppm.

## Implementation structure

- Finds the optimal magnetic susceptibility tensor that a user-supplied
- paramagnetic centre probability density must have in order to fit the
- PCS data supplied. Syntax:
- [chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)
- source_cube -paramagnetic centre probability density cube
- ranges -a six-element vector giving the extents of the
- probability density cube in Angstroms as
- [xmin xmax ymin ymax zmin zmax]
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is measured, in Angstroms.
- expt_pcs -pseudocontact shift in ppm at each nucleus.
- chi -optimised magnetic susceptibility tensor in cubic
