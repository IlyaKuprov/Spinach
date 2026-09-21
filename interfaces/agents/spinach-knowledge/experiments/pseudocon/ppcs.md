# experiments/pseudocon/ppcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/ppcs.m`
- Signature: `pcs=ppcs(nxyz,sxyz,chi)`
- Total lines: 80

## Purpose

Computes pseudocontact shift from a point electron centre at the nuclear coordinates supplied. Syntax: pred_pcs=ppcs(nxyz,mxyz,chi)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- chi -magnetic susceptibility tensor in cubic Angstroms
- as a 3x3 matrix, or its five unique components
- ordered as
- [chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is to be evaluated, in Angstroms.
- sxyz -susceptibility centre coordinates as [x y z], in
- Angstroms.
- Output:
- pcs -predicted pseudocontact shift (in ppm) at each of
- the nuclei.

## Implementation structure

- Computes pseudocontact shift from a point electron centre at the
- nuclear coordinates supplied. Syntax:
- pred_pcs=ppcs(nxyz,mxyz,chi)
- chi -magnetic susceptibility tensor in cubic Angstroms
- as a 3x3 matrix, or its five unique components
- ordered as
- [chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is to be evaluated, in Angstroms.
- sxyz -susceptibility centre coordinates as [x y z], in
- Angstroms.
- Output:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `chi()`, `grumble()`, `xyz2sph()`, `nxyz()`, `qform2sph()`, `chi2()`, `spher_harmon()`.
