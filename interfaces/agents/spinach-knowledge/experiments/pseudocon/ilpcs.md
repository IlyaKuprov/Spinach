# experiments/pseudocon/ilpcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/ilpcs.m`
- Signature: `[mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=ilpcs(nxyz,expt_pcs,ranks,mguess)`
- Total lines: 161

## Purpose

Fits experimental PCS data using the distributed paramagnetc centre model described in

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
[mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=...
ilpcs(nxyz,expt_pcs,ranks,mguess)
```

## Parameters / inputs

- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is to be evaluated, in Angstroms.
- expt_pcs -a column vector of experimental pseudocontact shifts
- in ppm
- ranks -row of multipole expansion ranks to be used in the
- fitting procedure
- mguess -guess value for the paramagnetic centre position,
- a three-element vector in Angstrom
- Output:
- mxyz -optimized paramagnetic centre coordinates as [x y z],
- in Angstroms.
- chi -optimized magnetic susceptibility tensor in cubic
- Angstroms.
- Ilm -{[],[]} cell array of numbers corresponding to the
- multipole moments defined in the paper cited above:
- for L=0, Ilm=N/2/sqrt(pi)
- for L=1, Ilm=[real(I11) I10 imag(I11)]
- for L=2, Ilm=[real(I22) real(I21) I20 imag(I21) imag(I22)]
- et cetera.
- pred_pcs -predicted pseudocontact shift (in ppm) at each of
- the nuclei.
- chi -optimized magnetic susceptibility tensor in cubic
- Angstroms.
- s_mxyz -standard deviations of paramagnetic centre
- coordinates as [x y z], in Angstroms.
- s_chi -standard deviations of magnetic susceptibility
- tensor elements in cubic Angstroms.
- s_Ilm -standard deviations of the multipole moments, arranged
- in the same order as the moments themselves.
- Note: a good initial guess for the paramagnetic centre location is
- essential for a successful fit.

## Implementation structure

- Fits experimental PCS data using the distributed paramagnetc centre
- model described in
- [mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=...
- ilpcs(nxyz,expt_pcs,ranks,mguess)
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is to be evaluated, in Angstroms.
- expt_pcs -a column vector of experimental pseudocontact shifts
- in ppm
- ranks -row of multipole expansion ranks to be used in the
- fitting procedure
- mguess -guess value for the paramagnetic centre position,
- a three-element vector in Angstrom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nonzeros()`, `onCleanup()`, `optimset()`, `lpcs()`, `multipack()`, `vec_res_sq()`, `fminunc()`, `jacobianest()`, `sum_res_sq()`, `inv()`.
