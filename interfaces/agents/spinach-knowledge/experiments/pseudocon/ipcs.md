# experiments/pseudocon/ipcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/ipcs.m`
- Signature: `[source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b]=ipcs(parameters,npoints,lambda)`
- Total lines: 602

## Purpose

Solves the inverse problem for pseudocontact shift by recovering the source term in the Kuprov equation using Tikhonov regularisation procedure. Syntax: [source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b]=... ipcs(parameters,nxyz,expt_pcs,chi,npoints,... lambda,margins,box_centre,box_size)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `myobj()`, `myhess()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- nxyz -nuclear coordinates as [x y z] with multiple rows
- at which PCS has been measured, in Angstroms
- expt_pcs -pseudocontact shift in ppm at each nucleus
- chi -electron magnetic susceptibility tensor, in units
- of Angstrom^3
- npoints -number of points in each dimension of the source
- cube, a positive integer greater than 10
- lambda -regularisation parameters, the first element is
- the coefficient in front of the contrast term
- and the second element is the coefficient in
- front of the Tikhonov term
- margins -a six-element vector specifying margins to take
- around the bounding box of the nuclear coordina-
- tes supplied, to account for the possibility that
- some unpaired electron may be located on the pe-
- riphery and require adequare margins
- box_centre -a three-element vector in Angstrom specifying
- the centre of the solution box
- box_size -a three-element vector in Angstrom specifying
- the size of the solution box
- equation -'poisson' to recover the right hand side of the
- Poisson's equation, 'kuprov' to recover the
- unpaired electron probability density
- gpu -set to 1 to enable GPU processing (much faster)

## Outputs

- source_cube -source term cube with dimensions ordered as
- [X Y Z]
- ranges -Cartesian axis extents for the source cube as
- [xmin xmax ymin ymax zmin zmax] in Angstroms
- pred_pcs -pseudocontact shifts produced by the source
- cube returned in the first parameter
- ls_err -least squares error in ppm^2
- reg_a -contrast penalty term
- reg_b -Tikhonov penalty term
- Note: for further information on the equations and algorithms used
- in this function see http://dx.doi.org/10.1039/C4CP03106G

## Implementation structure

- Solves the inverse problem for pseudocontact shift by recovering the
- source term in the Kuprov equation using Tikhonov regularisation
- procedure. Syntax:
- [source_cube,ranges,pred_pcs,err_ls,reg_a,reg_b]=...
- ipcs(parameters,nxyz,expt_pcs,chi,npoints,...
- lambda,margins,box_centre,box_size)
- nxyz -nuclear coordinates as [x y z] with multiple rows
- at which PCS has been measured, in Angstroms
- expt_pcs -pseudocontact shift in ppm at each nucleus
- chi -electron magnetic susceptibility tensor, in units
- of Angstrom^3
- npoints -number of points in each dimension of the source

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `set()`, `logical()`, `gpuDevice()`, `reset()`, `num2str()`, `ranges()`, `interpmat()`, `isfield()`, `false()`, `conmat()`, `soln_box()`, `nnz()`, `interpn()`, `strcmp()`, `guess()`.
