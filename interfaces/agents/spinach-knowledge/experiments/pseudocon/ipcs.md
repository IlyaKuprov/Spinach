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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 73-74: Check consistency; implemented by `grumble(parameters,npoints,lambda)`.
- Lines 76-77: Kill stupid ass figure defaults in R2025a and later; implemented by `set(groot,'defaultFigurePosition',[680 458 560 420])`.
- Lines 82-83: Allow GPUs Matlab has not yet seen; implemented by `parallel.gpu.enableCUDAForwardCompatibility(true)`.
- Lines 85-86: Validate GPU option; implemented by `if parameters.gpu`.
- Lines 90-91: Conserve GPU memory; implemented by `if parameters.gpu`.
- Lines 96-97: Set internal scaling parameters; implemented by `switch parameters.equation`.
- Lines 101-102: Pin the density roughly around [0,1]; implemented by `probdens_scaling=1.0e3`.
- Lines 108-109: Pin the density roughly around [0,1]; implemented by `probdens_scaling=1.0e4`.
- Lines 115-116: Set the GPU data types; implemented by `if parameters.gpu`.
- Lines 122-128: Choose appropriate axis ranges; implemented by `ranges=[min(parameters.xyz(:,1))-parameters.margins(1) max(parameters.xyz(:,1))+parameters.margins(2) min(parameters.xyz(:,2))-parameters.margins(3) max(parameters.xyz(:…`.
- Lines 130-131: Report to the user; implemented by `disp(['X axis extents: ' num2str(ranges([1 2])) ' Angstrom.'])`.
- Lines 136-137: Get the nuclear point sampling matrix; implemented by `disp('Building nuclear sampling matrix ')`.
- Lines 140-143: Get the coordinate grid; implemented by `[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),npoints), linspace(ranges(3),ranges(4),npoints), linspace(ranges(5),ranges(6),npoints))`.
- Lines 145-151: Get the solution box and apply distance cutoffs; implemented by `soln_box=(X>(parameters.box_cent(1)-parameters.box_size(1)/2))& (X<(parameters.box_cent(1)+parameters.box_size(1)/2))& (Y>(parameters.box_cent(2)-parameters.box_size(2)/…`.
- Lines 153-154: Confine the solution to the relevant layer; implemented by `if isfield(parameters,'confine')`.
- Lines 168-169: Pre-compute atom connectivity matrix; implemented by `if isfield(parameters,'xyz_all')`.
- Lines 173-174: Index the active space; implemented by `active_point_index=soln_box(:)`.
- Lines 176-177: Report to the user; implemented by `disp([num2str(nnz(active_point_index)) ' points in the probability density'])`.

### Control flow inferred from the code

- Line 86: conditional branch on `parameters.gpu`.
- Line 91: conditional branch on `parameters.gpu`.
- Line 97: dispatches on `parameters.equation`; cases `'kuprov'`, `'poisson'`, `'kuprov'`, `'poisson'`, `'kuprov'`, `'poisson'`.
- Line 116: conditional branch on `parameters.gpu`.
- Line 154: conditional branch on `isfield(parameters,'confine')`.
- Line 158: `parfor` loop over `m=1:size(parameters.xyz_all,1)`.
- Line 169: conditional branch on `isfield(parameters,'xyz_all')`.
- Line 180: conditional branch on `isfield(parameters,'guess')`.
- Line 190: conditional branch on `strcmp(parameters.equation,'kuprov'), guess(guess<0)=0; end`.
- Line 198: dispatches on `parameters.equation`; cases `'kuprov'`, `'poisson'`, `'kuprov'`, `'poisson'`.
- Line 225: conditional branch on `(parameters.sharpen>0)&&(~any(guess))`.
- Line 236: dispatches on `parameters.equation`; cases `'kuprov'`, `'poisson'`.
- Line 260: conditional branch on `parameters.gpu`.

### Key state/data transformations

- Lines 87: computes `parameters.gpu` using `parameters.gpu=logical(gpuDeviceCount)`.
- Lines 92: computes `G` using `G=gpuDevice(); reset(G)`.
- Lines 93: computes `G.CachePolicy` using `G.CachePolicy='minimum'`.
- Lines 102: computes `probdens_scaling` using `probdens_scaling=1.0e3`.
- Lines 103: computes `tikhonov_scaling` using `tikhonov_scaling=0.640/npoints`.
- Lines 104: computes `contrast_scaling` using `contrast_scaling=1.0e3/npoints^3`.
- Lines 117: computes `data_type` using `data_type='gpuArray'`.
- Lines 123-128: computes `ranges` using `ranges=[min(parameters.xyz(:,1))-parameters.margins(1) max(parameters.xyz(:,1))+parameters.margins(2) min(parameters.xyz(:,2))-parameters.margins(3) max(parameters.xyz(:…`.
- Lines 138: computes `P` using `P=interpmat([npoints npoints npoints],ranges,parameters.xyz)`.
- Lines 141-143: computes `[X,Y,Z]` using `[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),npoints), linspace(ranges(3),ranges(4),npoints), linspace(ranges(5),ranges(6),npoints))`.
- Lines 146-151: computes `soln_box` using `soln_box=(X>(parameters.box_cent(1)-parameters.box_size(1)/2))& (X<(parameters.box_cent(1)+parameters.box_size(1)/2))& (Y>(parameters.box_cent(2)-parameters.box_size(2)/…`.
- Lines 156: computes `reject_region` using `reject_region=false(size(soln_box))`.
- Lines 157: computes `accept_region` using `accept_region=false(size(soln_box))`.
- Lines 159-161: computes `distance` using `distance=sqrt((X-parameters.xyz_all(m,1)).^2+ (Y-parameters.xyz_all(m,2)).^2+ (Z-parameters.xyz_all(m,3)).^2)`.
- Lines 170: computes `conmatrix` using `conmatrix=conmat(parameters.xyz_all,1.6)`.
- Lines 174: computes `active_point_index` using `active_point_index=soln_box(:)`.
- Lines 184-186: computes `[X0,Y0,Z0]` using `[X0,Y0,Z0]=ndgrid(linspace(ranges(1),ranges(2),size(parameters.guess,1)), linspace(ranges(3),ranges(4),size(parameters.guess,2)), linspace(ranges(5),ranges(6),size(param…`.
- Lines 187: computes `guess` using `guess=probdens_scaling*interpn(X0,Y0,Z0,parameters.guess,X,Y,Z,'cubic')`.

### Local helper functions

- Line 268: `myobj()` — `function [err,grad_err,x]=myobj(x)`. Embed the active region
  - Representative operation: `den=zeros(npoints^3,1,data_type); den(active_point_index)=x`.
  - Representative operation: `fft_den=fftn(reshape(den,[npoints npoints npoints]))`.
- Line 414: `myhess()` — `function y=myhess(x,v)`. Preallocate the answer
  - Representative operation: `y=zeros(npoints^3,size(v,2),data_type)`.
  - Representative operation: `for n=1:size(v,2)`.
- Line 511: `grumble()` — `function grumble(parameters,npoints,lambda)`.
  - Representative operation: `if verLessThan('matlab','9.0')`.
  - Representative operation: `error('3D PCS reconstruction module requires 64-bit Matlab R2016a or later.')`.

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
