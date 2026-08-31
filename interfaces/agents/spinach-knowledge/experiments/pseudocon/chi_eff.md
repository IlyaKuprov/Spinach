# experiments/pseudocon/chi_eff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/chi_eff.m`
- Signature: `[chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)`
- Total lines: 90

## Purpose

Finds the optimal magnetic susceptibility tensor that a user-supplied paramagnetic centre probability density must have in order to fit the PCS data supplied. Syntax: [chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(source_cube,ranges,nxyz,expt_pcs)`.
- Lines 38-40: Normalise the density; implemented by `normint=(ranges(2)-ranges(1))*(ranges(4)-ranges(3))* (ranges(6)-ranges(5))*sum(sum(sum(source_cube)))/numel(source_cube)`.
- Lines 43-44: Run the optimisation; implemented by `theo_pcs=@(chi)kpcs(source_cube,[chi(1) chi(3) chi(4)`.
- Lines 52-53: Return parameters; implemented by `chi=[params(1) params(3) params(4)`.
- Lines 57-58: Isolate the second spherical rank component of chi; implemented by `[~,~,rank2]=mat2sphten(chi); chi=sphten2mat([],[],rank2)`.
- Lines 60-61: Back-calculate PCS values; implemented by `pred_pcs=kpcs(source_cube,chi,ranges,nxyz,'fft')`.

### Key state/data transformations

- Lines 39-40: computes `normint` using `normint=(ranges(2)-ranges(1))*(ranges(4)-ranges(3))* (ranges(6)-ranges(5))*sum(sum(sum(source_cube)))/numel(source_cube)`.
- Lines 41: computes `source_cube` using `source_cube=source_cube/normint`.
- Lines 44: computes `theo_pcs` using `theo_pcs=@(chi)kpcs(source_cube,[chi(1) chi(3) chi(4)`.
- Lines 47: computes `errfun` using `errfun=@(chi)norm(expt_pcs-theo_pcs(chi),'fro')^2`.
- Lines 48-49: computes `options` using `options=optimset('Display','iter','MaxIter',Inf,'MaxFunEvals',Inf, 'FinDiffType','central','UseParallel',true)`.
- Lines 50: computes `params` using `params=fminunc(errfun,[0.10 0.10 0.01 0.01 0.01],options)`.
- Lines 53: computes `chi` using `chi=[params(1) params(3) params(4)`.
- Lines 58: computes `[~,~,rank2]` using `[~,~,rank2]=mat2sphten(chi); chi=sphten2mat([],[],rank2)`.
- Lines 61: computes `pred_pcs` using `pred_pcs=kpcs(source_cube,chi,ranges,nxyz,'fft')`.

### Local helper functions

- Line 66: `grumble()` — `function grumble(source_cube,ranges,nxyz,expt_pcs)`.
  - Representative operation: `if (~isnumeric(source_cube))||(ndims(source_cube)~=3)`.
  - Representative operation: `error('source_cube must be a three-dimensional numerical array.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ranges()`, `kpcs()`, `chi()`, `theo_pcs()`, `optimset()`, `fminunc()`, `params()`, `mat2sphten()`, `sphten2mat()`, `ndims()`, `any()`, `source_cube()`.
