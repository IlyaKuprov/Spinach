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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(nxyz,mguess,expt_pcs)`.
- Lines 46-47: Suppress the method switch warning, restore on exit; implemented by `warn_state=warning('off','optim:fminunc:SwitchingMethod')`.
- Lines 50-52: Set optimisation options; implemented by `options=optimset('Display','iter','MaxIter',inf,'UseParallel',true, 'MaxFunEvals',inf,'FinDiffType','central')`.
- Lines 54-55: Vector of squared residuals; implemented by `vec_res_sq=@(x)(expt_pcs-ppcs(nxyz,x(1:3),x(4:8))).^2`.
- Lines 57-58: Vector of residuals; implemented by `vec_res=@(x)(expt_pcs-ppcs(nxyz,x(1:3),x(4:8)))`.
- Lines 60-61: Sum of squared residuals; implemented by `sum_res_sq=@(x)sum(vec_res_sq(x))`.
- Lines 63-64: Fit the point model; implemented by `p=fminunc(sum_res_sq,[mguess 0.1 0.1 0.1 0.1 0.1],options)`.
- Lines 66-67: Get the predicted PCS values; implemented by `pred_pcs=ppcs(nxyz,p(1:3),p(4:8))`.
- Lines 69-70: Compute the Jacobian at the optimal point; implemented by `jac=jacobianest(vec_res,p)`.
- Lines 72-73: Get the Studentized residual; implemented by `sdr=sqrt(sum_res_sq(p)/(numel(expt_pcs)-8))`.
- Lines 75-76: Get the standard deviations; implemented by `sp=sqrt(diag((sdr^2)*inv(jac'*jac)))'`.
- Lines 78-79: Return parameters and standard deviations; implemented by `mxyz=p(1:3); s_mxyz=sp(1:3)`.

### Key state/data transformations

- Lines 47: computes `warn_state` using `warn_state=warning('off','optim:fminunc:SwitchingMethod')`.
- Lines 48: computes `cleanup_obj` using `cleanup_obj=onCleanup(@()warning(warn_state))`.
- Lines 51-52: computes `options` using `options=optimset('Display','iter','MaxIter',inf,'UseParallel',true, 'MaxFunEvals',inf,'FinDiffType','central')`.
- Lines 55: computes `vec_res_sq` using `vec_res_sq=@(x)(expt_pcs-ppcs(nxyz,x(1:3),x(4:8))).^2`.
- Lines 58: computes `vec_res` using `vec_res=@(x)(expt_pcs-ppcs(nxyz,x(1:3),x(4:8)))`.
- Lines 61: computes `sum_res_sq` using `sum_res_sq=@(x)sum(vec_res_sq(x))`.
- Lines 64: computes `p` using `p=fminunc(sum_res_sq,[mguess 0.1 0.1 0.1 0.1 0.1],options)`.
- Lines 67: computes `pred_pcs` using `pred_pcs=ppcs(nxyz,p(1:3),p(4:8))`.
- Lines 70: computes `jac` using `jac=jacobianest(vec_res,p)`.
- Lines 73: computes `sdr` using `sdr=sqrt(sum_res_sq(p)/(numel(expt_pcs)-8))`.
- Lines 76: computes `sp` using `sp=sqrt(diag((sdr^2)*inv(jac'*jac)))'`.
- Lines 79: computes `mxyz` using `mxyz=p(1:3); s_mxyz=sp(1:3)`.
- Lines 80: computes `chi` using `chi=[p(4) p(5) p(6)`.
- Lines 83: computes `s_chi` using `s_chi=[sp(4) sp(5) sp(6)`.

### Local helper functions

- Line 90: `grumble()` — `function grumble(nxyz,mguess,expt_pcs)`.
  - Representative operation: `if iscell(nxyz)`.
  - Representative operation: `for n=1:numel(nxyz)`.

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
