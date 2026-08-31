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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 69-70: Check consistency; implemented by `grumble(nxyz,expt_pcs,ranks,mguess)`.
- Lines 72-73: Count the multipole variables; implemented by `n_mvars=sum(2*nonzeros(ranks)+1)`.
- Lines 75-76: Suppress the method switch warning, restore on exit; implemented by `warn_state=warning('off','optim:fminunc:SwitchingMethod')`.
- Lines 79-81: Set optimisation options; implemented by `options=optimset('Display','iter','MaxIter',inf,'UseParallel',true, 'MaxFunEvals',inf,'FinDiffType','central')`.
- Lines 83-85: Vector of squared residuals; implemented by `vec_res_sq=@(x)(expt_pcs-lpcs(nxyz,x(1:3),ranks, multipack(ranks,[0.5/sqrt(pi) x(9:end)]),x(4:8))).^2`.
- Lines 87-89: Vector of residuals; implemented by `vec_res=@(x)(expt_pcs-lpcs(nxyz,x(1:3),ranks, multipack(ranks,[0.5/sqrt(pi) x(9:end)]),x(4:8)))`.
- Lines 91-92: Sum of squared residuals; implemented by `sum_res_sq=@(x)sum(vec_res_sq(x))`.
- Lines 94-95: Run the fitting to the multipole model; implemented by `p=fminunc(sum_res_sq,[mguess 0.1 0.1 0.1 0.1 0.1 zeros(1,n_mvars)],options)`.
- Lines 97-98: Get metal coordinates; implemented by `mxyz=p(1:3)`.
- Lines 100-101: Get susceptibility; implemented by `chi=[p(4) p(5) p(6)`.
- Lines 105-106: Get multipoles; implemented by `Ilm=multipack(ranks,[0.5/sqrt(pi) p(9:end)])`.
- Lines 108-109: Get the predicted PCS values; implemented by `pred_pcs=lpcs(nxyz,p(1:3),ranks,Ilm,p(4:8))`.
- Lines 111-112: If necessary, compute the statistics; implemented by `if nargout>4`.
- Lines 114-115: Compute Jacobian at the optimal point; implemented by `jac=jacobianest(vec_res,p)`.
- Lines 117-118: Get the Studentized residual; implemented by `sdr=sqrt(sum_res_sq(p)/(numel(expt_pcs)-n_mvars-8))`.
- Lines 120-121: Get the standard deviations; implemented by `sp=sqrt(diag((sdr^2)*inv(jac'*jac)))'`.
- Lines 123-124: Get metal SD; implemented by `s_mxyz=sp(1:3)`.
- Lines 126-127: Get susceptibility SD; implemented by `s_chi=[sp(4) sp(5) sp(6)`.

### Control flow inferred from the code

- Line 112: conditional branch on `nargout>4`.

### Key state/data transformations

- Lines 73: computes `n_mvars` using `n_mvars=sum(2*nonzeros(ranks)+1)`.
- Lines 76: computes `warn_state` using `warn_state=warning('off','optim:fminunc:SwitchingMethod')`.
- Lines 77: computes `cleanup_obj` using `cleanup_obj=onCleanup(@()warning(warn_state))`.
- Lines 80-81: computes `options` using `options=optimset('Display','iter','MaxIter',inf,'UseParallel',true, 'MaxFunEvals',inf,'FinDiffType','central')`.
- Lines 84-85: computes `vec_res_sq` using `vec_res_sq=@(x)(expt_pcs-lpcs(nxyz,x(1:3),ranks, multipack(ranks,[0.5/sqrt(pi) x(9:end)]),x(4:8))).^2`.
- Lines 88-89: computes `vec_res` using `vec_res=@(x)(expt_pcs-lpcs(nxyz,x(1:3),ranks, multipack(ranks,[0.5/sqrt(pi) x(9:end)]),x(4:8)))`.
- Lines 92: computes `sum_res_sq` using `sum_res_sq=@(x)sum(vec_res_sq(x))`.
- Lines 95: computes `p` using `p=fminunc(sum_res_sq,[mguess 0.1 0.1 0.1 0.1 0.1 zeros(1,n_mvars)],options)`.
- Lines 98: computes `mxyz` using `mxyz=p(1:3)`.
- Lines 101: computes `chi` using `chi=[p(4) p(5) p(6)`.
- Lines 106: computes `Ilm` using `Ilm=multipack(ranks,[0.5/sqrt(pi) p(9:end)])`.
- Lines 109: computes `pred_pcs` using `pred_pcs=lpcs(nxyz,p(1:3),ranks,Ilm,p(4:8))`.
- Lines 115: computes `jac` using `jac=jacobianest(vec_res,p)`.
- Lines 118: computes `sdr` using `sdr=sqrt(sum_res_sq(p)/(numel(expt_pcs)-n_mvars-8))`.
- Lines 121: computes `sp` using `sp=sqrt(diag((sdr^2)*inv(jac'*jac)))'`.
- Lines 124: computes `s_mxyz` using `s_mxyz=sp(1:3)`.
- Lines 127: computes `s_chi` using `s_chi=[sp(4) sp(5) sp(6)`.
- Lines 132: computes `s_Ilm` using `s_Ilm=multipack(ranks,[0 sp(9:end)])`.

### Local helper functions

- Line 139: `grumble()` — `function grumble(nxyz,expt_pcs,L,mguess)`.
  - Representative operation: `if (~isnumeric(nxyz))||(~isreal(nxyz))||(size(nxyz,2)~=3)`.
  - Representative operation: `error('nxyz must be an Nx3 array of atomic coordinates.')`.

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
