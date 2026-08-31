# kernel/utilities/lcurve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/lcurve.m`
- Signature: `lam_opt=lcurve(lam,err,reg,mode)`
- Total lines: 165

## Purpose

L-curve analysis function. Syntax: lam_opt=lcurve(lam,err,reg,mode)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 47-48: Check consistency; implemented by `grumble(lam,err,reg,mode)`.
- Lines 50-51: Warn about a sweep that is not a trade-off curve, but do not bomb out; implemented by `if any(diff(err)<=0)||any(diff(reg)>=0)`.
- Lines 55-56: Move to logarithmic coordinates; implemented by `log_lam=log10(lam); log_err=log10(err); log_reg=log10(reg)`.
- Lines 58-59: Resample using quintic spline; implemented by `sp_err=spapi(optknt(log_lam,5),log_lam,log_err)`.
- Lines 65-66: Return to linear coordinates; implemented by `err=10.^log_err; reg=10.^log_reg; lam=10.^log_lam`.
- Lines 68-69: Plot the L-curve; implemented by `subplot(1,2,1); plot(err,reg,'b-')`.
- Lines 74-75: Get the derivatives; implemented by `switch mode`.
- Lines 79-80: Derivatives in logarithmic coordinates; implemented by `xp=fdvec(log_err,5,1); xpp=fdvec(log_err,5,2)`.
- Lines 83-84: Plot in logarithmic coordinates; implemented by `set(gca,'xscale','log'); set(gca,'yscale','log')`.
- Lines 88-89: Derivatives in linear coordinates; implemented by `xp=fdvec(err,5,1); xpp=fdvec(err,5,2)`.
- Lines 94-95: Complain and bomb out; implemented by `error('unknown differentiation mode.')`.
- Lines 99-100: Get the signed curvature; implemented by `kappa=(xp.*ypp-yp.*xpp)./((xp.^2+yp.^2).^(3/2))`.
- Lines 102-103: Plot the curvature; implemented by `subplot(1,2,2); plot(lam,kappa)`.
- Lines 109-110: Find the maximum curvature away from the unreliable stencil ends; implemented by `margin=ceil(numel(kappa)/40)`.
- Lines 113-114: A corner at the edge means it is outside the sampled interval; implemented by `if (index<=(2*margin))||(index>=(numel(kappa)-2*margin))`.
- Lines 118-119: Report the optimum point; implemented by `lam_opt=lam(index)`.

### Control flow inferred from the code

- Line 51: conditional branch on `any(diff(err)<=0)||any(diff(reg)>=0)`.
- Line 75: dispatches on `mode`; cases `'log'`, `'linear'`.
- Line 114: conditional branch on `(index<=(2*margin))||(index>=(numel(kappa)-2*margin))`.

### Key state/data transformations

- Lines 56: computes `log_lam` using `log_lam=log10(lam); log_err=log10(err); log_reg=log10(reg)`.
- Lines 59: computes `sp_err` using `sp_err=spapi(optknt(log_lam,5),log_lam,log_err)`.
- Lines 60: computes `sp_reg` using `sp_reg=spapi(optknt(log_lam,5),log_lam,log_reg)`.
- Lines 61: computes `log_err` using `log_err=fnval(linspace(min(log_lam),max(log_lam),1000),sp_err)`.
- Lines 62: computes `log_reg` using `log_reg=fnval(linspace(min(log_lam),max(log_lam),1000),sp_reg)`.
- Lines 66: computes `err` using `err=10.^log_err; reg=10.^log_reg; lam=10.^log_lam`.
- Lines 80: computes `xp` using `xp=fdvec(log_err,5,1); xpp=fdvec(log_err,5,2)`.
- Lines 81: computes `yp` using `yp=fdvec(log_reg,5,1); ypp=fdvec(log_reg,5,2)`.
- Lines 100: computes `kappa` using `kappa=(xp.*ypp-yp.*xpp)./((xp.^2+yp.^2).^(3/2))`.
- Lines 110: computes `margin` using `margin=ceil(numel(kappa)/40)`.
- Lines 111: computes `[~,index]` using `[~,index]=max(kappa((1+margin):(end-margin))); index=index+margin`.
- Lines 119: computes `lam_opt` using `lam_opt=lam(index)`.

### Local helper functions

- Line 126: `grumble()` — `function grumble(lam,err,reg,mode)`.
  - Representative operation: `if (~isnumeric(lam))||(~isreal(lam))||(any(~isfinite(lam)))|| (size(lam,1)~=1)||any(lam<=0)`.
  - Representative operation: `(size(lam,1)~=1)||any(lam<=0)`.

## Parameters / inputs

- lam -row vector of regularisation parameters, must
- be positive and in ascending order
- err -row vector of least squares errors, must be
- positive and increasing with lam
- reg -row vector of regularisation functional values,
- must be positive and decreasing with lam. This
- is the regularisation functional itself, not the
- penalty term of the error functional: when the
- optimiser reports lam*||L*x||^2, divide it by lam
- once before calling this function
- mode -'log' for logarithmic coordinates and 'linear'
- for linear ones; 'log' is recommended

## Outputs

- lam_opt -the regularisation parameter at the point
- of the maximum curvature of the L-curve
- Notes: the corner is the point of the greatest curvature, which for a
- smooth asymmetric bend is not the intersection of the asymptotes;
- the criterion locates the regularisation parameter to within a
- factor of a few and is not convergent in the zero noise limit
- (Vogel, SIAM J. Numer. Anal. 34, 1996). Treat the answer as an
- order of magnitude estimate and inspect the plotted curve.
- The curvature maximum must fall inside the sampled interval; if
- it falls on either end, the corner is outside the range and this
- function refuses to return the endpoint as an answer.
- This function requires the Curve Fitting Toolbox.

## Implementation structure

- L-curve analysis function. Syntax:
- lam_opt=lcurve(lam,err,reg,mode)
- lam -row vector of regularisation parameters, must
- be positive and in ascending order
- err -row vector of least squares errors, must be
- positive and increasing with lam
- reg -row vector of regularisation functional values,
- must be positive and decreasing with lam. This
- is the regularisation functional itself, not the
- penalty term of the error functional: when the
- optimiser reports lam*||L*x||^2, divide it by lam
- once before calling this function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `diff()`, `log10()`, `spapi()`, `optknt()`, `fnval()`, `subplot()`, `kxlabel()`, `kylabel()`, `fdvec()`, `set()`, `kappa()`, `lam()`, `err()`, `reg()`.
