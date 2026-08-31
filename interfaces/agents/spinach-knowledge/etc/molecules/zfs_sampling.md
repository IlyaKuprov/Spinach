# etc/molecules/zfs_sampling.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/zfs_sampling.m`
- Signature: `[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)`
- Total lines: 102

## Purpose

Gadolinium ZFS probability distribution function for DOTA-type ligand complexes in cryogenic water-methanol glasses. The para- meters match those given in Figure 5 of

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(npoints_d,npoints_e,tol)`.
- Lines 45-46: Generate Gauss-Legendre point set for D/D1; implemented by `[X,WX]=gaussleg(-2,2,npoints_d)`.
- Lines 48-49: Get the standard deviation for unit FWHM; implemented by `sigma=1/(2*sqrt(2*log(2)))`.
- Lines 51-52: Refract weights through a double Gaussian; implemented by `WX=WX.*(normpdf(X,-1,sigma)+normpdf(X,+1,sigma)); WX=WX/sum(WX)`.
- Lines 54-55: Plot the double Gaussian; implemented by `kfigure(); scale_figure([1.50 0.75]); subplot(1,2,1)`.
- Lines 60-61: Generate Gauss-legendre point set for E/D; implemented by `[Y,WY]=gaussleg(0,1/3,npoints_e)`.
- Lines 63-64: Refract weights through a quadratic function; implemented by `WY=WY.*(-(Y-0.25).^2+0.0625); WY=WY/sum(WY)`.
- Lines 66-67: Plot the quadratic function; implemented by `subplot(1,2,2); plot(Y,-(Y-0.25).^2+0.0625,'r-')`.
- Lines 71-72: Kron the weights; implemented by `D=kron(X,ones(size(Y)))`.
- Lines 76-77: Ignore small weights; implemented by `D(W<tol)=[]; E(W<tol)=[]; W(W<tol)=[]`.

### Key state/data transformations

- Lines 46: computes `[X,WX]` using `[X,WX]=gaussleg(-2,2,npoints_d)`.
- Lines 49: computes `sigma` using `sigma=1/(2*sqrt(2*log(2)))`.
- Lines 52: computes `WX` using `WX=WX.*(normpdf(X,-1,sigma)+normpdf(X,+1,sigma)); WX=WX/sum(WX)`.
- Lines 61: computes `[Y,WY]` using `[Y,WY]=gaussleg(0,1/3,npoints_e)`.
- Lines 64: computes `WY` using `WY=WY.*(-(Y-0.25).^2+0.0625); WY=WY/sum(WY)`.
- Lines 72: computes `D` using `D=kron(X,ones(size(Y)))`.
- Lines 73: computes `E` using `E=D.*kron(ones(size(X)),Y)`.
- Lines 74: computes `W` using `W=kron(WX,WY); W=W/sum(W)`.
- Lines 77: computes `D(W<tol)` using `D(W<tol)=[]; E(W<tol)=[]; W(W<tol)=[]`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(npoints_d,npoints_e,tol)`.
  - Representative operation: `if (~isnumeric(npoints_d))||(~isreal(npoints_d))|| (~isscalar(npoints_d))||(~isfinite(npoints_d))|| (npoints_d<5)||(mod(npoints_d,1)~=0)`.
  - Representative operation: `(~isscalar(npoints_d))||(~isfinite(npoints_d))|| (npoints_d<5)||(mod(npoints_d,1)~=0)`.

## Syntax

```matlab
[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)
```

## Parameters / inputs

- npoints_d -number of Gauss-Legendre quadrature
- points in D
- npoints_e -number of Gauss-Legendre quadrature
- points in E
- tol -tolerance for integration weights
- below which grid points are dropped

## Outputs

- D -a vector of D values at each integration
- grid point
- E -a vector of E values at each integration
- grid point
- W -a vector of weights for each integration
- grid point
- Notes: the function also creates a figure with the distributi-
- ons it has used for D and E parameters.

## Implementation structure

- Gadolinium ZFS probability distribution function for DOTA-type
- ligand complexes in cryogenic water-methanol glasses. The para-
- meters match those given in Figure 5 of
- [D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)
- npoints_d -number of Gauss-Legendre quadrature
- points in D
- npoints_e -number of Gauss-Legendre quadrature
- points in E
- tol -tolerance for integration weights
- below which grid points are dropped
- D -a vector of D values at each integration
- grid point

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gaussleg()`, `normpdf()`, `kfigure()`, `scale_figure()`, `subplot()`, `ktitle()`, `kxlabel()`, `kylabel()`, `isscalar()`.
