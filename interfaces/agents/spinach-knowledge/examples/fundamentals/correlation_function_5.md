# examples/fundamentals/correlation_function_5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/correlation_function_5.m`
- Signature: `correlation_function_5()`
- Total lines: 50

## Purpose

Computes the following rotational correlation function G(k,m,p,q)=<R(k,m)*R(p,q)> where R is the 3D Cartesian rotation matrix, using the Monte-Carlo method. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Set testing parameters; implemented by `sigma_iso=0.2; k=2; m=3; p=2; q=3`.
- Lines 17-18: Set number of points; implemented by `npoints=1e6; nlags=300`.
- Lines 20-21: Generate angle track; implemented by `angles=randn(3,npoints)`.
- Lines 23-24: Preallocate rotation matrix array; implemented by `R=zeros([3 3 npoints]); R(:,:,1)=eye(3)`.
- Lines 26-27: Loop over Monte-Carlo steps; implemented by `for n=2:npoints`.
- Lines 29-32: Generate a random rotation; implemented by `R_gen=[ 0 1 0; -1 0 0; 0 0 0]*angles(1,n)+ [ 0 0 1; 0 0 0; -1 0 0]*angles(2,n)+ [ 0 0 0; 0 0 1; 0 -1 0]*angles(3,n)`.
- Lines 37-39: Get Monte-Carlo correlation function; implemented by `[cf_mc,lags]=xcorr(squeeze(R(k,m,:)), squeeze(R(p,q,:)),nlags,'normalized')`.
- Lines 42-43: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=2:npoints`.

### Key state/data transformations

- Lines 15: computes `sigma_iso` using `sigma_iso=0.2; k=2; m=3; p=2; q=3`.
- Lines 18: computes `npoints` using `npoints=1e6; nlags=300`.
- Lines 21: computes `angles` using `angles=randn(3,npoints)`.
- Lines 24: computes `R` using `R=zeros([3 3 npoints]); R(:,:,1)=eye(3)`.
- Lines 30-32: computes `R_gen` using `R_gen=[ 0 1 0; -1 0 0; 0 0 0]*angles(1,n)+ [ 0 0 1; 0 0 0; -1 0 0]*angles(2,n)+ [ 0 0 0; 0 0 1; 0 -1 0]*angles(3,n)`.
- Lines 33: computes `R(:,:,n)` using `R(:,:,n)=R(:,:,n-1)*expm(sigma_iso*R_gen)`.
- Lines 38-39: computes `[cf_mc,lags]` using `[cf_mc,lags]=xcorr(squeeze(R(k,m,:)), squeeze(R(p,q,:)),nlags,'normalized')`.
- Lines 40: computes `cf_mc` using `cf_mc=(1/3)*ifftshift(cf_mc); lags=ifftshift(lags)`.

## Implementation structure

- Computes the following rotational correlation function
- G(k,m,p,q)=<R(k,m)*R(p,q)>
- where R is the 3D Cartesian rotation matrix, using the
- Monte-Carlo method.
- Calculation time: minutes.
- Set testing parameters
- Set number of points
- Generate angle track
- Preallocate rotation matrix array
- Loop over Monte-Carlo steps
- Generate a random rotation
- Get Monte-Carlo correlation function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `angles()`, `xcorr()`, `squeeze()`, `ifftshift()`, `kfigure()`, `lags()`, `cf_mc()`, `xlim()`, `kxlabel()`, `kylabel()`, `klegend()`.
