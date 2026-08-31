# examples/fundamentals/correlation_function_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/correlation_function_3.m`
- Signature: `correlation_function_3()`
- Total lines: 97

## Purpose

Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returned by Spinach kernel for the following correlation function: G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'> The sigma parameters refer to the rates of rotation and the four indices to the Wigner functions being correlated. Rhombic rotational diffusion tested here. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Set the testing parameters; implemented by `sigma_x=0.1; sigma_y=0.2; sigma_z=0.3`.
- Lines 21-22: Convert indices from [-L,L] to [1,2*L+1]; implemented by `k=L+1+k; m=L+1+m; p=L+1+p; q=L+1+q`.
- Lines 24-27: % Numerical Monte-Carlo calculation; implemented by `npoints=1e6; nlags=300`.
- Lines 26-27: Number of points and lags; implemented by `npoints=1e6; nlags=300`.
- Lines 29-30: Generate the angle track; implemented by `angles=randn(3,npoints)`.
- Lines 32-33: Preallocate and start DCM trajectory; implemented by `DCMT=zeros([3 3 npoints]); DCM=eye(3)`.
- Lines 35-36: Loop over Monte-Carlo steps; implemented by `for n=1:npoints`.
- Lines 38-39: Store the step; implemented by `DCMT(:,:,n)=DCM`.
- Lines 41-44: Generate and apply a random rotation; implemented by `R_gen=sigma_z*[ 0 1 0; -1 0 0; 0 0 0]*angles(1,n)+ sigma_y*[ 0 0 1; 0 0 0; -1 0 0]*angles(2,n)+ sigma_x*[ 0 0 0; 0 0 1; 0 -1 0]*angles(3,n)`.
- Lines 49-50: Convert DCMs into Wigner functions; implemented by `W=zeros([2*L+1, 2*L+1, npoints],'like',1i)`.
- Lines 53-54: Pull a DCM out; implemented by `DCM=squeeze(DCMT(:,:,n))`.
- Lines 56-57: Convert into Euler angles; implemented by `[alp,bet,gam]=dcm2euler(DCM)`.
- Lines 59-60: Convert into Wigner functions; implemented by `W(:,:,n)=wigner(L,alp,bet,gam)`.
- Lines 64-66: Get Monte-Carlo correlation function; implemented by `[cf_mc,lags]=xcorr(squeeze(W(k,m,:)), squeeze(W(p,q,:)),nlags,'normalized')`.
- Lines 69-72: % Analytical Spinach calculation; implemented by `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 71-72: Create a dummy spin system; implemented by `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 82-83: Get the analytical correlation function; implemented by `[weights,rates]=corrfun(spin_system,L,k,m,p,q)`.
- Lines 89-90: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:npoints`.
- Line 51: `parfor` loop over `n=1:size(DCMT,3)`.
- Line 85: `for` loop over `k=1:numel(weights{1})`.

### Key state/data transformations

- Lines 18: computes `sigma_x` using `sigma_x=0.1; sigma_y=0.2; sigma_z=0.3`.
- Lines 19: computes `L` using `L=2; k=-1; m=2; p=-1; q=2`.
- Lines 22: computes `k` using `k=L+1+k; m=L+1+m; p=L+1+p; q=L+1+q`.
- Lines 27: computes `npoints` using `npoints=1e6; nlags=300`.
- Lines 30: computes `angles` using `angles=randn(3,npoints)`.
- Lines 33: computes `DCMT` using `DCMT=zeros([3 3 npoints]); DCM=eye(3)`.
- Lines 39: computes `DCMT(:,:,n)` using `DCMT(:,:,n)=DCM`.
- Lines 42-44: computes `R_gen` using `R_gen=sigma_z*[ 0 1 0; -1 0 0; 0 0 0]*angles(1,n)+ sigma_y*[ 0 0 1; 0 0 0; -1 0 0]*angles(2,n)+ sigma_x*[ 0 0 0; 0 0 1; 0 -1 0]*angles(3,n)`.
- Lines 45: computes `DCM` using `DCM=DCM*expm(R_gen)`.
- Lines 50: computes `W` using `W=zeros([2*L+1, 2*L+1, npoints],'like',1i)`.
- Lines 57: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(DCM)`.
- Lines 60: computes `W(:,:,n)` using `W(:,:,n)=wigner(L,alp,bet,gam)`.
- Lines 65-66: computes `[cf_mc,lags]` using `[cf_mc,lags]=xcorr(squeeze(W(k,m,:)), squeeze(W(p,q,:)),nlags,'normalized')`.
- Lines 67: computes `cf_mc` using `cf_mc=(1/(2*L+1))*ifftshift(cf_mc); lags=ifftshift(lags)`.
- Lines 72: computes `sys.magnet` using `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 73: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 74: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 75: computes `inter.equilibrium` using `inter.equilibrium='zero'`.

## Implementation structure

- Computes rotational correlation functions using a Monte-Carlo method and
- compares them to the analytical results returned by Spinach kernel for
- the following correlation function:
- G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>
- The sigma parameters refer to the rates of rotation and the four indices
- to the Wigner functions being correlated. Rhombic rotational diffusion
- tested here.
- Calculation time: minutes.
- Set the testing parameters
- Convert indices from [-L,L] to [1,2*L+1]
- % Numerical Monte-Carlo calculation
- Number of points and lags

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `DCMT()`, `angles()`, `squeeze()`, `dcm2euler()`, `wigner()`, `xcorr()`, `ifftshift()`, `create()`, `basis()`, `corrfun()`, `kfigure()`, `lags()`, `cf_mc()`, `xlim()`, `kylabel()`, `kxlabel()`.
