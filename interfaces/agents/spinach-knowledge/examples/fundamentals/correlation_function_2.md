# examples/fundamentals/correlation_function_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/correlation_function_2.m`
- Signature: `correlation_function_2()`
- Total lines: 96

## Purpose

Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returned by Spinach kernel for the following correlation function: G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'> The sigma parameters refer to the rates of rotation and the four indices to the Wigner functions being correlated. Axial rotational diffusion tes- ted here. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Set testing parameters; implemented by `sigma_ax=0.1; sigma_eq=0.2; L=2; k=-1; m=2; p=-1; q=2`.
- Lines 20-21: Convert indices from [-L,L] to [1,2*L+1]; implemented by `k=L+1+k; m=L+1+m; p=L+1+p; q=L+1+q`.
- Lines 23-26: % Numerical Monte-Carlo calculation; implemented by `npoints=1e6; nlags=300`.
- Lines 25-26: Number of points and lags; implemented by `npoints=1e6; nlags=300`.
- Lines 28-29: Generate the angle track; implemented by `angles=randn(3,npoints)`.
- Lines 31-32: Preallocate and start DCM trajectory; implemented by `DCMT=zeros([3 3 npoints]); DCM=eye(3)`.
- Lines 34-35: Loop over Monte-Carlo steps; implemented by `for n=1:npoints`.
- Lines 37-38: Store the step; implemented by `DCMT(:,:,n)=DCM`.
- Lines 40-43: Generate and apply a random rotation; implemented by `R_gen=sigma_ax*[ 0 1 0; -1 0 0; 0 0 0]*angles(1,n)+ sigma_eq*[ 0 0 1; 0 0 0; -1 0 0]*angles(2,n)+ sigma_eq*[ 0 0 0; 0 0 1; 0 -1 0]*angles(3,n)`.
- Lines 48-49: Convert DCMs into Wigner functions; implemented by `W=zeros([2*L+1, 2*L+1, npoints],'like',1i)`.
- Lines 52-53: Pull a DCM out; implemented by `DCM=squeeze(DCMT(:,:,n))`.
- Lines 55-56: Convert into Euler angles; implemented by `[alp,bet,gam]=dcm2euler(DCM)`.
- Lines 58-59: Convert into Wigner functions; implemented by `W(:,:,n)=wigner(L,alp,bet,gam)`.
- Lines 63-65: Get Monte-Carlo correlation function; implemented by `[cf_mc,lags]=xcorr(squeeze(W(k,m,:)), squeeze(W(p,q,:)),nlags,'normalized')`.
- Lines 68-71: % Analytical Spinach calculation; implemented by `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 70-71: Create a dummy spin system; implemented by `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 81-82: Get the analytical correlation function; implemented by `[weights,rates]=corrfun(spin_system,L,k,m,p,q)`.
- Lines 88-89: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 35: `for` loop over `n=1:npoints`.
- Line 50: `parfor` loop over `n=1:size(DCMT,3)`.
- Line 84: `for` loop over `k=1:numel(weights{1})`.

### Key state/data transformations

- Lines 18: computes `sigma_ax` using `sigma_ax=0.1; sigma_eq=0.2; L=2; k=-1; m=2; p=-1; q=2`.
- Lines 21: computes `k` using `k=L+1+k; m=L+1+m; p=L+1+p; q=L+1+q`.
- Lines 26: computes `npoints` using `npoints=1e6; nlags=300`.
- Lines 29: computes `angles` using `angles=randn(3,npoints)`.
- Lines 32: computes `DCMT` using `DCMT=zeros([3 3 npoints]); DCM=eye(3)`.
- Lines 38: computes `DCMT(:,:,n)` using `DCMT(:,:,n)=DCM`.
- Lines 41-43: computes `R_gen` using `R_gen=sigma_ax*[ 0 1 0; -1 0 0; 0 0 0]*angles(1,n)+ sigma_eq*[ 0 0 1; 0 0 0; -1 0 0]*angles(2,n)+ sigma_eq*[ 0 0 0; 0 0 1; 0 -1 0]*angles(3,n)`.
- Lines 44: computes `DCM` using `DCM=DCM*expm(R_gen)`.
- Lines 49: computes `W` using `W=zeros([2*L+1, 2*L+1, npoints],'like',1i)`.
- Lines 56: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(DCM)`.
- Lines 59: computes `W(:,:,n)` using `W(:,:,n)=wigner(L,alp,bet,gam)`.
- Lines 64-65: computes `[cf_mc,lags]` using `[cf_mc,lags]=xcorr(squeeze(W(k,m,:)), squeeze(W(p,q,:)),nlags,'normalized')`.
- Lines 66: computes `cf_mc` using `cf_mc=(1/(2*L+1))*ifftshift(cf_mc); lags=ifftshift(lags)`.
- Lines 71: computes `sys.magnet` using `sys.magnet=0; sys.isotopes={'G'}`.
- Lines 72: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 73: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 74: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 75: computes `inter.tau_c` using `inter.tau_c={1./(3*[sigma_ax sigma_eq].^2)}`.

## Implementation structure

- Computes rotational correlation functions using a Monte-Carlo method and
- compares them to the analytical results returned by Spinach kernel for
- the following correlation function:
- G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>
- The sigma parameters refer to the rates of rotation and the four indices
- to the Wigner functions being correlated. Axial rotational diffusion tes-
- ted here.
- Calculation time: minutes.
- Set testing parameters
- Convert indices from [-L,L] to [1,2*L+1]
- % Numerical Monte-Carlo calculation
- Number of points and lags

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `DCMT()`, `angles()`, `squeeze()`, `dcm2euler()`, `wigner()`, `xcorr()`, `ifftshift()`, `create()`, `basis()`, `corrfun()`, `kfigure()`, `lags()`, `cf_mc()`, `xlim()`, `kylabel()`, `kxlabel()`.
