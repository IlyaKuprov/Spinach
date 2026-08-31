# experiments/slowpass.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/slowpass.m`
- Signature: `spectrum=slowpass(spin_system,parameters,H,R,K)`
- Total lines: 199

## Purpose

Slow passage detection -calculates spectrum values at the user- specified frequency positions using the Fourier transform of the Liouville -von Neumann equation. The biggest advantage over the fid+fft style detection is easy parallelization and the possibi- lity of getting spectrum values at specific frequencies without recalculating the entire free induction decay. Syntax: spectrum=slowpass(spin_system,parameters,H,

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 47-50: Get the frequency grid; implemented by `freq_grid=2*pi*linspace(parameters.sweep(1), parameters.sweep(2), parameters.npoints)'`.
- Lines 52-53: Preallocate the answer; implemented by `spectrum=zeros(size(freq_grid),'like',1i)`.
- Lines 55-56: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 58-59: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 61-62: Compute subspace projectors; implemented by `projectors=reduce(spin_system,L,parameters.coil)`.
- Lines 64-65: Loop over subspaces; implemented by `for k=1:numel(projectors)`.
- Lines 67-68: Project into the current subspace; implemented by `rho0_subs=projectors{k}'*parameters.rho0`.
- Lines 73-74: Run backslash on the GPU if instructed; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 76-77: Inform the user; implemented by `report(spin_system,'using GPU backslash path ')`.
- Lines 79-80: Move the objects to GPU; implemented by `rho0_subs=gpuArray(full(rho0_subs)); L_subs=gpuArray(full(L_subs))`.
- Lines 83-84: Run the calculation using backslash; implemented by `parfor n=1:numel(freq_grid)`.
- Lines 91-92: For large problems use GMRES; implemented by `if (size(rho0_subs,1)>5000)&&(size(rho0_subs,2)==1)`.
- Lines 94-95: Inform the user; implemented by `report(spin_system,'using preconditioned CPU GMRES path ')`.
- Lines 97-99: Get preconditioners; implemented by `[M1,M2]=ilu(1i*L_subs+1i*mean(freq_grid)*Id_subs, struct('type','crout','droptol',1e-3))`.
- Lines 104-105: MDCS diagnostics; implemented by `parallel_profiler_start`.
- Lines 107-108: Run the calculation in parallel; implemented by `parfor n=1:numel(freq_grid)`.
- Lines 110-112: Run using preconditioned GMRES; implemented by `spectrum(n)=spectrum(n)+coil_subs'*gmres(1i*L_subs+1i*freq_grid(n)*Id_subs, rho0_subs,10,1e-6,numel(rho0_subs),M1,M2)`.

### Control flow inferred from the code

- Line 65: `for` loop over `k=1:numel(projectors)`.
- Line 74: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 84: `parfor` loop over `n=1:numel(freq_grid)`.
- Line 92: conditional branch on `(size(rho0_subs,1)>5000)&&(size(rho0_subs,2)==1)`.
- Line 108: `parfor` loop over `n=1:numel(freq_grid)`.
- Line 128: `parfor` loop over `n=1:numel(freq_grid)`.

### Key state/data transformations

- Lines 48-50: computes `freq_grid` using `freq_grid=2*pi*linspace(parameters.sweep(1), parameters.sweep(2), parameters.npoints)'`.
- Lines 53: computes `spectrum` using `spectrum=zeros(size(freq_grid),'like',1i)`.
- Lines 56: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 59: computes `L` using `L=H+1i*R+1i*K`.
- Lines 62: computes `projectors` using `projectors=reduce(spin_system,L,parameters.coil)`.
- Lines 68: computes `rho0_subs` using `rho0_subs=projectors{k}'*parameters.rho0`.
- Lines 69: computes `coil_subs` using `coil_subs=projectors{k}'*parameters.coil`.
- Lines 70: computes `L_subs` using `L_subs=projectors{k}'*L*projectors{k}`.
- Lines 71: computes `Id_subs` using `Id_subs=speye(size(L_subs))`.
- Lines 85: computes `spectrum_subs` using `spectrum_subs=dot(coil_subs,((1i*L_subs+1i*freq_grid(n)*Id_subs)\rho0_subs))`.
- Lines 86: computes `spectrum(n)` using `spectrum(n)=spectrum(n)+gather(spectrum_subs)`.
- Lines 98-99: computes `[M1,M2]` using `[M1,M2]=ilu(1i*L_subs+1i*mean(freq_grid)*Id_subs, struct('type','crout','droptol',1e-3))`.
- Lines 100-102: computes `report(spin_system,['nnz(L)` using `report(spin_system,['nnz(L)=' num2str(nnz(L_subs)) ', nnz(M1)=' num2str(nnz(M1)) ', nnz(M2)=' num2str(nnz(M2))])`.
- Lines 101-102: computes `', nnz(M1)` using `', nnz(M1)=' num2str(nnz(M1)) ', nnz(M2)=' num2str(nnz(M2))])`.
- Lines 102: computes `', nnz(M2)` using `', nnz(M2)=' num2str(nnz(M2))])`.
- Lines 145-146: computes `sample_rate` using `sample_rate=abs(parameters.sweep(2)-parameters.sweep(1))* parameters.npoints/(parameters.npoints-1)`.

### Local helper functions

- Line 154: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the spectrum of the system with the specified
- starting state detected on the specified coil
- state within the frequency interval requested
- Note: relaxation must be present in the system dynamics, or the
- matrix inversion operation would fail to converge. The re-
- laxation matrix R must *not* be thermalized.

## Implementation structure

- Slow passage detection -calculates spectrum values at the user-
- specified frequency positions using the Fourier transform of the
- Liouville -von Neumann equation. The biggest advantage over the
- fid+fft style detection is easy parallelization and the possibi-
- lity of getting spectrum values at specific frequencies without
- recalculating the entire free induction decay. Syntax:
- spectrum=slowpass(spin_system,parameters,H,R,K)
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sim2liouv()`, `reduce()`, `speye()`, `ismember()`, `report()`, `gpuArray()`, `dot()`, `freq_grid()`, `spectrum()`, `gather()`, `ilu()`, `nnz()`, `num2str()`, `gmres()`, `ismatrix()`.
