# experiments/esr_hyperfine/endor_cw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_cw.m`
- Signature: `fid=endor_cw(spin_system,parameters,H,R,K)`
- Total lines: 97

## Purpose

Fast approximate simulation of isotropic continuous-wave ENDOR pulse sequence -essentially an NMR spectrum weighted by hyper- fine couplings is recorded. Syntax: fid=endor_cw(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 33-34: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 36-37: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 39-40: Nuclear pulse operator; implemented by `Sp=operator(spin_system,'L+','nuclei'); Sy=(Sp-Sp')/2i`.
- Lines 42-43: The initial state is a sum of nuclear Lz states weighted by HFCs; implemented by `rho=zeros([size(L,1) 1],'like',1i)`.
- Lines 53-54: Detect the nuclei; implemented by `coil=state(spin_system,'L+','nuclei','cheap')`.
- Lines 56-57: Evolution time step; implemented by `timestep=1/parameters.sweep`.
- Lines 59-60: Pulse the nuclei; implemented by `rho=step(spin_system,Sy,rho,pi/2)`.
- Lines 62-63: Acquire on the nuclei; implemented by `fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable')`.
- Lines 65-66: Frequency symmetrization; implemented by `fid=real(fid)`.

### Control flow inferred from the code

- Line 44: `for` loop over `n=find(cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes))`.
- Line 45: `for` loop over `k=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes))`.

### Key state/data transformations

- Lines 31: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 37: computes `L` using `L=H+1i*R+1i*K`.
- Lines 40: computes `Sp` using `Sp=operator(spin_system,'L+','nuclei'); Sy=(Sp-Sp')/2i`.
- Lines 43: computes `rho` using `rho=zeros([size(L,1) 1],'like',1i)`.
- Lines 46-47: computes `amplitude` using `amplitude=trace(spin_system.inter.coupling.matrix{n,k})/3+ trace(spin_system.inter.coupling.matrix{k,n})/3`.
- Lines 54: computes `coil` using `coil=state(spin_system,'L+','nuclei','cheap')`.
- Lines 57: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 63: computes `fid` using `fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 71: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of FID points to be computed
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay whose Fourier transform
- approximates a CW ENDOR spectrum

## Implementation structure

- Fast approximate simulation of isotropic continuous-wave ENDOR
- pulse sequence -essentially an NMR spectrum weighted by hyper-
- fine couplings is recorded. Syntax:
- fid=endor_cw(spin_system,parameters,H,R,K)
- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of FID points to be computed
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay whose Fourier transform
- approximates a CW ENDOR spectrum
- Move into adjoint representation if needed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `cellfun()`, `strncmp()`, `state()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`.
