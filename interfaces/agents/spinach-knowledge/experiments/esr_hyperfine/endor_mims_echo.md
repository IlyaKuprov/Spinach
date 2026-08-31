# experiments/esr_hyperfine/endor_mims_echo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_mims_echo.m`
- Signature: `stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)`
- Total lines: 145

## Purpose

Stimulated echo diagnostics for the Mims ENDOR sequence. Syntax: stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 50-51: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 53-54: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 56-57: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 59-60: Ideal pulse operators on all electrons; implemented by `Ex=operator(spin_system,'Lx',parameters.electrons)`.
- Lines 63-64: Ideal initial and detection states; implemented by `rho0=state(spin_system,'Lz',parameters.electrons)`.
- Lines 67-68: Ideal pi/2 pulse on the electrons; implemented by `rho=step(spin_system,Ex,rho0,pi/2)`.
- Lines 70-71: Run stimulated echo delay; implemented by `rho=step(spin_system,L,rho,parameters.tau)`.
- Lines 73-74: Ideal pi/2 pulse on the electrons; implemented by `rho=step(spin_system,Ex,rho,pi/2)`.
- Lines 76-77: Delay corresponding to the missing nuclear pulse; implemented by `rho=evolution(spin_system,L,[],rho,parameters.n_dur,1,'final')`.
- Lines 79-80: Ideal pi/2 pulse on the electrons; implemented by `rho=step(spin_system,Ey,rho,-pi/2)`.
- Lines 82-85: Digitise the stimulated echo; implemented by `stim_echo=evolution(spin_system,L,coil,rho, 2*parameters.tau/parameters.nsteps, parameters.nsteps,'observable')`.

### Key state/data transformations

- Lines 51: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 57: computes `L` using `L=H+1i*R+1i*K`.
- Lines 60: computes `Ex` using `Ex=operator(spin_system,'Lx',parameters.electrons)`.
- Lines 61: computes `Ey` using `Ey=operator(spin_system,'Ly',parameters.electrons)`.
- Lines 64: computes `rho0` using `rho0=state(spin_system,'Lz',parameters.electrons)`.
- Lines 65: computes `coil` using `coil=state(spin_system,'L+',parameters.electrons)`.
- Lines 68: computes `rho` using `rho=step(spin_system,Ex,rho0,pi/2)`.
- Lines 83-85: computes `stim_echo` using `stim_echo=evolution(spin_system,L,coil,rho, 2*parameters.tau/parameters.nsteps, parameters.nsteps,'observable')`.

### Local helper functions

- Line 90: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.spins -working spins, normally {'E'}; spe-
- cify multiplicity if electron spin
- is not 1/2, for example {'7E'} for
- gadolinium
- parameters.electrons -a vector of integers specifying
- which spins in sys.isotopes are
- electrons
- parameters.tau -the delay between the first two
- 90-degree pulses of the Mims
- ENDOR sequence, seconds; 200e-9
- is typical
- parameters.n_dur -duration of the nuclear pulse,
- which is NOT ACTUALLY APPLIED
- within this pulse sequence,
- seconds; 50e-6 is typical
- parameters.nsteps -number of time steps to make
- in the detection period, which
- runs from 0 to 2*paramters.tau
- H -Hamiltonian matrix, received from the context
- function, normally powder() in this case
- R -relaxation superoperator, received from the context
- function, normally powder() in this case
- K -kinetics superoperator, received from the context
- function, normally powder() in this case

## Outputs

- stim_echo -stimulated echo seen in the Mims ENDOR sequ-
- ence in the absence of the nuclear RF pulse

## Implementation structure

- Stimulated echo diagnostics for the Mims ENDOR sequence. Syntax:
- stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)
- parameters.spins -working spins, normally {'E'}; spe-
- cify multiplicity if electron spin
- is not 1/2, for example {'7E'} for
- gadolinium
- parameters.electrons -a vector of integers specifying
- which spins in sys.isotopes are
- electrons
- parameters.tau -the delay between the first two
- 90-degree pulses of the Mims
- ENDOR sequence, seconds; 200e-9

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `state()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `iscell()`, `ischar()`, `isscalar()`, `isrow()`, `any()`.
