# experiments/esr_hyperfine/endor_mims_ideal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_mims_ideal.m`
- Signature: `endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)`
- Total lines: 196

## Purpose

Mims ENDOR sequence with ideal electron pulses. Syntax: endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 54-55: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 57-58: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 60-61: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 63-64: Ideal pulse operators on all electrons; implemented by `Ex=operator(spin_system,'Lx',parameters.electrons)`.
- Lines 67-68: Ideal initial and detection states; implemented by `rho0=state(spin_system,'Lz',parameters.electrons)`.
- Lines 71-72: Nutation frequencies at the specified RF B1 field; implemented by `nut_freqs=-spin_system.inter.gammas*parameters.rf_b1_field`.
- Lines 74-75: Nuclear pulse operators; implemented by `Nx=sparse(0); Ny=sparse(0)`.
- Lines 78-79: Gamma-weighted because isotopes may differ; implemented by `Nx=Nx+nut_freqs(n)*operator(spin_system,'Lx',n)`.
- Lines 84-85: Ideal pi/2 pulse on the electrons; implemented by `rho=step(spin_system,Ex,rho0,pi/2)`.
- Lines 87-88: Run stimulated echo delay; implemented by `rho=step(spin_system,L,rho,parameters.tau)`.
- Lines 90-91: Ideal pi/2 pulse on the electrons; implemented by `rho=step(spin_system,Ex,rho,pi/2)`.
- Lines 93-94: Kron up to make array over radiofrequencies; implemented by `rho=kron(rho,ones(1,numel(parameters.n_frq)))`.
- Lines 96-97: Loop over the radiofrequency array; implemented by `parfor n=1:numel(parameters.n_frq)`.
- Lines 99-103: Blast the nuclei and subtract the background; implemented by `rho(:,n)=shaped_pulse_af(spin_system,L,Nx,Ny,rho(:,n),parameters.n_frq(n),1, parameters.n_dur,0,parameters.n_rnk,'expm')- shaped_pulse_af(spin_system,L,Nx,Ny,rho(:,n),pa…`.
- Lines 107-108: Ideal pi/2 pulse on the electrons; implemented by `rho=step(spin_system,Ey,rho,-pi/2)`.
- Lines 113-114: Detect echo tips; implemented by `endor_spec=coil'*rho`.

### Control flow inferred from the code

- Line 76: `for` loop over `n=parameters.nuclei`.
- Line 97: `parfor` loop over `n=1:numel(parameters.n_frq)`.

### Key state/data transformations

- Lines 55: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 61: computes `L` using `L=H+1i*R+1i*K`.
- Lines 64: computes `Ex` using `Ex=operator(spin_system,'Lx',parameters.electrons)`.
- Lines 65: computes `Ey` using `Ey=operator(spin_system,'Ly',parameters.electrons)`.
- Lines 68: computes `rho0` using `rho0=state(spin_system,'Lz',parameters.electrons)`.
- Lines 69: computes `coil` using `coil=state(spin_system,'L+',parameters.electrons)`.
- Lines 72: computes `nut_freqs` using `nut_freqs=-spin_system.inter.gammas*parameters.rf_b1_field`.
- Lines 75: computes `Nx` using `Nx=sparse(0); Ny=sparse(0)`.
- Lines 80: computes `Ny` using `Ny=Ny+nut_freqs(n)*operator(spin_system,'Ly',n)`.
- Lines 85: computes `rho` using `rho=step(spin_system,Ex,rho0,pi/2)`.
- Lines 100-103: computes `rho(:,n)` using `rho(:,n)=shaped_pulse_af(spin_system,L,Nx,Ny,rho(:,n),parameters.n_frq(n),1, parameters.n_dur,0,parameters.n_rnk,'expm')- shaped_pulse_af(spin_system,L,Nx,Ny,rho(:,n),pa…`.

### Local helper functions

- Line 119: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
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
- seconds; 50e-6 is typical
- parameters.n_frq -nuclear pulse frequency offsets
- parameters.rf_b1_field -RF B1 field strength
- parameters.n_rnk -nuclear pulse grid rank
- parameters.nuclei -a vector of integers specifying
- which spins in sys.isotopes are
- nuclei to irradiate
- H -Hamiltonian matrix, received from the context
- function, normally powder() in this case
- R -relaxation superoperator, received from the context
- function, normally powder() in this case
- K -kinetics superoperator, received from the context
- function, normally powder() in this case

## Outputs

- endor_spec -Mims ENDOR spectrum, a vector of the same
- size as parameters.n_frq

## Implementation structure

- Mims ENDOR sequence with ideal electron pulses. Syntax:
- endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)
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

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `state()`, `nut_freqs()`, `step()`, `rho()`, `shaped_pulse_af()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `iscell()`, `ischar()`, `isscalar()`, `isrow()`.
