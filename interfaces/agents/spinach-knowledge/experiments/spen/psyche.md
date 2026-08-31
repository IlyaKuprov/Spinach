# experiments/spen/psyche.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/psyche.m`
- Signature: `fid=psyche(spin_system,parameters,H,R,K,G,F)`
- Total lines: 266

## Purpose

PSYCHE pure-shift NMR pulse sequence. Syntax: fid=psyche_1d(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,F,G)`.
- Lines 56-57: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 59-60: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 64-67: Get chirp pulse waveform; implemented by `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.duration, parameters.bandwidth,parameters.smfactor, parameters.chirptype)`.
- Lines 69-70: Compute chirp RF amplitudes; implemented by `if strcmp(parameters.chirptype,'saltire')`.
- Lines 77-78: Calibrate chirps; implemented by `norm_factor=max(Cx)`.
- Lines 82-83: Apply the first pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 85-87: Run the first half of the t1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep1/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 89-90: Select "+1" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{'1H',+1}})`.
- Lines 92-94: Evolve for the delta period; implemented by `rho_stack=evolution(spin_system,L+parameters.g_amp*G{1},[], rho_stack,parameters.delta,1,'final')`.
- Lines 96-97: Apply the hard 180 pulse; implemented by `rho_stack=step(spin_system,Lx,rho_stack,pi)`.
- Lines 103-104: Select "-1" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{'1H',-1}})`.
- Lines 106-107: Apply the 1st pulse of the PSYCHE element; implemented by `durations=ones(size(Cx))*parameters.duration/numel(Cx)`.
- Lines 111-112: Select "0" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{'1H',0}})`.
- Lines 114-116: Apply the 2nd pulse of the PSYCHE element; implemented by `rho_stack=shaped_pulse_xy(spin_system,L+parameters.g_amp*G{1}, {Lx,Ly},{Cx,-Cy},durations,rho_stack,'expv-pwc')`.
- Lines 121-123: Run the second half of the t1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep1/2, parameters.npoints(1)-1,'refocus')`.
- Lines 125-127: Run the F2 evolution; implemented by `fid=evolution(spin_system,L,parameters.coil,rho_stack, parameters.timestep2,parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 70: conditional branch on `strcmp(parameters.chirptype,'saltire')`.

### Key state/data transformations

- Lines 57: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 60: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 61: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 62: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 65-67: computes `[Cx,Cy]` using `[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.duration, parameters.bandwidth,parameters.smfactor, parameters.chirptype)`.
- Lines 71: computes `rfbeta` using `rfbeta=(parameters.beta/360)*sqrt(2*parameters.bandwidth/(parameters.duration))`.
- Lines 73: computes `q_beta` using `q_beta=-(2*log(cosd(parameters.beta)/2+1/2))/pi`.
- Lines 78: computes `norm_factor` using `norm_factor=max(Cx)`.
- Lines 79: computes `Cx` using `Cx=Cx/norm_factor; Cy=Cy/norm_factor`.
- Lines 83: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 86-87: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep1/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 107: computes `durations` using `durations=ones(size(Cx))*parameters.duration/numel(Cx)`.
- Lines 126-127: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho_stack, parameters.timestep2,parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 132: `grumble()` — `function grumble(spin_system,parameters,H,R,K,F,G)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.rho0 initial state
- parameters.coil detection state
- parameters.spins nuclei on which the sequence runs
- parameters.g_amp gradient amplitude (T/m)
- parameters.dims size of the sample (m)
- parameters.npts number of discretization points in the grid
- parameters.sweep spectral range (Hz)
- parameters.npoints number of points in the sweep
- parameters.zerofill number of points for the zero filling
- parameters.diff diffusion constant (m^2/s)
- H Fokker-Planck Hamiltonian, received
- from the imaging context
- R Fokker-Planck relaxation superoperator,
- received from the imaging context
- K Fokker-Planck kinetics superoperator,
- received from the imaging context
- G Fokker-Planck gradient superoperators,
- received from the imaging context
- F Fokker-Planck diffusion and flow super-
- operator, received from the context

## Outputs

- fid -a PSYCHE free induction decay as a 2D array

## Implementation structure

- PSYCHE pure-shift NMR pulse sequence. Syntax:
- fid=psyche_1d(spin_system,parameters,H,R,K,G,F)
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.spins nuclei on which the sequence runs
- parameters.g_amp gradient amplitude (T/m)
- parameters.dims size of the sample (m)
- parameters.npts number of discretization points in the grid
- parameters.sweep spectral range (Hz)
- parameters.npoints number of points in the sweep
- parameters.zerofill number of points for the zero filling
- parameters.diff diffusion constant (m^2/s)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `strcmp()`, `cosd()`, `step()`, `evolution()`, `coherence()`, `shaped_pulse_xy()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`.
