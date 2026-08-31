# experiments/spen/psycosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/psycosy.m`
- Signature: `fid=psycosy(spin_system,parameters,H,R,K,G,F)`
- Total lines: 253

## Purpose

Alan Kenwright's spatially encoded COSY sequence described in fid=psycosy(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 57-58: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,F,G)`.
- Lines 60-61: Coherent evolution timestep; implemented by `parameters.delta=1/(4*parameters.sweep)`.
- Lines 63-64: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 66-67: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 71-73: Get chirp pulse qaveform; implemented by `[Cx,Cy]=chirp_pulse(parameters.sal_npt,parameters.sal_dur, parameters.sal_swp,parameters.sal_smf,'saltire')`.
- Lines 75-76: Normalize chirp amplitude; implemented by `norm_factor=max(Cx); Cx=Cx/norm_factor; Cy=Cy/norm_factor`.
- Lines 78-79: Q factor for saltire flip angle; implemented by `q_beta=-(2*log(cosd(parameters.sal_ang)/2+1/2))/pi`.
- Lines 81-82: RF field strength of saltire chirp based on Q factor (Hz); implemented by `rfbeta=sqrt(parameters.sal_swp*q_beta/(2*pi*parameters.sal_dur))`.
- Lines 84-85: Calibrate chirps; implemented by `Cx=2*pi*rfbeta*Cx; Cy=2*pi*rfbeta*Cy`.
- Lines 87-88: Apply the first pulse; implemented by `rho_stack=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 90-92: Run the first half of the t1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,0.5/parameters.sweep, parameters.npoints(1)-1,'trajectory')`.
- Lines 94-95: Select "+1" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}})`.
- Lines 97-98: Apply the hard 180 pulse; implemented by `rho_stack=step(spin_system,Lx,rho_stack,pi)`.
- Lines 100-101: Select "-1" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},-1}})`.
- Lines 103-104: Apply the 1st pulse of the PSYCHE element; implemented by `durations=ones(size(Cx))*parameters.sal_dur/numel(Cx)`.
- Lines 108-110: Run the diffusion time evolution; implemented by `rho_stack=evolution(spin_system,L+parameters.gamp*G{1},[],rho_stack, parameters.sal_del-2*parameters.sal_dur,1,'final')`.
- Lines 112-113: Select "0" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},0}})`.
- Lines 115-117: Apply the 2nd pulse of the PSYCHE element; implemented by `rho_stack=shaped_pulse_xy(spin_system,L+parameters.gamp*G{1}, {Lx,Ly},{Cx,-Cy},durations,rho_stack,'expv-pwc')`.

### Key state/data transformations

- Lines 61: computes `parameters.delta` using `parameters.delta=1/(4*parameters.sweep)`.
- Lines 64: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 67: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 68: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 69: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 72-73: computes `[Cx,Cy]` using `[Cx,Cy]=chirp_pulse(parameters.sal_npt,parameters.sal_dur, parameters.sal_swp,parameters.sal_smf,'saltire')`.
- Lines 76: computes `norm_factor` using `norm_factor=max(Cx); Cx=Cx/norm_factor; Cy=Cy/norm_factor`.
- Lines 79: computes `q_beta` using `q_beta=-(2*log(cosd(parameters.sal_ang)/2+1/2))/pi`.
- Lines 82: computes `rfbeta` using `rfbeta=sqrt(parameters.sal_swp*q_beta/(2*pi*parameters.sal_dur))`.
- Lines 85: computes `Cx` using `Cx=2*pi*rfbeta*Cx; Cy=2*pi*rfbeta*Cy`.
- Lines 88: computes `rho_stack` using `rho_stack=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 104: computes `durations` using `durations=ones(size(Cx))*parameters.sal_dur/numel(Cx)`.
- Lines 136-137: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep, parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 142: `grumble()` — `function grumble(spin_system,parameters,H,R,K,F,G)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix mixing time, seconds
- parameters.gamp gradient amplitude, T/m
- parameters.sal_ang flip angle of the saltire chirp (degrees)
- parameters.sal_dur pulse width of saltire chirp (s)
- parameters.sal_del chirp pulse gradient duration (s)
- parameters.sal_swp sweep width of saltire chirp (Hz)
- parameters.sal_npt number of points in the saltire chirp
- parameters.sal_smf saltire chirp smoothing factor
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

- fid -two-dimensional free induction decay

## Implementation structure

- Alan Kenwright's spatially encoded COSY sequence described in
- fid=psycosy(spin_system,parameters,H,R,K,G,F)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix mixing time, seconds
- parameters.gamp gradient amplitude, T/m
- parameters.sal_ang flip angle of the saltire chirp (degrees)
- parameters.sal_dur pulse width of saltire chirp (s)
- parameters.sal_del chirp pulse gradient duration (s)
- parameters.sal_swp sweep width of saltire chirp (Hz)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `cosd()`, `step()`, `evolution()`, `coherence()`, `shaped_pulse_xy()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`, `ischar()`.
