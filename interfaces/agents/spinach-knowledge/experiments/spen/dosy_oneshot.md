# experiments/spen/dosy_oneshot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/dosy_oneshot.m`
- Signature: `fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)`
- Total lines: 226

## Purpose

One-shot DOSY pulse sequence. Syntax: fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 55-56: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,F,G)`.
- Lines 58-59: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 61-62: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 65-66: Apply the first 90 degree pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 68-69: Select coherence -1; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},-1}})`.
- Lines 71-73: Evolve under a positive gradient; implemented by `rho=evolution(spin_system,L+(1+parameters.kappa)*parameters.g_amp*G{1}, [],rho,parameters.g_dur/2,1,'final')`.
- Lines 75-76: Evolve under the first gradient stabilization delay; implemented by `rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final')`.
- Lines 78-79: Apply the second 180 degree pulse; implemented by `rho=step(spin_system,Ly,rho,pi)`.
- Lines 81-82: Select coherence +1; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 84-86: Evolve under a negative gradient; implemented by `rho=evolution(spin_system,L-(1-parameters.kappa)*parameters.g_amp*G{1},[], rho,parameters.g_dur/2,1,'final')`.
- Lines 88-89: Evolve under the second gradient stabilization delay; implemented by `rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final')`.
- Lines 91-92: Apply the third 90 degree pulse; implemented by `rho=step(spin_system,Ly,rho,pi/2)`.
- Lines 94-95: Select coherence 0; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},0}})`.
- Lines 97-99: Evolve under a positive gradient; implemented by `rho=evolution(spin_system,L-(2*parameters.kappa)*parameters.g_amp*G{1},[], rho,parameters.g_dur/2,1,'final')`.
- Lines 104-106: Run the diffusion time evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.del-4*(parameters.g_dur/2) -4*parameters.g_stab_del,1,'final')`.
- Lines 108-110: Evolve under a negative gradient; implemented by `rho=evolution(spin_system,L-(2*parameters.kappa)*parameters.g_amp*G{1},[], rho,parameters.g_dur/2,1,'final')`.
- Lines 115-116: Apply the second 90 degree pulse; implemented by `rho=step(spin_system,Ly,rho,pi/2)`.
- Lines 118-120: Evolve under a positive gradient; implemented by `rho=evolution(spin_system,L+(1+parameters.kappa)*parameters.g_amp*G{1},[], rho,parameters.g_dur/2,1,'final')`.

### Key state/data transformations

- Lines 59: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 62: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 63: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 66: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 136-137: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho,1/parameters.sweep, parameters.npoints-1,'observable')`.

### Local helper functions

- Line 142: `grumble()` — `function grumble(spin_system,parameters,H,R,K,F,G)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- parameters.rho0 initial state
- parameters.coil detection state
- parameters.spins nuclei on which the sequence runs
- parameters.g_amp gradient amplitude for diffusion
- encoding (T/m)
- parameters.g_dur pulse width of the gradient for diffusion
- encoding (s)
- parameters.kappa unbalancing factor to unbalance the bipolar
- gradients in the ratio (1+kappa):(1-kappa)
- parameters.g_stab_del gradient stabilization delay (s)
- parameters.del diffusion delay, seconds
- parameters.dims size of the sample (m)
- parameters.npts number of discretization points in the grid
- parameters.npoints number of points in the acquired signal
- parameters.sweep acquisition sweep width, Hz
- H Fokker-Planck Hamiltonian
- R Fokker-Planck relaxation superoperator
- K Fokker-Planck kinetics superoperator
- G Fokker-Planck gradient superoperators
- F Fokker-Planck diffusion and flow
- superoperator

## Outputs

- fid -free induction decay

## Implementation structure

- One-shot DOSY pulse sequence. Syntax:
- fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.spins nuclei on which the sequence runs
- parameters.g_amp gradient amplitude for diffusion
- encoding (T/m)
- parameters.g_dur pulse width of the gradient for diffusion
- encoding (s)
- parameters.kappa unbalancing factor to unbalance the bipolar
- gradients in the ratio (1+kappa):(1-kappa)
- parameters.g_stab_del gradient stabilization delay (s)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `coherence()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`.
