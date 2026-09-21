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
