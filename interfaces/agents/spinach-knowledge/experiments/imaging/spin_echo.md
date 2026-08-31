# experiments/imaging/spin_echo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/spin_echo.m`
- Signature: `fid=spin_echo(spin_system,parameters,H,R,K,G,F)`
- Total lines: 124

## Purpose

Spin echo pulse sequence. Syntax: fid=spin_echo(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.g_amp -the amplitude of gradient T/m parameters.g_step_dur -time step duration parameters.g_n_steps -number of time steps

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 28-29: Assemble the Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 31-32: Make pulse operators; implemented by `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 35-36: Hard 90-degree pulse; implemented by `rho=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 38-40: Evolution under the X gradient; implemented by `rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.g_step_dur, parameters.g_n_steps,'final')`.
- Lines 42-43: Hard 180-degree pulse; implemented by `rho=step(spin_system,Hy,rho,pi)`.
- Lines 45-48: Detection under the X gradient; implemented by `fid=evolution(spin_system,L+parameters.g_amp*G{1},parameters.coil,rho, parameters.g_step_dur, 2*parameters.g_n_steps,'observable')`.

### Key state/data transformations

- Lines 29: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 32: computes `Hp` using `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 33: computes `Hy` using `Hy=kron(speye(parameters.npts),(Hp-Hp')/2i)`.
- Lines 36: computes `rho` using `rho=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 46-48: computes `fid` using `fid=evolution(spin_system,L+parameters.g_amp*G{1},parameters.coil,rho, parameters.g_step_dur, 2*parameters.g_n_steps,'observable')`.

### Local helper functions

- Line 53: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Outputs

- fid -the time domain echo signal

## Implementation structure

- Spin echo pulse sequence. Syntax:
- fid=spin_echo(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.g_amp -the amplitude of gradient T/m
- parameters.g_step_dur -time step duration
- parameters.g_n_steps -number of time steps
- fid -the time domain echo signal
- Check consistency
- Assemble the Liouvillian
- Make pulse operators
- Hard 90-degree pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `isscalar()`.
