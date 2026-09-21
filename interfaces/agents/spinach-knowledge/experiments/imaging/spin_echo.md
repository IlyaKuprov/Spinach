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
