# experiments/nmr_liquids/mqs_refocus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/mqs_refocus.m`
- Signature: `fid=mqs_refocus(spin_system,parameters,H,R,K)`
- Total lines: 205

## Purpose

Multiple quantum correlation pulse sequence with refocusing, as described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=mqs_refocus(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of fid points
- parameters.spins {F1 F2} nuclei, in this case:
- {'1H','1H'}
- parameters.angle flip angle for the final
- pulse, radians
- parameters.mqorder coherence orders to select,
- a two-element integer array
- parameters.delay_1 first evolution delay, seconds
- parameters.delay_2 second evolution delay, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian superoperator, provided by the context function
- R -relaxation superoperator, provided by the context function
- K -kinetics superoperator, provided by the context function

## Outputs

- fid -2D free induction decay for amplitude-mode processing
- Note: this implementation is homonuclear in practice and uses
- exact analytical coherence-order projection rather than an
- explicit phase cycle.

## Implementation structure

- Multiple quantum correlation pulse sequence with refocusing,
- as described in:
- fid=mqs_refocus(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of fid points
- parameters.spins {F1 F2} nuclei, in this case:
- {'1H','1H'}
- parameters.angle flip angle for the final
- pulse, radians
- parameters.mqorder coherence orders to select,
- a two-element integer array
- parameters.delay_1 first evolution delay, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
