# experiments/nmr_liquids/mqs.m

- Signature: `fid=mqs(spin_system,parameters,H,R,K)`

## Purpose

2D multiple-quantum NMR pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=mqs(spin_system,parameters,H,R,K)
This function should be invoked through liquid.m context,
which would provide H, R, and K.
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of points
- parameters.spins working spins, e.g. {'1H','1H'}
- parameters.angle flip angle, radians
- parameters.delay J-coupling evolution delay, seconds
- parameters.mqorder coherence order to select
- parameters.rho0 initial state
- parameters.coil detection state

## Outputs

- fid -2D magnitude-mode free induction decay
- Note: this is the non-refocused multiple-quantum/MaxQ variant.
- Use mqs_refocus.m when post-mixing refocusing is required.

## Implementation structure

- 2D multiple-quantum NMR pulse sequence from:
- fid=mqs(spin_system,parameters,H,R,K)
- This function should be invoked through liquid.m context,
- which would provide H, R, and K.
- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of points
- parameters.spins working spins, e.g. {'1H','1H'}
- parameters.angle flip angle, radians
- parameters.delay J-coupling evolution delay, seconds
- parameters.mqorder coherence order to select
- parameters.rho0 initial state
- parameters.coil detection state
