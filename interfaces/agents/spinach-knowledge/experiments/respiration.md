# experiments/respiration.m

- Signature: `fid=respiration(spin_system,parameters,H,R,K)`

## Purpose

RESPIRATION cross-polarisation method described in the paper from the Aarhus group (http://dx.doi.org/10.1021/jz3000905).

## Physical / mathematical content

The Liouvillian is `L=H+1i*R+1i*K`. The pulse operators are formed for `parameters.spins{1}` and `parameters.spins{2}`.

## Numerical / algorithmic content

Starting from `parameters.rho0`, each loop applies a +X pulse and a −X pulse on the first specified spin, with each pulse implemented by `step` for duration `1/(2*parameters.rate)`; it then applies an ideal `parameters.theta` pulse on both spins. After the loops, the first specified spin is decoupled and `evolution` acquires the FID using `parameters.coil`, dwell time `1/parameters.sweep`, and `parameters.npoints-1` intervals.

## Parameters / inputs

- `parameters.sweep` — sweep width, Hz
- `parameters.npoints` — number of points in the FID
- `parameters.rho0` — initial state
- `parameters.coil` — detection state
- `parameters.nloops` — number of RESPIRATION loops
- `parameters.theta` — angle of the ideal pulse, applied at the end of each loop
- `parameters.rate` — RESPIRATION pulse train rate, Hz
- `parameters.spins` — working spins, e.g. {'1H','13C'}
- `H` — Hamiltonian matrix, received from the context function
- `R` — relaxation superoperator, received from the context function
- `K` — kinetics superoperator, received from the context function

## Outputs

- `fid` — free induction decay detected using `parameters.coil`

## Implementation structure

The implementation checks the input consistency with `grumble`, constructs the pulse operators, applies the RESPIRATION loop, then decouples the first specified spin and acquires the FID with `evolution`.

## References

- [Aarhus group paper](http://dx.doi.org/10.1021/jz3000905)
- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=respiration.m)
