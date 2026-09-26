# experiments/cpmg.m

- Signature: `fid=cpmg(spin_system,parameters,H,R,K)`

## Purpose

CPMG echo train with detection. Syntax: fid=cpmg(spin_system,parameters,H,R,K)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_op -pulse operator
- parameters.nloops -number of CPMG loops
- parameters.timestep -time step
- parameters.npoints -number of steps per half-echo

## Outputs

- fid -free induction decay throughout the sequence

## Implementation structure

- CPMG echo train with detection. Syntax:
- fid=cpmg(spin_system,parameters,H,R,K)
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_op -pulse operator
- parameters.nloops -number of CPMG loops
- parameters.timestep -time step
- parameters.npoints -number of steps per half-echo
- fid -free induction decay throughout the sequence
- Check consistency
- Project the operator
- Compose Liouvillian
