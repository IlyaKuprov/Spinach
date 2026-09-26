# experiments/nmr_solids/mqmas.m

- Signature: `fid=mqmas(spin_system,parameters,H,R,K)`

## Purpose

Rotor-synchronous MQMAS pulse sequence. Syntax: fid=mqmas(spin_system,parameters,H,R,K) This function should normally be called using singlerot.m context that would provide H, R, and K.

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.pulse_dur -duration of each pulse, a two-
- element vector, seconds
- parameters.pulse_amp -amplitude of each pulse, a two-
- element vector, rad/s
- parameters.mq_order -MQMAS coherence order
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- + the parameters required by the singlerot.m
- context function that will provide H, R, and K

## Outputs

- fid -2D amplitude mode free induction decay

## Header notes

- parameters.sweep should be a positive real scalar
- equal to abs(parameters.rate), this is because
- this pulse sequence is stroboscopic with respect
- to the rotor period; both dimensions are sampled
- at that sweep width

## Implementation structure

- Rotor-synchronous MQMAS pulse sequence. Syntax:
- fid=mqmas(spin_system,parameters,H,R,K)
- This function should normally be called using singlerot.m context
- that would provide H, R, and K.
- parameters.pulse_dur -duration of each pulse, a two-
- element vector, seconds
- parameters.pulse_amp -amplitude of each pulse, a two-
- element vector, rad/s
- parameters.mq_order -MQMAS coherence order
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- + the parameters required by the singlerot.m
