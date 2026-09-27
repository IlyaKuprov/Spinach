# experiments/hyperpol/xixdnp.m

- Signature: `contact_curve=xixdnp(spin_system,parameters,H,R,K)`

## Purpose

TPPM DNP and its special case X-inverse-X (XiX) DNP experiment from (https://doi.org/10.1021/jacs.1c09900). Syntax (call from powder context): contact_curve=xixdnp(spin_system,parameters,H,R,K)

## Physical / mathematical content

Each XiX block applies two phase-dependent microwave pulses to the state, then records the overlap with the coil detection state. The sequence is not a steady-state or MAS calculation.

## Numerical / algorithmic content

The function forms `L=H+1i*R+1i*K` and records `hdot` after each XiX block to build the contact curve.

## Parameters / inputs

- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_dur -pulse duration, seconds
- parameters.nloops -number of XiX DNP blocks
- parameters.phase -phase of the second pulse
- Output:
- contact_curve -time dependence of the coil state

## Implementation structure

The function validates inputs, constructs electron control operators, initializes the curve from `hdot(parameters.coil,parameters.rho0)`, then applies each XiX block and records its coil-state overlap.
