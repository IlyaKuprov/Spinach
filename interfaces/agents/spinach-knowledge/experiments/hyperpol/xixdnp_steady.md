# experiments/hyperpol/xixdnp_steady.m

- Signature: `dnp=xixdnp_steady(spin_system,parameters,H,R,K)`

## Purpose

TPPM DNP and its special case X-inverse-X (XiX) DNP experiment from (https://doi.org/10.1021/jacs.1c09900), steady state ver- sion. Call from powder context. Syntax: dnp=xixdnp_steady(spin_system,parameters,H,R,K)

## Physical / mathematical content

At each microwave offset, the function builds phase-dependent pulse and delay propagators from `L=H+1i*R+1i*K`, then computes the steady-state detected signal.

## Numerical / algorithmic content

The steady-state density is obtained with `steady(...,'newton')`; the function evaluates the coil overlap for each requested resonance offset.

## Parameters / inputs

- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.coil -detection state vector
- parameters.pulse_dur -pulse duration, seconds
- parameters.phase -phase of each second pulse in radians
- parameters.nloops -number of XiX/TPPM DNP blocks, the
- calculation is faster when this is
- an integer power of 2
- parameters.shot_spacing -delay between microwave irradiation
- periods, seconds
- parameters.addshift -shift of the centre of the field
- profile, Hz
- parameters.el_offs -microwave resonance offsets, a vector
- of frequencies in Hz
- Output:
- dnp -steady state observable on the de-
- tection state vector as a function
- of microwave resonance offset

## Implementation structure

The function validates inputs, constructs electron control operators, loops over `el_offs`, cleans the propagators, solves for the steady state, and stores `parameters.coil'*rho` in `dnp`.
