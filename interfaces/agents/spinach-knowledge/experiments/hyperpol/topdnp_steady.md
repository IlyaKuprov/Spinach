# experiments/hyperpol/topdnp_steady.m

- Signature: `dnp=topdnp_steady(spin_system,parameters,H,R,K)`

## Purpose

A steady-state implementation of the time-optimised pulsed DNP experiment (https://doi.org/10.1126/sciadv.aav6909), called from a powder context.

## Physical / mathematical content

For each microwave offset, the function adds rotating-frame and +X irradiation terms to `L=H+1i*R+1i*K`, then forms the TOP DNP pulse-and-delay propagator. It computes the detected coil overlap at steady state.

## Numerical / algorithmic content

The function repeats the cleaned TOP DNP block propagator `nloops` times, applies the shot-spacing delay, and obtains the steady state with `steady(...,'newton')`.

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
- parameters.delay_dur -delay duration, seconds
- parameters.nloops -number of XiX/TPPM DNP blocks,
- must be an integer power of 2
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

The function validates inputs, constructs electron operators, loops over resonance offsets, and computes each offset's steady-state signal; propagators may be exponentiated on the GPU when enabled.
