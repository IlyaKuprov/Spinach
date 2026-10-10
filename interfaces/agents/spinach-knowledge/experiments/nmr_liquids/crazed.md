# experiments/nmr_liquids/crazed.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/crazed.m) · [Spinach Wiki: crazed.m](https://spindynamics.org/wiki/index.php?title=crazed.m)

- Signature: `fid=crazed(spin_system,parameters,H,R,K)`

## Purpose

CRAZED pulse sequence, implemented as an ideal analytical coherence-pathway version of the sequence described in [DOI 10.1126/science.8266096](https://doi.org/10.1126/science.8266096). The source represents gradient selection by explicit coherence projections; it does not model a spatial gradient waveform.

## Inputs

- `parameters.sweep`: positive scalar sweep width in Hz.
- `parameters.npoints`: two positive integers, giving the F1 and F2 point counts.
- `parameters.spins`: a one-element cell array naming the isotope used by the sequence, for example `{'1H'}` or `{'13C'}`.
- `parameters.angle`: finite real second-pulse angle in radians.
- `parameters.rho0`: numeric initial state with a row dimension matching the Liouville-space dimension.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the source combines them as `H + 1i*R + 1i*K`. The implementation requires the `sphten-liouv` formalism.

## Sequence and coherence selection

The source applies a 90-degree y pulse to the selected spin, evolves the F1 trajectory with timestep `1/parameters.sweep`, projects onto the +2 double-quantum coherence branch, applies the second y pulse with angle `parameters.angle`, then projects onto the +1 single-quantum branch. Detection uses the `L+` state of the same spin during F2 evolution, with the same timestep and the F2 point count. These projections are the source's analytical representation of gradient selection.

## Output and scope

- `fid`: two-dimensional free induction decay with F1 and F2 evolution as above.

This is a description of the parameterised sequence implementation, not a measured spectrum or a claim that a simulation was run and validated.
