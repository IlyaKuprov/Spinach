# experiments/nmr_solids/pdsd.m

- Signature: `fid=pdsd(spin_system,parameters,H,R,K)`

## Purpose

A simplified 2D PDSD experiment with NOESY-type quadrature detection and a four-step phase cycle, called from the `singlerot` context.

## Parameters / inputs

- `spin_system` — Spinach spin system object.
- `parameters.sweep` — Sweep width in Hz; sets the evolution timestep to `1/parameters.sweep`.
- `parameters.npoints` — Two-element vector giving the number of complex points in the indirect and direct dimensions.
- `parameters.tmix` — Mixing time in seconds.
- `parameters.rate` — MAS rate in Hz, used to set proton irradiation power during mixing.
- `parameters.spc_dim` — Spatial dimension of the MAS problem, received from the context function.
- `H`, `R`, `K` — Hamiltonian, relaxation, and kinetics superoperators received from the context function.

## Outputs

- `fid.cos`, `fid.sin` — Quadrature components of the 2D PDSD spectrum.

## Implementation summary

The sequence starts from a `13C` `Ly` state, omitting cross-polarisation. It evolves the indirect dimension under proton decoupling, applies the phase-cycled second pulse, evolves for `parameters.tmix` with proton irradiation, and applies a third pulse. After the proton subspace is removed, direct-dimension evolution and `13C` detection occur under proton decoupling. Differences between paired phase-cycle signals form `fid.cos` and `fid.sin` to eliminate axial peaks.

Source reference: <https://spindynamics.org/wiki/index.php?title=pdsd.m>