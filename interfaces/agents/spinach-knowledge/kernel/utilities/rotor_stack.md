# kernel/utilities/rotor_stack.m

- Signature: `[L,rotor_phases]=rotor_stack(spin_system,parameters,assumptions)`

## Purpose

Returns a rotor stack of Liouvillians or Hamiltonians for traditional-style calculations of magic-angle spinning (MAS) dynamics. Syntax: `L=rotor_stack(spin_system,parameters,assumptions)`.

## Parameters / inputs

- `parameters.axis`: Spinning axis as a normalized three-element row vector.
- `parameters.offset`: Transmitter offsets in Hz, corresponding element by element to `parameters.spins`. The header comment calls this a cell array; the executable validation instead rejects nonnumeric values and requires a nonempty numeric array of the same length as `parameters.spins`. The validation’s error text describes the array as real.
- `parameters.spins`: Nonempty cell array of spin identifiers that the offsets refer to, for example `{'1H','13C'}`. Each identifier must refer to a spin present in the system.
- `parameters.max_rank`: Maximum harmonic rank retained in the solution, approximately equal to the number of spinning sidebands in the spectrum. Increase it until the result converges. The implementation constructs `2*parameters.max_rank+1` rotor ticks.
- `parameters.rframes`: Cell array of numerical rotating-frame specifications. For example, `{{'13C',2},{'14N',3}}` requests a second-order transformation for carbon-13 and a third-order transformation for nitrogen-14. Use an empty cell array when no such transformations are needed.
- `parameters.orientation`: Initial orientation of the spin system at rotor phase zero, as three Euler angles in radians.
- `parameters.masframe`: Frame in which rotations are applied. `'magnet'` specifies an initial orientation in the laboratory frame and requires three-angle powder grids; `'rotor'` specifies an initial orientation in the rotor frame and requires two-angle powder grids.
- `assumptions`: Assumption set used to generate the Hamiltonian and validate numerical rotating frames, regardless of prior assumptions on the input object. It is passed to `assume` before Hamiltonian construction. Spins named in `parameters.rframes` must remain in the laboratory frame under this set; `rotframe` rejects already-rotating spins. See `assume.m`.

## Numerical-frame caveats

Numerical frames on the carrier-free `se_dnp_h+`, `se_dnp_h-`, and `se_dnp_h0` components are not implemented. Component stacks with empty `parameters.rframes` remain valid.

## Outputs

- `L`: Cell array of Hamiltonian or Liouvillian matrices, one per rotor tick.
- `rotor_phases`: Rotor phase at each tick, in radians.

Relaxation and chemical kinetics are not included.

Contact: `ilya.kuprov@weizmann.ac.il`

<https://spindynamics.org/wiki/index.php?title=rotor_stack.m>