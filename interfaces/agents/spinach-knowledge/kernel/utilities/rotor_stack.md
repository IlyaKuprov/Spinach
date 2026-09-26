# kernel/utilities/rotor_stack.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rotor_stack.m`
- Signature: `[L,rotor_phases]=rotor_stack(spin_system,parameters,assumptions)`
- Total lines: 256

## Purpose

Returns a rotor stack of Liouvillians or Hamiltonians. The stack is needed for the traditional style calculation of MAS dynamics. Syntax: L=rotor_stack(spin_system,parameters,assumptions)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

The explicit `assumptions` argument governs Hamiltonian construction and numerical rotating-frame transformations alike, independently of any prior `assume` call on the input object. Nonempty `parameters.rframes` requires laboratory-frame assumptions on the transformed spins. Numerical frames on the carrier-free `se_dnp_h+`, `se_dnp_h-`, and `se_dnp_h0` components are not implemented; empty-frame component stacks remain valid.

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.axis -spinning axis, given as a normalized
- 3-element vector
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.spins -a cell array giving the spins that
- the offsets refer to, e.g. {'1H','13C'}
- parameters.max_rank -maximum harmonic rank to retain in
- the solution (increase till conver-
- gence is achieved, approximately
- equal to the number of spinning si-
- debands in the spectrum)
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.orientation -the orientation of the spin system
- at rotor phase zero, a vector of
- three Euler angles in radians.
- parameters.masframe -the frame in which the rotations
- are applied. The possibilities are:
- 'magnet' -the initial orientation in the lab frame
- (three-angle powder grids will be required)
- 'rotor' -the initial orientation in the rotor frame
- (two-angle powder grids will be required)
- assumptions -assumption set to be used in generating the
- Hamiltonian and numerical rotating-frame validation, regardless of the input object's prior assumptions. The transformed spins must remain in the laboratory frame under this set; already-rotating spins are rejected by `rotframe`. See `assume.m`.

## Outputs

- L -a cell array of Hamiltonian or Liouvillian matrices,
- one for each tick of the rotor.
- rotor_phases -rotor phases at each tick, radians
- Note: relaxation and chemical kinetics are not included.

## Header notes

The spinning axis is a normalised three-vector; transmitter offsets are in hertz for the entries in parameters.spins. Increase the retained harmonic rank until the requested rotor-stack result is converged.
