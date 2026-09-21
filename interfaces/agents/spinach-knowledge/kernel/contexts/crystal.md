# kernel/contexts/crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/crystal.m`
- Signature: `answer=crystal(spin_system,pulse_sequence,parameters,assumptions)`
- Total lines: 259

## Purpose

Single-crystal interface to pulse sequences. Generates a Liouvillian superoperator and passes it on to the pulse sequence function, which should be supplied as a handle. Syntax: answer=crystal(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pulse_sequence -a function handle to one of the pulse se-
- quences located in the experiments folder
- assumptions -a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters.spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.orientation -a row vector of the three Euler angles
- (in radians) giving the orientation of
- the system relative to the input orien-
- tation.
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.needs -a cell array of strings specifying additional
- information required by the sequence:
- 'zeeman_op' -Zeeman part of the Hamiltonian
- in the laboratory frame, to be placed into
- parameters.hzeeman and sent to pulse sequence
- 'aniso_eq' -thermal equilibrium is recomputed
- using the full anisotropic Hamiltonian at the
- current orientation, and sent to the pulse
- sequence in parameters.rho0 subfield
- parameters.* -additional subfields may be required by your
- pulse sequence -check its documentation page
- The parameters structure is passed to the pulse sequence with the follo-
- wing additional parameters set:
- parameters.spc_dim -matrix dimension for the spatial
- dynamics subspace (1 in this case)
- parameters.spn_dim -matrix dimension for the spin
- dynamics subspace

## Outputs

- this function returns whatever it is that the pulse sequence returns
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.

## Implementation structure

- Single-crystal interface to pulse sequences. Generates a Liouvillian
- superoperator and passes it on to the pulse sequence function, which
- should be supplied as a handle. Syntax:
- answer=crystal(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -a function handle to one of the pulse se-
- quences located in the experiments folder
- assumptions -a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters.spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- parameters.offset -a cell array giving transmitter off-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `ismember()`, `equilibrium()`, `frqoffset()`, `orientation()`, `clear()`, `carrier()`, `rotframe()`, `relaxation()`, `kinetics()`, `pulse_sequence()`.
