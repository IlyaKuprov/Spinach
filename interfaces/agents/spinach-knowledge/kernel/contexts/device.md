# kernel/contexts/device.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/device.m`
- Signature: `answer=device(spin_system,pulse_sequence,parameters,assumptions)`
- Total lines: 276

## Purpose

Spin-boson device interface to pulse sequences. Generates the evolution generators for a device containing spins and bosonic modes at a fixed orientation of the spin subsystem, and passes them to the pulse sequence function, which should be supplied as a handle. Syntax: answer=device(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spin species that
- the pulse sequence works on, in the order
- of channels, e.g. {'E'}; may be omitted
- when no spin channels are needed
- parameters.offset -transmitter offsets in Hz on each of the
- spin species listed in parameters.spins
- parameters.mode_offset -detuning offsets in Hz, one for each
- bosonic mode in the order of declaration;
- the transmitter sign convention of para-
- meters.offset applies: each offset enters
- the Hamiltonian as minus the offset times
- the number operator of its mode, so that a
- positive offset lowers the mode frequency
- parameters.decouple -a cell array of spin species to be wiped
- from the evolution generators and from the
- initial state, e.g. {'1H'}; the default is
- an empty cell array, meaning no decoupling
- parameters.orientation -Euler angles (ZYZ active convention, ra-
- dians) giving the orientation of the spin
- subsystem; bosonic terms are not affected;
- the default is [0 0 0]
- parameters.rframes -numerical rotating frame specification for
- spin species, e.g. {{'E',2}}, as described
- in the header of rotframe.m
- parameters.needs -a cell array of strings specifying additional
- information required by the sequence:
- 'rho_eq' -thermal equilibrium state at the
- system temperature, with the Bose-Einstein
- populations of the bosonic modes included;
- this is placed into parameters.rho0
- parameters.* -additional subfields may be required by
- the pulse sequence -check its docs
- assumptions -'labframe', 'cavity', or 'spin-phonon';
- see the header of assume.m for details

## Outputs

- answer -whatever it is that the pulse sequence returns.
- Note: the system must contain at least one bosonic mode; pure spin
- systems belong with liquid.m, crystal.m, and powder.m contexts.
- Note: dissipative bosonic modes require a Liouville space formalism;
- coherent simulations may also use zeeman-hilb.

## Implementation structure

- Spin-boson device interface to pulse sequences. Generates the evolution
- generators for a device containing spins and bosonic modes at a fixed
- orientation of the spin subsystem, and passes them to the pulse sequence
- function, which should be supplied as a handle. Syntax:
- answer=device(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spin species that
- the pulse sequence works on, in the order
- of channels, e.g. {'E'}; may be omitted
- when no spin channels are needed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `assume()`, `hamiltonian()`, `orientation()`, `relaxation()`, `kinetics()`, `ismember()`, `report()`, `equilibrium()`, `frqoffset()`, `num2str()`, `mode_list()`, `operator()`, `carrier()`.
