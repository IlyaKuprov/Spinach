# kernel/contexts/powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/powder.m`
- Signature: `[answer,sph_grid]=powder(spin_system,pulse_sequence,...`
- Total lines: 441

## Purpose

Static powder interface to pulse sequences. Generates a Liouvillian superoperator, the initial state, the coil state, then passes them on to the pulse sequence function. Syntax: [answer,sph_grid]=powder(spin_system,pulse_sequence,... parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `defaults()`, `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spins that the
- pulse sequence works on, in the order
- of channels, e.g. {'1H','13C'}
- parameters.offset -a cell array giving transmitter offsets
- in Hz on each of the spins listed in
- parameters.spins
- parameters.grid -name of the spherical averaging grid
- file (see the grids directory in the
- kernel).
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.needs -a cell array of strings specifying ad-
- ditional information required by the
- sequence:
- 'zeeman_op' -Zeeman part of the Hami-
- ltonian in the laboratory frame, to be
- placed into parameters.hzeeman and sent
- to the pulse sequence
- 'iso_eq' -thermal equilibrium is com-
- computed using the isotropic part of
- the Hamiltonian, and sent to the pulse
- sequence via parameters.rho0
- 'aniso_eq' -thermal equilibrium is re-
- computed using the full anisotropic Ha-
- miltonian at each orientation, and sent
- to pulse sequence via parameters.rho0
- parameters.rho0 -initial state; may be a function handle
- that depends on the three Euler angles
- in ZYZ active convention
- parameters.serial -if set to true, disables automatic pa-
- rallelisation
- parameters.sum_up -if set to false, causes the pulse sequ-
- ence output at each orientation to be
- returned instead of the powder average
- parameters.* -additional subfields may be required by
- the pulse sequence -check its documen-
- tation page
- assumptions -context-specific assumptions ('nmr', 'epr',
- 'labframe', etc.) -see the pulse sequence
- header for information on this setting.

## Outputs

- answer -powder average of whatever it is that the pulse
- sequence returns; if parameters.sum_up is set to
- false, a cell array of outputs at each orienta-
- tion is returned
- sph_grid -powder averaging grid data structure with three
- Euler angles and weights for each point
- Note: THIS IS FOR STATIC POWDERS -use singlerot for MAS simulations.
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.
- Note: the function supports parallel processing via Matlab's Distri-
- buted Computing Toolbox -different system orientations are eva-
- luated on different labs.

## Implementation structure

- Static powder interface to pulse sequences. Generates a Liouvillian
- superoperator, the initial state, the coil state, then passes them
- on to the pulse sequence function. Syntax:
- [answer,sph_grid]=powder(spin_system,pulse_sequence,...
- parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spins that the
- pulse sequence works on, in the order
- of channels, e.g. {'1H','13C'}
- parameters.offset -a cell array giving transmitter offsets

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `ismember()`, `hamiltonian()`, `assume()`, `equilibrium()`, `kinetics()`, `frqoffset()`, `carrier()`, `load()`, `isfield()`, `num2str()`, `afterEach()`, `ticBytes()`.
