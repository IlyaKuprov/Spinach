# kernel/contexts/floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/floquet.m`
- Signature: `[answer,sph_grid]=floquet(spin_system,pulse_sequence,...`
- Total lines: 433

## Purpose

Floquet magic angle spinning context. Generates a Liouvillian super- operator and passes it on to the pulse sequence function, which sho- uld be supplied as a handle. Syntax: [answer,sph_grid]=floquet(spin_system,pulse_sequence,... parameters,assumptions) where pulse sequence is a function handle to one of the pulse sequences located in the experiments directory, assumptions is a string that would be passed to assume

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- answer -the poweder average or a cell array ofwhatever it is
- that the pulse sequence returns
- sph_grid -spherical grid used ithe calculation
- Note: the choice of the rank depends on the spinning rate (the slower
- the spinning, the greater ranks are required). The rank is appro-
- ximately equal to the number of spinning sidebands.
- Note: the state projector assumes a powder --single crystal MAS is not
- currently supported.
- Note: perturbative corrections to the rotating frame transformation are
- not supported -use singlerot.m instead.
- Note: the spinning sense matches singlerot.m -the same parameters.rate
- produces the same powder result in both contexts.
- Note: the function supports parallel processing via Matlab's Distri-
- buted Computing Toolbox -different system orientations are eva-
- luated on different labs.

## Implementation structure

- Floquet magic angle spinning context. Generates a Liouvillian super-
- operator and passes it on to the pulse sequence function, which sho-
- uld be supplied as a handle. Syntax:
- [answer,sph_grid]=floquet(spin_system,pulse_sequence,...
- parameters,assumptions)
- where pulse sequence is a function handle to one of the pulse sequences
- located in the experiments directory, assumptions is a string that would
- be passed to assume.m when the Hamiltonian is built and parameters is a
- structure with the following subfields:
- parameters.rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `any()`, `cellfun()`, `num2str()`, `frqoffset()`, `ismember()`, `isfield()`, `equilibrium()`, `relaxation()`, `kinetics()`, `load()`.
