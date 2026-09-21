# kernel/contexts/doublerot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/doublerot.m`
- Signature: `[answer,sph_grid]=doublerot(spin_system,pulse_sequence,...`
- Total lines: 551

## Purpose

Double angle spinning context. In Liouville space, this wrapper builds the Fokker-Planck evolution generator and passes it on to the pulse se- quence function, which should be supplied as a handle. In Hilbert space, this wrapper builds the stack of spin Hamiltonians, one for each pair of rotor phases on the two-rotor phase grid, and hands that stack to the pulse sequence. Syntax: [answer,sph_grid]=doublerot(spin_syst

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- answer -the poweder average or a cell array ofwhatever it is
- that the pulse sequence returns
- sph_grid -spherical grid used ithe calculation
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.
- Note: the state projector assumes a powder --single crystal DOR is not
- currently supported.
- Note: the function supports parallel processing via Matlab's Distri-
- buted Computing Toolbox -different system orientations are eva-
- luated on different labs.

## Implementation structure

- Double angle spinning context. In Liouville space, this wrapper builds
- the Fokker-Planck evolution generator and passes it on to the pulse se-
- quence function, which should be supplied as a handle. In Hilbert space,
- this wrapper builds the stack of spin Hamiltonians, one for each pair of
- rotor phases on the two-rotor phase grid, and hands that stack to the
- pulse sequence. Syntax:
- [answer,sph_grid]=doublerot(spin_system,pulse_sequence,...
- parameters,assumptions)
- where pulse sequence is a function handle to one of the pulse sequences
- located in the experiments directory, assumptions is a string that would
- be passed to assume.m when the Hamiltonian is built and parameters is a
- structure with the following subfields:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `frqoffset()`, `ismember()`, `isfield()`, `equilibrium()`, `num2str()`, `fourdif()`, `speye()`, `cart2sph()`, `carrier()`, `relaxation()`.
