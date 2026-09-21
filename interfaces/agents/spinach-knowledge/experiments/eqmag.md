# experiments/eqmag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/eqmag.m`
- Signature: `magn=eqmag(spin_system,parameters)`
- Total lines: 120

## Purpose

Computes the molar magnetization vector at the thermal equilibrium at the temperature specified in inter.temperature and magnetic field spe- cified in sys.magnet (assumed to be along the Z-axis), averaged over system orientations using the spherical grid specified. Syntax: magn=eqmag(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.grid -spherical grid for averaging

## Outputs

- magn -molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
- Note: the use of bas.formalism='zeeman-hilb' is required.
- Note: Spinach uses NMR convention for the exchange coupling: exchange
- interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
- where J is in Hz.

## Implementation structure

- Computes the molar magnetization vector at the thermal equilibrium at
- the temperature specified in inter.temperature and magnetic field spe-
- cified in sys.magnet (assumed to be along the Z-axis), averaged over
- system orientations using the spherical grid specified. Syntax:
- magn=eqmag(spin_system,parameters)
- parameters.grid -spherical grid for averaging
- magn -molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
- Note: the use of bas.formalism='zeeman-hilb' is required.
- Note: Spinach uses NMR convention for the exchange coupling: exchange
- interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
- where J is in Hz.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gtensorof()`, `operator()`, `hamiltonian()`, `assume()`, `load()`, `euler2dcm()`, `alphas()`, `betas()`, `gammas()`, `g_rot()`, `equilibrium()`, `strcmp()`, `isfield()`, `ischar()`.
