# kernel/utilities/sphten2zeeman.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sphten2zeeman.m`
- Signature: `P=sphten2zeeman(spin_system)`
- Total lines: 77

## Purpose

Returns a matrix that converts state vectors written in the spherical tensor basis set used by Spinach into state vectors written in the Zeeman basis set in Liouville space. Syntax: P=sphten2zeeman(spin_system)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -main Spinach data structure using
- sphten-liouv formalism and inclu-
- ding basis set information

## Outputs

- P -projector matrix that is to be used in the fol-
- lowing way:
- rho_zeeman=P*rho_sphten
- Note: the projector need not be square and may be huge.

## Implementation structure

- Returns a matrix that converts state vectors written in the
- spherical tensor basis set used by Spinach into state vectors
- written in the Zeeman basis set in Liouville space. Syntax:
- P=sphten2zeeman(spin_system)
- spin_system -main Spinach data structure using
- sphten-liouv formalism and inclu-
- ding basis set information
- P -projector matrix that is to be used in the fol-
- lowing way:
- rho_zeeman=P*rho_sphten
- Note: the projector need not be square and may be huge.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `irr_sph_ten()`, `rho()`, `strcmp()`.
