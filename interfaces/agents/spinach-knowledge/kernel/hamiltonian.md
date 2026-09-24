# kernel/hamiltonian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/hamiltonian.m`
- Signature: `[I,Q]=hamiltonian(spin_system,operator_type)`
- Total lines: 1478

## Purpose

Hamiltonian operator or superoperator and its rotational decomposi- tion. Descriptor and operator generation are parallelised. Syntax: [I,Q]=hamiltonian(spin_system,operator_type)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- With `ham_cache` enabled, the cache identity includes the giant-spin coefficients and retention strengths, so full and Zeeman-only requests do not reuse each other's anisotropic Hamiltonian.

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `mode_quads()`, `spin_facts()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- in Liouville space, operator_type can be set to
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- in Hilbert space this parameter is ignored.

## Outputs

- I -rotationally invariant part of the Hamiltonian
- Q -irreducible components of the anisotropic part,
- use orientation.m to get the full Hamiltonian
- at each specific orientation
- Note: the code has a few rather eccentric blocks that bring the
- memory footprint to the absolute minimum and work around
- the sparse matrix addition efficiency problem.
- Note: bosonic mode terms declared in inter.modes are added to the
- invariant part I at their input orientation; the Q part and
- the orientation machinery refer to the spin subsystem only,
- for which rotations are well-defined.

## Implementation structure

- Hamiltonian operator or superoperator and its rotational decomposi-
- tion. Descriptor and operator generation are parallelised. Syntax:
- [I,Q]=hamiltonian(spin_system,operator_type)
- in Liouville space, operator_type can be set to
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- in Hilbert space this parameter is ignored.
- I -rotationally invariant part of the Hamiltonian
- Q -irreducible components of the anisotropic part,
- use orientation.m to get the full Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `report()`, `opL()`, `opS()`, `num2str()`, `isotropic()`, `mat2sphten()`, `lam_zeeman()`, `phi_zeeman()`, `phi_quad()`, `table()`, `clear()`, `cellfun()`, `pair_list()`, `lam_coupling()`.
