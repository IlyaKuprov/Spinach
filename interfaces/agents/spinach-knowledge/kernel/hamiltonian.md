# kernel/hamiltonian.m

- Signature: `[I,Q]=hamiltonian(spin_system,operator_type)`

## Purpose

Hamiltonian operator or superoperator and its rotational decomposi- tion. Descriptor and operator generation are parallelised. Syntax: [I,Q]=hamiltonian(spin_system,operator_type)

## Physical / mathematical content

- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- With `ham_cache` enabled, the cache identity includes the giant-spin coefficients and retention strengths, so full and Zeeman-only requests do not reuse each other's anisotropic Hamiltonian.

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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

## Header notes

Liouville requests support left, right, commutation, and anticommutation superoperators; Hilbert requests ignore operator_type. The invariant component and anisotropic components are returned separately, and orientation supplies the latter at a chosen geometry.
