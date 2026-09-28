# kernel/hamiltonian.m

- Signature: `[I,Q]=hamiltonian(spin_system,operator_type)`

## Purpose

Hamiltonian operator or superoperator and its rotational decomposition. Descriptor and operator generation are parallelised. Syntax: [I,Q]=hamiltonian(spin_system,operator_type)

## Physical / mathematical content

- Quadrupolar interactions are represented by second-rank anisotropic components.

## Numerical / algorithmic content

- With `ham_cache` enabled, the cache hash includes the Hamiltonian descriptor, operator type, isotope and basis hashes, and giant-spin data; when modes are present, it also includes mode data and base frequencies.

- Descriptor and operator generation are parallelised; sparse-matrix assembly includes measures to reduce memory use.

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
