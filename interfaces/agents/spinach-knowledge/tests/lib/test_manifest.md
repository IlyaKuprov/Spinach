# tests/lib/test_manifest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_manifest.m`
- Signature: `manifest=test_manifest()`
- Total lines: 130

## Purpose

Returns Spinach regression test metadata. Syntax: manifest=test_manifest()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `test_entry()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-118: Hand-written physically motivated tests; implemented by `manifest=[ test_entry('kernel/pauli_spin_half_algebra','Spin-half angular momentum algebra','test_pauli_spin_half_algebra'), test_entry('kernel/pauli_spin_one_algebra','…`.

### Key state/data transformations

- Lines 14-118: computes `manifest` using `manifest=[ test_entry('kernel/pauli_spin_half_algebra','Spin-half angular momentum algebra','test_pauli_spin_half_algebra'), test_entry('kernel/pauli_spin_one_algebra','…`.

### Local helper functions

- Line 123: `test_entry()` — `function entry=test_entry(id,name,function_name)`. Build one manifest entry
  - Representative operation: `entry.id=id`.
  - Representative operation: `entry.name=name`.

## Outputs

- manifest -structure array with test identifiers and functions

## Implementation structure

- Returns Spinach regression test metadata. Syntax:
- manifest=test_manifest()
- manifest -structure array with test identifiers and functions
- Hand-written physically motivated tests
- Build one manifest entry

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `test_entry()`.
