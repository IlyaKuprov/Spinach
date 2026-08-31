# kernel/conventions/transforms/hartree2joule.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/hartree2joule.m`
- Signature: `energy=hartree2joule(energy)`
- Total lines: 43

## Purpose

Converts Hartree energy units into J/mol. A Hartree is twice the ground state ionisation energy of the hydrogen atom. Syntax: energy=hartree2joule(energy)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(energy)`.
- Lines 26-27: Perform the conversion; implemented by `energy=2625499.62*energy`.

### Key state/data transformations

- Lines 27: computes `energy` using `energy=2625499.62*energy`.

### Local helper functions

- Line 32: `grumble()` — `function grumble(energy)`. "We only need to be lucky once. You need to be lucky every time."
  - Representative operation: `if (~isnumeric(energy))||(~isreal(energy))`.
  - Representative operation: `error('the argument must be an array of real numbers.')`.

## Parameters / inputs

- energy -a numerical array of energies in
- Hartree units

## Outputs

- energy -a numerical array of energies in
- Joules

## Implementation structure

- Converts Hartree energy units into J/mol. A Hartree is twice the
- ground state ionisation energy of the hydrogen atom. Syntax:
- energy=hartree2joule(energy)
- energy -a numerical array of energies in
- Hartree units
- Joules
- Check consistency
- Perform the conversion
- Consistency enforcement
- "We only need to be lucky once. You need to be
- lucky every time."
- The IRA to Margaret Thatcher, after

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
