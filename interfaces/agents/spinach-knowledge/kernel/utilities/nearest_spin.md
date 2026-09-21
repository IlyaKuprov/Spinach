# kernel/utilities/nearest_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/nearest_spin.m`
- Signature: `[k,d]=nearest_spin(spin_system,n)`
- Total lines: 64

## Purpose

Returns the index of the nearest spin to the one speci- fied. Only spins for which Cartesian coordinates are available are considered. Syntax: k=nearest_spin(spin_system,n)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- n -index of the spin in question

## Outputs

- k -index of the nearest spin
- d -distance to the nearest spin, Angstrom

## Implementation structure

- Returns the index of the nearest spin to the one speci-
- fied. Only spins for which Cartesian coordinates are
- available are considered. Syntax:
- k=nearest_spin(spin_system,n)
- n -index of the spin in question
- k -index of the nearest spin
- d -distance to the nearest spin, Angstrom
- Check consistency
- Starting point
- Find the nearest spin
- Catch pathological cases
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
