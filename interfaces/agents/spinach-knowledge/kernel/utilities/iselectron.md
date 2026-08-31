# kernel/utilities/iselectron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/iselectron.m`
- Signature: `verdict=iselectron(spin_spec)`
- Total lines: 42

## Purpose

Returns true if the particle is an electron. Syntax: verdict=iselectron(spin_spec)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `spin()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Check consistency; implemented by `grumble(spin_spec)`.
- Lines 22-23: A simple matching check; implemented by `if ismember(spin_spec(1),{'E'})`.

### Control flow inferred from the code

- Line 23: conditional branch on `ismember(spin_spec(1),{'E'})`.

### Key state/data transformations

- Lines 24: computes `verdict` using `verdict=true()`.

### Local helper functions

- Line 32: `grumble()` — `function grumble(spin_spec)`. Лучше умереть героем чем жить пидорасом. Евгений Пригожин
  - Representative operation: `if ~ischar(spin_spec)`.
  - Representative operation: `error('spin_spec must be a character string.')`.

## Parameters / inputs

- spin_spec -a Spinach particle specification

## Outputs

- verdict -true for an electron, false otherwise

## Implementation structure

- Returns true if the particle is an electron. Syntax:
- verdict=iselectron(spin_spec)
- spin_spec -a Spinach particle specification
- verdict -true for an electron, false otherwise
- Check consistency
- A simple matching check
- Consistency enforcement
- Лучше умереть героем чем жить пидорасом.
- Евгений Пригожин

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `spin_spec()`, `true()`, `false()`, `ischar()`, `spin()`.
