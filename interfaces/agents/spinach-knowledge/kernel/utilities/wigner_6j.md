# kernel/utilities/wigner_6j.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/wigner_6j.m`
- Signature: `w=wigner_6j(j1,j2,j3,j4,j5,j6)`
- Total lines: 90

## Purpose

Wigner 6j-symbols. Syntax: w=wigner_6j(j1,j2,j3,j4,j5,j6) If physically inadmissible indices are supplied, a zero is returned. Order of elements: / j1 j2 j3 \ \ j4 j5 j6 /

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(j1,j2,j3,j4,j5,j6)`.
- Lines 28-29: Start from zero; implemented by `w=0`.
- Lines 31-32: Use the definition; implemented by `for m1=-j1:j1`.
- Lines 39-40: Get the power index for -1; implemented by `power_of_minus=(j1-m1)+(j2-m2)+(j3-m3)+(j4-m4)+(j5-m5)+(j6-m6)`.
- Lines 42-43: Get the screen; implemented by `screen=(m1+m2-m3==0)&(-m1+m5+m6==0)&(m4-m5+m3==0)&(-m4-m2-m6==0)`.
- Lines 45-46: Add to the total; implemented by `if screen`.

### Control flow inferred from the code

- Line 32: `for` loop over `m1=-j1:j1`.
- Line 33: `for` loop over `m2=-j2:j2`.
- Line 34: `for` loop over `m3=-j3:j3`.
- Line 35: `for` loop over `m4=-j4:j4`.
- Line 36: `for` loop over `m5=-j5:j5`.
- Line 37: `for` loop over `m6=-j6:j6`.
- Line 46: conditional branch on `screen`.

### Key state/data transformations

- Lines 29: computes `w` using `w=0`.
- Lines 40: computes `power_of_minus` using `power_of_minus=(j1-m1)+(j2-m2)+(j3-m3)+(j4-m4)+(j5-m5)+(j6-m6)`.

### Local helper functions

- Line 63: `grumble()` — `function grumble(j1,j2,j3,j4,j5,j6)`.
  - Representative operation: `if (~isnumeric(j1))||(~isnumeric(j3))||(~isnumeric(j5))|| (~isnumeric(j2))||(~isnumeric(j4))||(~isnumeric(j6))`.
  - Representative operation: `(~isnumeric(j2))||(~isnumeric(j4))||(~isnumeric(j6))`.

## Parameters / inputs

- j1-j6 -integers arranged in the order shown above

## Outputs

- w -the resulting 6j-symbol

## Implementation structure

- Wigner 6j-symbols. Syntax:
- w=wigner_6j(j1,j2,j3,j4,j5,j6)
- If physically inadmissible indices are supplied, a zero is
- returned. Order of elements:
- / j1 j2 j3 \
- \ j4 j5 j6 /
- j1-j6 -integers arranged in the order shown above
- w -the resulting 6j-symbol
- Check consistency
- Start from zero
- Use the definition
- Get the power index for -1

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `wigner_3j()`.
