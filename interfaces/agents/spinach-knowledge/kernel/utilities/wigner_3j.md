# kernel/utilities/wigner_3j.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/wigner_3j.m`
- Signature: `w=wigner_3j(j1,m1,j2,m2,j3,m3)`
- Total lines: 67

## Purpose

Calculates Wigner 3j-symbols. Syntax: w=wigner_3j(j1,m1,j2,m2,j3,m3) If physically inadmissible indices are supplied, a zero is returned. Order of elements: /j1 j2 j3\ \m1 m2 m3/

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(j1,m1,j2,m2,j3,m3)`.
- Lines 30-31: Call Clebsch-Gordan coefficients; implemented by `w=(((-1)^(-m3+j1+j2))/sqrt(2*j3+1))*clebsch_gordan(j3,-m3,j1,m1,j2,m2)`.

### Key state/data transformations

- Lines 31: computes `w` using `w=(((-1)^(-m3+j1+j2))/sqrt(2*j3+1))*clebsch_gordan(j3,-m3,j1,m1,j2,m2)`.

### Local helper functions

- Line 36: `grumble()` — `function grumble(j1,m1,j2,m2,j3,m3)`.
  - Representative operation: `if (~isnumeric(j1))||(~isnumeric(j2))||(~isnumeric(j3))|| (~isnumeric(m1))||(~isnumeric(m2))||(~isnumeric(m3))`.
  - Representative operation: `(~isnumeric(m1))||(~isnumeric(m2))||(~isnumeric(m3))`.

## Parameters / inputs

- j1-j3 -integers arranged in the order shown above
- m1-m3 -integers arranged in the order shown above

## Outputs

- w -the resulting 3j-symbol

## Implementation structure

- Calculates Wigner 3j-symbols. Syntax:
- w=wigner_3j(j1,m1,j2,m2,j3,m3)
- If physically inadmissible indices are supplied, a zero is
- returned. Order of elements:
- /j1 j2 j3\
- \m1 m2 m3/
- j1-j3 -integers arranged in the order shown above
- m1-m3 -integers arranged in the order shown above
- w -the resulting 3j-symbol
- Check consistency
- Call Clebsch-Gordan coefficients
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `clebsch_gordan()`.
