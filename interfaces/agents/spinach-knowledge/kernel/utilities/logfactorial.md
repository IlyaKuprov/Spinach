# kernel/utilities/logfactorial.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/logfactorial.m`
- Signature: `lf=logfactorial(n)`
- Total lines: 44

## Purpose

Logarithm of the factorial function. Avoids complications with factorials of large numbers overflowing 64-bit numbers. Syntax: lf=logfactorial(n)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(n)`.
- Lines 27-28: Use the built-in log(gamma(n)) function; implemented by `lf=gammaln(n+1)`.

### Key state/data transformations

- Lines 28: computes `lf` using `lf=gammaln(n+1)`.

### Local helper functions

- Line 33: `grumble()` — `function grumble(n)`. For ten years or so, my name was "that jerk". But that was a promotion. Before, I was "Who's he?"
  - Representative operation: `if (~isnumeric(n))||(~isreal(n))|| any(n(:)<0)||any(mod(n(:),1)~=0)`.
  - Representative operation: `any(n(:)<0)||any(mod(n(:),1)~=0)`.

## Parameters / inputs

- n -non-negative integer number

## Outputs

- lf -logarithm of the factorial of n
- Notes: double precision overflow is a persistent problem with
- Clebsch-Gordan coefficients and other objects that in-
- volve factorials.

## Implementation structure

- Logarithm of the factorial function. Avoids complications with
- factorials of large numbers overflowing 64-bit numbers. Syntax:
- lf=logfactorial(n)
- n -non-negative integer number
- lf -logarithm of the factorial of n
- Clebsch-Gordan coefficients and other objects that in-
- volve factorials.
- Check consistency
- Use the built-in log(gamma(n)) function
- Consistency enforcement
- For ten years or so, my name was "that jerk". But
- that was a promotion. Before, I was "Who's he?"

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gammaln()`, `any()`.
