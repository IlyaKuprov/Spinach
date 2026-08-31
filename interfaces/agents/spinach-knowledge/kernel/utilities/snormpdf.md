# kernel/utilities/snormpdf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/snormpdf.m`
- Signature: `p=snormpdf(x,mu,sigma,alpha)`
- Total lines: 55

## Purpose

Azzalini's skew normal distribution. Syntax: p=snormpdf(x,mu,sigma,alpha)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(x,mu,sigma,alpha)`.

### Key state/data transformations

- Lines 30: computes `p` using `p=2*normpdf(x,mu,sigma).*normcdf(alpha*x,alpha*mu,sigma)`.

### Local helper functions

- Line 35: `grumble()` — `function grumble(x,mu,sigma,alpha)`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))`.
  - Representative operation: `error('x must be a real numeric array.')`.

## Parameters / inputs

- x -an array of real numbers
- mu -expectation value of the normal distribution
- sigma -standard deviation of the normal distribution
- alpha -skew factor, a real number

## Outputs

- p -an array of probability densities,
- same shape as x

## Implementation structure

- Azzalini's skew normal distribution. Syntax:
- p=snormpdf(x,mu,sigma,alpha)
- x -an array of real numbers
- mu -expectation value of the normal distribution
- sigma -standard deviation of the normal distribution
- alpha -skew factor, a real number
- p -an array of probability densities,
- same shape as x
- Check consistency
- Equation 2 in http://www.jstor.org/stable/4615982
- Consistency enforcement
- The smallest minority on earth is the individual. Those who deny

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `normpdf()`, `normcdf()`, `isscalar()`.
