# kernel/conventions/transforms/g2freq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/g2freq.m`
- Signature: `f=g2freq(g,B)`
- Total lines: 47

## Purpose

Converts g-tensor units into electron Zeeman frequency units. Syntax: f=g2freq(g,B)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(g,B)`.
- Lines 25-26: Get the free electron carrier frequency; implemented by `omega=B*spin('E')/(2*pi)`.
- Lines 28-29: Scale the frequency; implemented by `f=g.*omega/2.0023193043622`.

### Key state/data transformations

- Lines 26: computes `omega` using `omega=B*spin('E')/(2*pi)`.
- Lines 29: computes `f` using `f=g.*omega/2.0023193043622`.

### Local helper functions

- Line 34: `grumble()` — `function grumble(g,B)`. "Authors are listed in order of degree of belief in the central thesis."
  - Representative operation: `if (~isnumeric(g))||(~isreal(g))||any(~isfinite(g),'all')`.
  - Representative operation: `error('g must be a finite real numeric array.')`.

## Parameters / inputs

- g -g-values, scalar or array
- B -magnetic field in Tesla

## Outputs

- f -frequency in Hz

## Implementation structure

- Converts g-tensor units into electron Zeeman frequency
- units. Syntax:
- f=g2freq(g,B)
- g - g-values, scalar or array
- B - magnetic field in Tesla
- f - frequency in Hz
- Check consistency
- Get the free electron carrier frequency
- Scale the frequency
- Consistency enforcement
- "Authors are listed in order of degree of belief in
- the central thesis."

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `any()`, `isscalar()`.
