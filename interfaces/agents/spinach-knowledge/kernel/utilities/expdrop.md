# kernel/utilities/expdrop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/expdrop.m`
- Signature: `drop=expdrop(from,to,duration,npoints,drop_rate)`
- Total lines: 67

## Purpose

Exponential drop function. Produces an exponential fall-off from a specified value to a specified value with the specified rate and the number of points. Syntax: drop=expdrop(from,to,duration,npoints,drop_rate)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(from,to,duration,npoints,drop_rate)`.
- Lines 33-34: Get the exponential drop parameters; implemented by `B=(from-to)/(1-exp(-drop_rate*duration)); A=from-B`.
- Lines 36-37: Compute the drop; implemented by `drop=A+B*exp(-drop_rate*linspace(0,duration,npoints))`.

### Key state/data transformations

- Lines 34: computes `B` using `B=(from-to)/(1-exp(-drop_rate*duration)); A=from-B`.
- Lines 37: computes `drop` using `drop=A+B*exp(-drop_rate*linspace(0,duration,npoints))`.

### Local helper functions

- Line 42: `grumble()` — `function grumble(from,to,duration,npoints,drop_rate)`.
  - Representative operation: `if (~isnumeric(npoints))||(~isreal(npoints))|| (numel(npoints)~=1)||(npoints<1)||(mod(npoints,1)~=0)`.
  - Representative operation: `(numel(npoints)~=1)||(npoints<1)||(mod(npoints,1)~=0)`.

## Parameters / inputs

- from -the value to drop from
- to -the value to drop to
- duration -drop duration, seconds
- npoints -the number of discretisation points
- in the drop
- drop_rate -exponential drop rate, Hz

## Outputs

- drop -a row vector with the fall-off

## Implementation structure

- Exponential drop function. Produces an exponential fall-off
- from a specified value to a specified value with the specified
- rate and the number of points. Syntax:
- drop=expdrop(from,to,duration,npoints,drop_rate)
- from -the value to drop from
- to -the value to drop to
- duration -drop duration, seconds
- npoints -the number of discretisation points
- in the drop
- drop_rate -exponential drop rate, Hz
- drop -a row vector with the fall-off
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
